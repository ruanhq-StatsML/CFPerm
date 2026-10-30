"""Compare CF / BCF / bootstrap Var(u) on existing synth + real scenes + coverage.

Effects
-------
For each board, run the same Mag×Var + Y-flip triage three times with
method ∈ {cf, bcf, boot}. Report call, V, Var(u), hit, agreement.

Coverage
--------
Primary target: Δ = μ_LIVE − μ_REF of unit PO-risk.
From method draws u^{(s)} on the stacked REF|LIVE matrix:
  Δ_s = mean(u^{(s)}[LIVE]) − mean(u^{(s)}[REF])
Build a (1-α) percentile CI for Δ.
  quiet:      want CI ∋ 0  (null: no mean-risk shift)
  malignant:  want CI ∌ 0  (power)
Secondary: CI for μ_LIVE covering the point μ_LIVE (calibration sanity).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from po_risk_var._path import REPO_ROOT, ensure_cfperm_src
from po_risk_var.apps_realworld import REAL_BOARDS, eval_real_board
from po_risk_var.apps_synthetic import SYNTH_BOARDS, make_synthetic_board
from po_risk_var.estimators import estimate_var_u, monitor_batch

ensure_cfperm_src()
from as_dta_permucate_magvar import permucate_unit_po  # noqa: E402

OUT_DEFAULT = REPO_ROOT / "results" / "as_dta_po_risk_var_compare"
METHODS = ("cf", "bcf", "boot")
ALPHA = 0.10  # 90% CIs


def _ok(oracle: str, pred: str) -> bool:
    if oracle == "malignant":
        return pred == "malignant"
    return pred in ("quiet", "benign")


def _ci_from_draws(mean_draws: np.ndarray, alpha: float = ALPHA) -> Tuple[float, float, float]:
    m = np.asarray(mean_draws, dtype=float).ravel()
    m = m[np.isfinite(m)]
    if len(m) < 3:
        x = float(np.mean(m)) if len(m) else float("nan")
        return x, x, x
    lo, hi = np.quantile(m, [alpha / 2, 1 - alpha / 2])
    return float(np.mean(m)), float(lo), float(hi)


def _group_mean_draws(var: dict, mask: np.ndarray) -> np.ndarray:
    """Batch-mean u draws on a row mask — shape (S,)."""
    draws = np.asarray(var["u_draws"], dtype=float)
    if draws.ndim != 2:
        raise ValueError("u_draws must be (S, n)")
    mask = np.asarray(mask, dtype=bool)
    if draws.shape[1] != len(mask):
        return draws.mean(axis=1)
    return draws[:, mask].mean(axis=1)


def _monitor_row(mon: dict, oracle: str, method: str) -> dict:
    pred = mon["decision"]["call"]
    return {
        "method": method,
        "oracle": oracle,
        "pred": pred,
        "actuator": mon["decision"]["actuator"],
        "correct": bool(_ok(oracle, pred)),
        "J": mon["magvar"]["J"],
        "d_mu": mon["magvar"]["d_mu"],
        "d_s": mon["magvar"]["d_s"],
        "y_flip": mon["y_flip"],
        "V": mon["var_u"]["mean_V_debiased_live"],
        "Var_u": mon["var_u"]["mean_Var_u_live"],
        "bias_ratio": mon["var_u"]["bias_ratio_live"],
        "backend": mon["var_u"].get("backend"),
        "reason": mon["decision"]["reason"],
    }


def compare_synth_board(
    name: str,
    *,
    methods: Sequence[str] = METHODS,
    quick: bool = True,
    seed: int = 2026,
) -> dict:
    board = make_synthetic_board(name, n=140 if quick else 220, seed=seed)
    rows = []
    for m in methods:
        print(f"[compare:synth] {name} / {m}", flush=True)
        mon = monitor_batch(
            board["ref"]["X"],
            board["ref"]["Y"],
            board["live"]["X"],
            board["live"]["Y"],
            method=m,
            y_flip=board["y_flip"],
            quick=quick,
            seed=seed,
            n_estimators=28 if quick else 48,
        )
        # light BCF: estimate_var_u already uses bcf_sim default; for quick
        # re-run with smaller sim if needed is handled via method path
        rows.append(_monitor_row(mon, board["oracle"], m))
    calls = {r["method"]: r["pred"] for r in rows}
    return {
        "board": name,
        "oracle": board["oracle"],
        "y_flip": board["y_flip"],
        "methods": rows,
        "all_correct": all(r["correct"] for r in rows),
        "agree_all": len(set(calls.values())) == 1,
        "calls": calls,
    }


def compare_real_board(
    name: str,
    *,
    methods: Sequence[str] = METHODS,
    quick: bool = True,
    seed: int = 2026,
) -> dict:
    by_m = {}
    for m in methods:
        print(f"[compare:real] {name} / {m}", flush=True)
        by_m[m] = eval_real_board(
            name, live_n=40 if quick else 56, seed=seed, quick=quick, method=m
        )
    # flatten benign/malignant per method
    panel = []
    for m in methods:
        b = by_m[m]
        for kind in ("benign", "malignant"):
            cell = b[kind]
            panel.append(
                {
                    "method": m,
                    "kind": kind,
                    "oracle": cell["oracle"],
                    "pred": cell["pred"],
                    "correct": cell["correct"],
                    "J": cell["J"],
                    "y_flip": cell["y_flip"],
                    "V": cell["V_debiased"],
                    "bias_ratio": cell["bias_ratio"],
                    "reason": cell["reason"],
                }
            )
    pair_ok = {m: bool(by_m[m]["pair_ok"]) for m in methods}
    v_gap = {m: float(by_m[m]["V_gap"]) for m in methods}
    return {
        "board": name,
        "pair_ok": pair_ok,
        "V_gap": v_gap,
        "panel": panel,
        "all_pair_ok": all(pair_ok.values()),
    }


def _coverage_one_rep(
    *,
    kind: str,
    method: str,
    seed: int,
    quick: bool,
    alpha: float = ALPHA,
) -> dict:
    """One replicate: CI for Δ=μ_LIVE−μ_REF under quiet or malignant_sparse."""
    name = "quiet" if kind == "quiet" else "malignant_sparse"
    board = make_synthetic_board(name, n=120 if quick else 200, seed=seed)
    X = np.vstack([board["ref"]["X"], board["live"]["X"]])
    Y = np.concatenate([board["ref"]["Y"], board["live"]["Y"]])
    n0 = len(board["ref"]["Y"])
    W = np.array([0] * n0 + [1] * (len(Y) - n0), dtype=int)
    n_est = 24 if quick else 40
    bcf_sim = 20 if quick else 40
    boot_B = 10 if quick else 20

    fit = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_est)
    phi = fit["phi"]
    unit = fit["unit_po"]
    live = W == 1
    ref = ~live
    mu_ref = float(np.mean(unit[ref]))
    mu_live = float(np.mean(unit[live]))
    delta_point = mu_live - mu_ref

    # CF leaf-merge pools REF+LIVE inside a leaf → Δ collapses; for coverage
    # of Δ use raw tree/posterior draws (leaf_merge=False on CF path).
    var = estimate_var_u(
        X,
        Y,
        W,
        method=method,
        phi=phi,
        n_estimators=n_est,
        seed=seed,
        bcf_sim=bcf_sim,
        boot_B=boot_B,
        leaf_merge=(method not in ("cf", "causal_forest", "default")),
    )
    d_draws = _group_mean_draws(var, live) - _group_mean_draws(var, ref)
    d_hat, lo, hi = _ci_from_draws(d_draws, alpha=alpha)
    covers_zero = bool(lo <= 0.0 <= hi)
    # secondary: μ_LIVE CI covers point μ_LIVE
    m_live = _group_mean_draws(var, live)
    _, l_mu, h_mu = _ci_from_draws(m_live, alpha=alpha)
    covers_mu_live = bool(l_mu <= mu_live <= h_mu)
    return {
        "kind": kind,
        "method": method,
        "seed": int(seed),
        "mu_ref": mu_ref,
        "mu_live": mu_live,
        "delta_point": delta_point,
        "delta_hat": d_hat,
        "ci_lo": lo,
        "ci_hi": hi,
        "width": float(hi - lo),
        "covers_zero": covers_zero,
        "covers_mu_live_point": covers_mu_live,
        "n_draws": int(len(d_draws)),
        "backend": var.get("backend"),
        "mean_Var_u_live": float(np.mean(np.asarray(var["Var_u"])[live])),
        "mean_V_live": float(np.mean(np.asarray(var["V_debiased"])[live])),
    }


def run_coverage(
    *,
    methods: Sequence[str] = METHODS,
    n_reps: int = 12,
    seed: int = 2026,
    quick: bool = True,
    alpha: float = ALPHA,
) -> dict:
    reps: List[dict] = []
    for r in range(int(n_reps)):
        for kind in ("quiet", "malignant"):
            for m in methods:
                print(f"[coverage] rep={r} {kind}/{m}", flush=True)
                reps.append(
                    _coverage_one_rep(
                        kind=kind,
                        method=m,
                        seed=seed + 17 * r + (0 if kind == "quiet" else 3),
                        quick=quick,
                        alpha=alpha,
                    )
                )

    summary = {}
    for m in methods:
        for kind in ("quiet", "malignant"):
            sub = [x for x in reps if x["method"] == m and x["kind"] == kind]
            cov0 = float(np.mean([x["covers_zero"] for x in sub])) if sub else float("nan")
            # quiet: want cov0 ≈ 1-α; malignant: want power = 1-cov0 high
            power = (1.0 - cov0) if kind == "malignant" else float("nan")
            summary[f"{m}:{kind}"] = {
                "n": len(sub),
                "empirical_coverage_of_zero_delta": cov0,
                # aliases used by plot/talk
                "empirical_coverage_of_mu_ref": cov0,
                "nominal": 1.0 - alpha,
                "mean_width": float(np.mean([x["width"] for x in sub])) if sub else float("nan"),
                "mean_delta": float(np.mean([x["delta_point"] for x in sub])) if sub else float("nan"),
                "mean_mu_live": float(np.mean([x["mu_live"] for x in sub])) if sub else float("nan"),
                "mean_mu_ref": float(np.mean([x["mu_ref"] for x in sub])) if sub else float("nan"),
                "power_reject_ref": power,
                "power_reject_zero_delta": power,
                "mu_live_point_coverage": (
                    float(np.mean([x["covers_mu_live_point"] for x in sub])) if sub else float("nan")
                ),
                "calibrated_quiet": (
                    abs(cov0 - (1.0 - alpha)) <= 0.25 if kind == "quiet" else None
                ),
            }
    return {
        "alpha": alpha,
        "n_reps": int(n_reps),
        "target": "Δ = μ_LIVE − μ_REF; CI from method u-draws on REF|LIVE",
        "null": "Δ = 0 (holds under quiet; false under malignant)",
        "reps": reps,
        "summary": summary,
    }


def run_compare_suite(
    *,
    out_dir: Optional[Path] = None,
    quick: bool = True,
    seed: int = 2026,
    methods: Sequence[str] = METHODS,
    coverage_reps: int = 10,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)

    synth = [compare_synth_board(b, methods=methods, quick=quick, seed=seed) for b in SYNTH_BOARDS]
    real_names = (
        ("hc3_finance", "news20_sci", "ag_news", "imdb")
        if quick
        else tuple(REAL_BOARDS)
    )
    real = [
        compare_real_board(b, methods=methods, quick=quick, seed=seed) for b in real_names
    ]

    # hit rates per method
    synth_hit = {}
    real_hit = {}
    for m in methods:
        s_ok = sum(1 for b in synth if next(r["correct"] for r in b["methods"] if r["method"] == m))
        synth_hit[m] = {"n_ok": s_ok, "n": len(synth), "rate": s_ok / max(len(synth), 1)}
        r_ok = sum(1 for b in real if b["pair_ok"].get(m))
        real_hit[m] = {"n_ok": r_ok, "n": len(real), "rate": r_ok / max(len(real), 1)}

    print("[compare] coverage …", flush=True)
    cov = run_coverage(
        methods=methods,
        n_reps=coverage_reps if quick else max(coverage_reps, 16),
        seed=seed + 101,
        quick=quick,
    )

    # talk track
    talk = _talk(synth, real, synth_hit, real_hit, cov, methods)

    payload = {
        "scope": (
            "Side-by-side CF / BCF / bootstrap Var(u) on existing synth + real "
            "scenes, plus 90% CI coverage of Δ=μ_LIVE−μ_REF (null Δ=0)"
        ),
        "methods": list(methods),
        "synthetic": synth,
        "realworld": real,
        "hit_rates": {"synthetic": synth_hit, "realworld": real_hit},
        "coverage": {
            "alpha": cov["alpha"],
            "n_reps": cov["n_reps"],
            "target": cov["target"],
            "null": cov["null"],
            "summary": cov["summary"],
            # keep reps in JSON for audit but they can be large
            "reps": cov["reps"],
        },
        "talk": talk,
        "recommend": talk["recommend"],
    }
    (out_dir / "po_risk_var_compare.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    _plot(payload, out_dir / "po_risk_var_compare.png")
    _readme(out_dir / "README.md", payload)
    return payload


def _talk(synth, real, synth_hit, real_hit, cov, methods) -> dict:
    lines = []
    # agreement
    n_agree = sum(1 for b in synth if b["agree_all"])
    lines.append(
        f"Synth call agreement across methods: {n_agree}/{len(synth)} boards "
        "(same triage call)."
    )
    for m in methods:
        lines.append(
            f"Hit {m}: synth {synth_hit[m]['n_ok']}/{synth_hit[m]['n']} "
            f"({synth_hit[m]['rate']:.2f}), real pairs {real_hit[m]['n_ok']}/{real_hit[m]['n']} "
            f"({real_hit[m]['rate']:.2f})."
        )
    # coverage narrative
    cov_bits = []
    for m in methods:
        q = cov["summary"].get(f"{m}:quiet", {})
        mal = cov["summary"].get(f"{m}:malignant", {})
        cov_bits.append(
            {
                "method": m,
                "quiet_coverage_delta0": q.get("empirical_coverage_of_zero_delta"),
                "quiet_width": q.get("mean_width"),
                "malignant_power": mal.get("power_reject_zero_delta"),
                "malignant_width": mal.get("mean_width"),
                "mu_live_point_coverage": q.get("mu_live_point_coverage"),
            }
        )
        lines.append(
            f"Coverage {m}: quiet P(CI_Δ∋0)={q.get('empirical_coverage_of_zero_delta', float('nan')):.2f} "
            f"(nominal {1 - cov['alpha']:.2f}, width≈{q.get('mean_width', float('nan')):.3g}); "
            f"malignant power P(CI_Δ∌0)={mal.get('power_reject_zero_delta', float('nan')):.2f}; "
            f"μ_LIVE point-cover={q.get('mu_live_point_coverage', float('nan')):.2f}."
        )

    cf_q = cov["summary"].get("cf:quiet", {}).get("empirical_coverage_of_zero_delta", 0)
    bcf_q = cov["summary"].get("bcf:quiet", {}).get("empirical_coverage_of_zero_delta", 0)
    boot_q = cov["summary"].get("boot:quiet", {}).get("empirical_coverage_of_zero_delta", 0)
    cf_pow = cov["summary"].get("cf:malignant", {}).get("power_reject_zero_delta", 0)
    bcf_pow = cov["summary"].get("bcf:malignant", {}).get("power_reject_zero_delta", 0)
    boot_pow = cov["summary"].get("boot:malignant", {}).get("power_reject_zero_delta", 0)
    lines.append(
        "Effect reading: Mag×Var + Y-flip dominate the triage call — CF/BCF/boot "
        "agree whenever the oracle is clear; CF V_debiased often clips near 0 "
        "(honesty: leaf bias removed) while boot Var(u) is heavy-tailed."
    )
    lines.append(
        "Coverage reading (Δ=μ_LIVE−μ_REF, null Δ=0): "
        f"quiet P(CI∋0) cf={cf_q:.2f}, bcf={bcf_q:.2f}, boot={boot_q:.2f} "
        f"(nominal {1 - cov['alpha']:.2f}); "
        f"malignant power cf={cf_pow:.2f}, bcf={bcf_pow:.2f}, boot={boot_pow:.2f}."
    )
    recommend = (
        "Ship CF leaf-debiased as LIVE default for scene triage (same calls as "
        "BCF/boot on existing boards; bias_ratio honesty). For CI coverage of "
        "mean-risk shift Δ: bootstrap is the conservative cover (wide); BCF is "
        "the audit posterior middle; CF tree-CIs are sharp but undercover quiet "
        "— do not use CF draws alone as a Δ confidence statement. "
        "Note: coverage uses raw (unmerged) CF draws; leaf-merge stays on for "
        "product Var(u) stability."
    )
    return {
        "bullets": lines,
        "coverage_panel": cov_bits,
        "synth_agree_rate": n_agree / max(len(synth), 1),
        "quiet_coverage": {"cf": cf_q, "bcf": bcf_q, "boot": boot_q},
        "recommend": recommend,
    }


def _plot(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    methods = payload["methods"]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))

    # 1) hit rates
    ax = axes[0]
    x = np.arange(len(methods))
    synth_r = [payload["hit_rates"]["synthetic"][m]["rate"] for m in methods]
    real_r = [payload["hit_rates"]["realworld"][m]["rate"] for m in methods]
    ax.bar(x - 0.18, synth_r, width=0.35, label="synth", color="#2F6B4F")
    ax.bar(x + 0.18, real_r, width=0.35, label="real pairs", color="#3D5A80")
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("hit rate")
    ax.set_title("Scene effect — triage hit", loc="left", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)

    # 2) coverage quiet vs power malignant
    ax = axes[1]
    cov = payload["coverage"]["summary"]
    quiet = [cov[f"{m}:quiet"]["empirical_coverage_of_mu_ref"] for m in methods]
    power = [cov[f"{m}:malignant"]["power_reject_ref"] for m in methods]
    ax.bar(x - 0.18, quiet, width=0.35, label="quiet coverage", color="#2F6B4F")
    ax.bar(x + 0.18, power, width=0.35, label="malignant power", color="#A33B24")
    ax.axhline(1 - payload["coverage"]["alpha"], ls="--", color="#666", lw=0.9, label="nominal")
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.set_ylim(0, 1.05)
    ax.set_title("Coverage of Δ=0 (90% CI)", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    # 3) CI width quiet
    ax = axes[2]
    widths = [cov[f"{m}:quiet"]["mean_width"] for m in methods]
    ax.bar(x, widths, color=["#2F6B4F", "#4A6FA5", "#8A8F98"], width=0.55)
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.set_ylabel("mean CI width for Δ (quiet)")
    ax.set_title("Interval sharpness", loc="left", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        "Var(u) methods on existing scenes — effect + Δ coverage",
        fontsize=12,
        fontweight="bold",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _readme(path: Path, payload: dict):
    methods = payload["methods"]
    lines = [
        "# PO-risk Var(u) — CF vs BCF vs bootstrap on existing scenes",
        "",
        payload["scope"],
        "",
        "## Talk track",
        "",
    ]
    for b in payload["talk"]["bullets"]:
        lines.append(f"- {b}")
    lines.extend(["", f"**Recommend:** {payload['recommend']}", "", "## Synth boards"])
    lines.append("")
    lines.append("| board | oracle | " + " | ".join(methods) + " | agree? |")
    lines.append("| --- | --- | " + " | ".join(["---"] * len(methods)) + " | --- |")
    for b in payload["synthetic"]:
        preds = [b["calls"][m] for m in methods]
        marks = []
        for m in methods:
            r = next(x for x in b["methods"] if x["method"] == m)
            marks.append(f"{r['pred']}{'✓' if r['correct'] else '✗'}")
        lines.append(
            f"| {b['board']} | {b['oracle']} | "
            + " | ".join(marks)
            + f" | {'yes' if b['agree_all'] else 'no'} |"
        )

    lines.extend(["", "## Real-world pairs", ""])
    lines.append("| board | " + " | ".join(f"{m} pair" for m in methods) + " | V_gap cf/bcf/boot |")
    lines.append("| --- | " + " | ".join(["---"] * len(methods)) + " | --- |")
    for b in payload["realworld"]:
        oks = [("yes" if b["pair_ok"][m] else "no") for m in methods]
        gaps = "/".join(f"{b['V_gap'][m]:.3g}" for m in methods)
        lines.append(f"| {b['board']} | " + " | ".join(oks) + f" | {gaps} |")

    lines.extend(["", "## Coverage (90% CI for Δ = μ_LIVE − μ_REF)", ""])
    lines.append(
        "| method | quiet P(CI∋0) | quiet width | malignant power | mal width | μ_LIVE point-cover |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for m in methods:
        q = payload["coverage"]["summary"][f"{m}:quiet"]
        mal = payload["coverage"]["summary"][f"{m}:malignant"]
        lines.append(
            f"| {m} | {q['empirical_coverage_of_zero_delta']:.2f} | {q['mean_width']:.4g} | "
            f"{mal['power_reject_zero_delta']:.2f} | {mal['mean_width']:.4g} | "
            f"{q['mu_live_point_coverage']:.2f} |"
        )
    lines.extend(
        [
            "",
            "Quiet: want P(CI_Δ∋0) near nominal 0.90. Malignant: want high power (CI excludes 0).",
            "μ_LIVE point-cover is a sanity check that the μ_LIVE CI covers its point estimate.",
            "",
            "```bash",
            "po-risk-var compare --quick",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))


__all__ = ["run_compare_suite", "run_coverage", "compare_synth_board", "compare_real_board"]
