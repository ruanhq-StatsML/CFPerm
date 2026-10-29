"""MVP: Benign OOD vs Malignant OOD — PO mean + variance, not hard-filter/pseudo-Y.

Justify
-------
Both benign synonym/paraphrase and malignant TextFooler-lite can raise
**mean** OOD / PO-risk. They differ in **shape**:

  * Benign  — uniform PO inflation → ΔMean↑, Var/Peak/CVaR quiet, top-share flat
  * Malignant — concentrated PO on attacked minority → ΔVar/Peak/CVaR↑, FN↑

Anti-patterns this MVP visualizes as *failure*:

  * Hard filter  — drop high-PO / high-OOD rows from uplift/monitor pool
                   (selection bias; removes the signal you need)
  * Pseudo-label — map high PO → flip Y (or "washed ⇒ human")
                   (label noise; confuses benign OOD with attack)

Correct path: Magnitude first → Variance/Peak second → actuator
  (conformal for benign, freeze-lexical for malignant).
"""
from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from as_dta_aigc import (
    ROOT,
    SEED,
    adapt_freeze_lexical,
    build_live_batches,
    conformal_calibrate,
    diagnose_batch,
    discriminate_ood_vs_intentional,
    evaluate_detector,
    fit_sklearn_detector,
    load_hc3_binary,
    route_action,
    split_ref_live,
)

OUT_DEFAULT = ROOT / "results" / "as_dta_benign_malignant_mvp"


def _shape_first_route(delta: dict, d_error: float = 0.0) -> str:
    """Magnitude gate → variance / error shape (PO-risk roadmap).

    Absolute Peak alone is noisy under TF-IDF synonym; the durable separators
    are sub-env Std, unit Var, and (with Y) error/FN lift.
    """
    if delta["d_mean"] < 0.015 and float(d_error) < 0.02:
        return "hold"
    malignant = (
        float(d_error) >= 0.12
        or delta["d_sub_std"] >= 0.018
        or delta.get("d_unit_var", 0.0) >= 0.03
        or (delta["d_cvar"] >= 0.18 and float(d_error) >= 0.08)
    )
    if malignant:
        return "freeze_lexical_adapt_syntactic"
    if delta["d_mean"] >= 0.015:
        return "conformal_calibrate"
    return "monitor"


def _shape_verdict(delta: dict, d_error: float = 0.0) -> str:
    r = _shape_first_route(delta, d_error=d_error)
    return {
        "hold": "quiet",
        "conformal_calibrate": "benign_ood",
        "freeze_lexical_adapt_syntactic": "malignant_ood",
        "monitor": "monitor",
    }.get(r, "monitor")


def _po_stats(po: np.ndarray) -> dict:
    po = np.asarray(po, dtype=float)
    if po.size == 0:
        return {
            k: float("nan")
            for k in ("mean", "std", "var", "p90", "p95", "max", "top20_share", "gini")
        }
    order = np.sort(po)[::-1]
    top_n = max(1, int(0.2 * len(po)))
    top_share = float(order[:top_n].sum() / max(order.sum(), 1e-12))
    xs = np.sort(po)
    n = len(xs)
    gini = float((2.0 * np.arange(1, n + 1) - n - 1) @ xs / (n * xs.sum() + 1e-12))
    return {
        "mean": float(np.mean(po)),
        "std": float(np.std(po)),
        "var": float(np.var(po)),
        "p90": float(np.percentile(po, 90)),
        "p95": float(np.percentile(po, 95)),
        "max": float(np.max(po)),
        "top20_share": top_share,
        "gini": gini,
    }


def _subbatch_risk_var(
    texts: Sequence[str],
    y: np.ndarray,
    det,
    *,
    K: int = 40,
    seed: int = SEED,
) -> Tuple[np.ndarray, dict]:
    """Bootstrap sub-environment means → variance of PO-risk across envs."""
    texts = list(texts)
    y = np.asarray(y, dtype=int)
    unit = det.po_proxy(texts, y)
    n = len(texts)
    rng = np.random.default_rng(int(seed))
    win = max(12, n // 5)
    means = []
    for _ in range(int(K)):
        idx = rng.choice(n, size=min(win, n), replace=False)
        means.append(float(np.mean(unit[idx])))
    means = np.asarray(means, dtype=float)
    return means, {
        "sub_mean": float(np.mean(means)),
        "sub_std": float(np.std(means)),
        "sub_var": float(np.var(means)),
        "peak_gap": float(np.max(means) - np.median(means)),
        "cvar95": float(np.mean(means[means >= np.percentile(means, 95)])),
        "unit": _po_stats(unit),
        "raw_sub_means": means.tolist(),
        "unit_po": unit.tolist(),
    }


def _anti_hard_filter(po: np.ndarray, y: np.ndarray, drop_frac: float = 0.25) -> dict:
    """Drop top-PO rows (as if 'filter OOD out of the monitor')."""
    po = np.asarray(po, dtype=float)
    y = np.asarray(y, dtype=int)
    n = len(po)
    k = max(1, int(drop_frac * n))
    keep = np.ones(n, dtype=bool)
    keep[np.argsort(po)[-k:]] = False
    y_k, po_k = y[keep], po[keep]
    return {
        "name": "hard_filter_drop_top_po",
        "drop_frac": float(drop_frac),
        "n_before": int(n),
        "n_after": int(keep.sum()),
        "mean_po_before": float(np.mean(po)),
        "mean_po_after": float(np.mean(po_k)) if keep.any() else float("nan"),
        "pos_rate_before": float(np.mean(y)),
        "pos_rate_after": float(np.mean(y_k)) if keep.any() else float("nan"),
        "selection_bias_d_pos_rate": float(np.mean(y_k) - np.mean(y)) if keep.any() else float("nan"),
        "why_bad": (
            "Removes the high-PO tail that *is* the malignant signal; "
            "also shifts class mix → selection bias for uplift / R_PO."
        ),
    }


def _anti_pseudo_label(det, texts: Sequence[str], y_true: np.ndarray) -> dict:
    """Pseudo-label: high PO ⇒ flip Y (treating uncertainty as 'wrong label')."""
    y_true = np.asarray(y_true, dtype=int)
    texts = list(texts)
    po = det.po_proxy(texts, y_true)
    thr = float(np.percentile(po, 75))
    y_pseudo = y_true.copy()
    flip = po >= thr
    y_pseudo[flip] = 1 - y_pseudo[flip]
    p = det.predict_proba(texts)
    thr_d = float(getattr(det, "threshold", 0.5))
    yhat = (p >= thr_d).astype(int)
    err_true = float(1.0 - (yhat == y_true).mean())
    err_pseudo = float(1.0 - (yhat == y_pseudo).mean())
    correct = yhat == y_true
    noise_on_correct = float(np.mean(flip & correct)) if flip.any() else 0.0
    return {
        "name": "pseudo_label_flip_high_po",
        "po_thr": thr,
        "flip_rate": float(np.mean(flip)),
        "err_vs_true_Y": err_true,
        "err_vs_pseudo_Y": err_pseudo,
        "noise_on_correct_rate": noise_on_correct,
        "label_agreement": float(np.mean(y_pseudo == y_true)),
        "why_bad": (
            "High PO on benign OOD is uniform calibration stress, not wrong Y. "
            "Flipping invents label noise; malignant FN rows get 'fixed' on paper "
            "while true Y still fails — fake uplift."
        ),
    }


def run_mvp(
    max_pairs: int = 900,
    live_n: int = 220,
    K: int = 36,
    seed: int = SEED,
    out_dir: Optional[Path] = None,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)

    data = load_hc3_binary(max_pairs=max_pairs, seed=seed)
    split = split_ref_live(data, live_frac=0.35, seed=seed)
    X_tr, y_tr = split["train"]
    X_live, y_live = split["live"]

    n_fit = min(len(X_tr), 1200)
    det = fit_sklearn_detector(X_tr[:n_fit], y_tr[:n_fit], seed=seed)

    texts_l = list(X_live[:live_n])
    y_l = np.asarray(y_live[:live_n], dtype=int)
    texts_r = list(X_tr[:live_n])
    y_r = np.asarray(y_tr[:live_n], dtype=int)

    # Benign = uniform mild synonym; malignant = mixture attack on a subset
    # (heteroscedastic peaks). This is the separable geometry.
    live = build_live_batches(
        texts_l,
        y_l,
        seed=seed,
        synonym_rate=0.12,
        attack_rate=0.40,
        attack_max=max(24, live_n // 3),
        detector=det,
    )
    batches: Dict[str, Tuple[list, np.ndarray]] = {
        "ref": (texts_r, y_r),
        "benign_ood": (
            list(live["benign_synonym"]["texts"]),
            np.asarray(live["benign_synonym"]["labels"], dtype=int),
        ),
        "malignant_ood": (
            list(live["malignant_attack"]["texts"]),
            np.asarray(live["malignant_attack"]["labels"], dtype=int),
        ),
    }

    panel: Dict[str, dict] = {}
    for name, (tx, yy) in batches.items():
        _sub_means, stats = _subbatch_risk_var(
            tx, yy, det, K=K, seed=seed + (abs(hash(name)) % 10000)
        )
        diag = diagnose_batch(tx, yy, det, K_subbatches=K, seed=seed + 7, kind=name)
        disc = discriminate_ood_vs_intentional(
            det,
            texts_r,
            y_r,
            tx,
            yy,
            K_subbatches=max(20, K // 2),
            seed=seed,
            label=name,
        )
        metrics = evaluate_detector(det, tx, yy)
        panel[name] = {
            "diag": diag,
            "stats": stats,
            "disc": disc,
            "route": "hold",
            "n": len(tx),
            "eval": asdict(metrics),
        }

    ref_diag = panel["ref"]["diag"]
    for name in ("benign_ood", "malignant_ood"):
        delta = {
            "d_mean": panel[name]["diag"].mean - ref_diag.mean,
            "d_peak": panel[name]["diag"].peak_gap - ref_diag.peak_gap,
            "d_cvar": panel[name]["diag"].cvar - ref_diag.cvar,
            "d_sub_var": panel[name]["stats"]["sub_var"] - panel["ref"]["stats"]["sub_var"],
            "d_sub_std": panel[name]["stats"]["sub_std"] - panel["ref"]["stats"]["sub_std"],
            "d_unit_var": panel[name]["stats"]["unit"]["var"] - panel["ref"]["stats"]["unit"]["var"],
            "d_top20": (
                panel[name]["stats"]["unit"]["top20_share"]
                - panel["ref"]["stats"]["unit"]["top20_share"]
            ),
        }
        panel[name]["delta"] = delta
        d_err = float(panel[name]["disc"]["delta"].get("d_error", 0.0))
        panel[name]["delta"]["d_error"] = d_err
        # Legacy mean-first router (can over-call malignant when ΔMean is large)
        panel[name]["route_legacy"] = route_action(panel[name]["diag"], ref=ref_diag)
        # Shape-first: Magnitude → Variance/error (subStd, unitVar, FN/error with Y)
        panel[name]["route"] = _shape_first_route(delta, d_error=d_err)
        panel[name]["verdict_shape"] = _shape_verdict(delta, d_error=d_err)

    anti = {}
    for name in ("benign_ood", "malignant_ood"):
        po = np.asarray(panel[name]["stats"]["unit_po"], dtype=float)
        tx, yy = batches[name]
        anti[name] = {
            "hard_filter": _anti_hard_filter(po, yy, drop_frac=0.25),
            "pseudo_label": _anti_pseudo_label(det, tx, yy),
        }

    det_c, _meta_c = conformal_calibrate(
        det, batches["benign_ood"][0], batches["benign_ood"][1], target_fpr=0.08
    )
    det_f = adapt_freeze_lexical(
        det, batches["malignant_ood"][0], batches["malignant_ood"][1], seed=seed
    )
    remedies = {}
    for tag, d_after, batch_key in (
        ("conformal_on_benign", det_c, "benign_ood"),
        ("freeze_on_malignant", det_f, "malignant_ood"),
    ):
        tx, yy = batches[batch_key]
        before = evaluate_detector(det, tx, yy)
        after = evaluate_detector(d_after, tx, yy)
        remedies[tag] = {
            "batch": batch_key,
            "acc_before": float(before.accuracy),
            "acc_after": float(after.accuracy),
            "auc_before": float(before.auc),
            "auc_after": float(after.auc),
        }

    def _pack(name: str) -> dict:
        p = panel[name]
        return {
            "route": p["route"],
            "route_legacy": p.get("route_legacy", p["route"]),
            "verdict": p.get("verdict_shape", p["disc"]["verdict"]),
            "verdict_disc": p["disc"]["verdict"],
            "actuator_suggested": p["disc"]["actuator"],
            "diag": {
                "mean": p["diag"].mean,
                "cvar": p["diag"].cvar,
                "peak_gap": p["diag"].peak_gap,
                "kurtosis": p["diag"].kurtosis,
            },
            "po_unit": p["stats"]["unit"],
            "po_subenv": {
                "sub_mean": p["stats"]["sub_mean"],
                "sub_std": p["stats"]["sub_std"],
                "sub_var": p["stats"]["sub_var"],
                "peak_gap": p["stats"]["peak_gap"],
                "cvar95": p["stats"]["cvar95"],
            },
            "delta_vs_ref": p.get("delta"),
            "disc_delta": p["disc"]["delta"],
            "y_live": p["disc"]["y_metrics_live"],
            "raw_sub_means": p["stats"]["raw_sub_means"],
            "unit_po": p["stats"]["unit_po"],
            "eval": p["eval"],
        }

    payload = {
        "mvp": "benign_vs_malignant_ood_po_variance",
        "justify": {
            "rule": (
                "Magnitude (ΔMean PO) first; if quiet → hold. "
                "If Mean↑, check Variance/Peak/CVaR/top-share: "
                "flat → benign OOD → conformal; spiked → malignant → freeze-lexical."
            ),
            "separable": True,
            "anti_patterns": [
                "hard_filter_drop_top_po — selection bias, kills malignant signal",
                "pseudo_label_flip_high_po — label noise; fake metric lifts",
            ],
            "same": "Both can raise mean PO / OOD score",
            "diff": "Benign = uniform; Malignant = heavy-tail / high subenv variance",
        },
        "batches": {k: _pack(k) for k in ("ref", "benign_ood", "malignant_ood")},
        "anti_patterns": anti,
        "remedies": remedies,
    }

    png = out_dir / "benign_vs_malignant_ood.png"
    plot_mvp(payload, png)
    (out_dir / "benign_vs_malignant_ood.json").write_text(json.dumps(payload, indent=2))
    _write_readme(out_dir / "README.md", payload, png.name)
    return payload


def plot_mvp(payload: dict, path: Path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    colors = {"ref": "#5B6B7A", "benign_ood": "#2F6B4F", "malignant_ood": "#A33B24"}
    labels = {"ref": "REF", "benign_ood": "Benign OOD", "malignant_ood": "Malignant OOD"}

    fig = plt.figure(figsize=(13.2, 9.6))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.32)
    keys = ["ref", "benign_ood", "malignant_ood"]

    # 1) Mean × subenv-std
    ax = fig.add_subplot(gs[0, 0])
    for k in keys:
        b = payload["batches"][k]
        ax.scatter(
            b["diag"]["mean"],
            b["po_subenv"]["sub_std"],
            s=160,
            c=colors[k],
            label=labels[k],
            edgecolors="white",
            linewidths=0.8,
            zorder=3,
        )
        ax.annotate(
            labels[k],
            (b["diag"]["mean"], b["po_subenv"]["sub_std"]),
            fontsize=8,
            xytext=(7, 5),
            textcoords="offset points",
        )
    ax.set_xlabel("Mean PO-risk (magnitude)")
    ax.set_ylabel("Sub-env Std of PO (variance)")
    ax.set_title("Magnitude × Variance — the separator", loc="left", fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7.5)

    # 2) Peak + CVaR
    ax = fig.add_subplot(gs[0, 1])
    x = np.arange(len(keys))
    peaks = [payload["batches"][k]["diag"]["peak_gap"] for k in keys]
    cvars = [payload["batches"][k]["diag"]["cvar"] for k in keys]
    ax.bar(x - 0.18, peaks, width=0.36, color=[colors[k] for k in keys], label="Peak Gap")
    ax.bar(x + 0.18, cvars, width=0.36, color=[colors[k] for k in keys], alpha=0.45, label="CVaR95")
    ax.set_xticks(x)
    ax.set_xticklabels([labels[k] for k in keys], fontsize=8)
    ax.set_ylabel("PO shape")
    ax.set_title("Peaks: Peak Gap & CVaR", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    # 3) Unit PO histogram
    ax = fig.add_subplot(gs[0, 2])
    for k in keys:
        po = np.asarray(payload["batches"][k]["unit_po"], dtype=float)
        ax.hist(po, bins=24, density=True, alpha=0.45, color=colors[k], label=labels[k])
    ax.set_xlabel("Unit PO = 1 − p̂(Y|x)")
    ax.set_ylabel("density")
    ax.set_title("Unit PO mass — uniform vs heavy tail", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # 4) Sub-env strip
    ax = fig.add_subplot(gs[1, 0])
    for i, k in enumerate(keys):
        vals = np.asarray(payload["batches"][k]["raw_sub_means"], dtype=float)
        jitter = (np.random.default_rng(i + 3).random(len(vals)) - 0.5) * 0.18
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=14, alpha=0.55, c=colors[k])
        ax.hlines(np.median(vals), i - 0.25, i + 0.25, colors="black", lw=1.4)
    ax.set_xticks(range(3))
    ax.set_xticklabels([labels[k] for k in keys], fontsize=8)
    ax.set_ylabel("Bootstrap sub-env mean PO")
    ax.set_title("Sub-environment spread (OOD variance)", loc="left", fontweight="bold")
    ax.grid(True, axis="y", alpha=0.3)

    # 5) Grouped Δ metrics — the justify panel
    ax = fig.add_subplot(gs[1, 1])
    metrics = ["d_mean", "d_sub_std", "d_unit_var", "d_error"]
    metric_labs = ["ΔMean", "ΔSubStd", "ΔUnitVar", "ΔError"]
    xb = np.arange(len(metrics))
    w = 0.36
    ben_d = payload["batches"]["benign_ood"]["delta_vs_ref"]
    mal_d = payload["batches"]["malignant_ood"]["delta_vs_ref"]
    # scale sub_std / unit_var for visibility beside mean/error
    scale = {"d_mean": 1.0, "d_sub_std": 8.0, "d_unit_var": 4.0, "d_error": 1.0}
    ax.bar(
        xb - w / 2,
        [ben_d.get(m, 0.0) * scale[m] for m in metrics],
        width=w,
        color=colors["benign_ood"],
        label="Benign",
    )
    ax.bar(
        xb + w / 2,
        [mal_d.get(m, 0.0) * scale[m] for m in metrics],
        width=w,
        color=colors["malignant_ood"],
        label="Malignant",
    )
    ax.set_xticks(xb)
    ax.set_xticklabels(metric_labs, fontsize=8)
    ax.set_ylabel("Δ vs REF (SubStd×8, UnitVar×4)")
    ax.set_title("Justify separators (Mean alone is not enough)", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    # 6) Anti-pattern panel
    ax = fig.add_subplot(gs[1, 2])
    ax.axis("off")
    ax.set_title("Anti-patterns vs correct route", loc="left", fontweight="bold")
    mal = payload["anti_patterns"]["malignant_ood"]
    ben = payload["anti_patterns"]["benign_ood"]
    lines = [
        "HARD FILTER (drop top-25% PO)",
        "  malignant Δpos_rate: %+.3f" % mal["hard_filter"]["selection_bias_d_pos_rate"],
        "  benign    Δpos_rate: %+.3f" % ben["hard_filter"]["selection_bias_d_pos_rate"],
        "  → selection bias; kills attack tail",
        "",
        "PSEUDO-LABEL (flip high-PO Y)",
        "  malignant flip_rate: %.2f" % mal["pseudo_label"]["flip_rate"],
        "  noise_on_correct:   %.2f" % mal["pseudo_label"]["noise_on_correct_rate"],
        "  label_agreement:    %.2f" % mal["pseudo_label"]["label_agreement"],
        "  → fake metric lift, true Y still wrong",
        "",
        "CORRECT (this MVP)",
        "  benign  → %s / %s"
        % (
            payload["batches"]["benign_ood"]["route"],
            payload["batches"]["benign_ood"]["verdict"],
        ),
        "  malign  → %s / %s"
        % (
            payload["batches"]["malignant_ood"]["route"],
            payload["batches"]["malignant_ood"]["verdict"],
        ),
    ]
    ax.text(
        0.02,
        0.98,
        "\n".join(lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7.8,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F7F4EF", edgecolor="#CCC"),
    )

    fig.suptitle(
        "Benign vs Malignant OOD MVP — PO magnitude + variance (not hard-filter / pseudo-Y)",
        fontsize=12.5,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.004,
        "Same: Mean PO can rise. Diff: sub-env variance / Peak / CVaR / top-share spike only under malignant injection.",
        ha="center",
        fontsize=8,
        style="italic",
        color="#444",
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _write_readme(path: Path, payload: dict, png_name: str):
    b = payload["batches"]
    lines = [
        "# Benign vs Malignant OOD MVP (PO mean + variance)",
        "",
        "## Justify",
        "",
        payload["justify"]["rule"],
        "",
        "- **Same:** " + payload["justify"]["same"],
        "- **Diff:** " + payload["justify"]["diff"],
        "",
        "## Separation panel",
        "",
        "| batch | mean PO | subStd | peak | CVaR | verdict | route |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for k in ("ref", "benign_ood", "malignant_ood"):
        x = b[k]
        lines.append(
            "| %s | %.4f | %.4f | %.4f | %.4f | %s | %s |"
            % (
                k,
                x["diag"]["mean"],
                x["po_subenv"]["sub_std"],
                x["diag"]["peak_gap"],
                x["diag"]["cvar"],
                x["verdict"],
                x["route"],
            )
        )
    lines.extend(
        [
            "",
            "## Anti-patterns",
            "",
            "Do **not** hard-filter high-PO rows or pseudo-label flips from PO.",
            "See dashboard right panel for selection-bias / label-noise numbers.",
            "",
            f"![]({png_name})",
            "",
            "```bash",
            "python3 scripts/run_as_dta_benign_malignant_mvp.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
