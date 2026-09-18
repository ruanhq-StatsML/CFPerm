#!/usr/bin/env python3
"""Obs-PO hard-reweight v7: hard-support + beijing-class packMSE gate.

Logic we buy
------------
1. **Hard-rank always** — obs PO ranks hard rows (Spearman / P@20%).
   Mechanism that matches: ``hard_support`` (top-k boost scaled by λ).
2. **Pack MSE only under beijing-class drift** — chase all-row next-MSE
   only when drift ≫ mild (gate≈0.45). Calm rejects stay uniform for packMSE.
3. Primary claim metric = **hard-subset next-MSE**; packMSE is conditional.
4. **v5**: CV-MSE picks PO^power × temper (or family) per reject; add log1p map.
5. **v6**: hard-objective CV; soft power grid; dual/blend hard+qrt under beijing.
6. **v7**: blend mix scan {0.25,0.5,0.75}; dual with n_recent=2.

  PYTHONPATH=. python3 scripts/run_agod_obs_po_weights.py \\
    --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT stocks_IWM waymo_proxy
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor

from agod.hard_rank_metrics import hard_rank_metrics
from agod.obs_po_cv import SOFT_POWERS, cv_select_power
from agod.obs_po_weights import (
    BEIJING_DRIFT_GATE,
    adaptive_temper,
    blend_hard_qrt_weights,
    drift_intensity,
    gated_obs_po_weights,
    hard_subset_mask,
    is_beijing_class_drift,
)
from agod.online_rfperm import fit_online_rfperm, update_online_rfperm
from agod.po_refit import build_recent_ood_windows, refit_po_on_windows
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

# v4 modes: hard-claim path + beijing packMSE path + v3 baselines
MODES: Tuple[str, ...] = (
    "uniform",
    "gated_hard_adapt",  # hard claim default
    "dual",              # v6 winner: mild→hard, BJ→soft CV
    "dual_r2",           # NEW: dual + n_recent=2
    "blend_25",          # NEW: mix scan (more qrt)
    "blend_50",          # mix=0.5 (v6 blend_bj)
    "blend_75",          # NEW: mix scan (more hard_support)
)
COLORS = {
    "uniform": "#4C566A",
    "gated_hard_adapt": "#88C0D0",
    "dual": "#5E81AC",
    "dual_r2": "#81A1C1",
    "blend_25": "#EBCB8B",
    "blend_50": "#D08770",
    "blend_75": "#BF616A",
}


def ensure_loaders() -> None:
    for t in ("MSFT", "IWM", "AAPL", "SPY", "QQQ"):
        key = f"stocks_{t}"
        if key not in LOADERS:
            LOADERS[key] = (
                lambda tk: (lambda root, max_n=20000: load_stocks(root, tk, max_n))
            )(t)


def load_xy(name: str, root: Path, max_n: int = 20000):
    ensure_loaders()
    out = LOADERS[name](root, max_n=max_n)
    return np.asarray(out[0], float), np.asarray(out[1], float).ravel()


def pca_fit(X, d, seed):
    X = np.asarray(X, float)
    if d and X.shape[1] > d:
        return PCA(n_components=d, random_state=seed).fit_transform(X).astype(np.float32)
    return X.astype(np.float32)


def make_stream(X, y, bs: int, n_batches: int):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, len(X), bs)]


def fit_rf(X, y, w, seed):
    m = RandomForestRegressor(
        n_estimators=40, max_depth=8, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    m.fit(X, y, sample_weight=w)
    return m


def mode_spec(mode: str) -> dict:
    """kind / obs map / drift gate / n_recent."""
    if mode == "uniform":
        return {"kind": "fixed", "obs": "uniform", "adapt": False, "drift_gate": 0.0, "n_recent": None}
    if mode == "gated_hard_adapt":
        return {
            "kind": "fixed",
            "obs": "hard_support",
            "adapt": True,
            "drift_gate": 0.20,
            "n_recent": None,
        }
    if mode == "dual":
        return {"kind": "dual", "n_recent": None, "bj_objective": "all"}
    if mode == "dual_r2":
        return {"kind": "dual", "n_recent": 2, "bj_objective": "all"}
    if mode == "blend_25":
        return {"kind": "blend", "drift_gate": BEIJING_DRIFT_GATE, "mix": 0.25, "n_recent": None}
    if mode == "blend_50":
        return {"kind": "blend", "drift_gate": BEIJING_DRIFT_GATE, "mix": 0.50, "n_recent": None}
    if mode == "blend_75":
        return {"kind": "blend", "drift_gate": BEIJING_DRIFT_GATE, "mix": 0.75, "n_recent": None}
    raise KeyError(mode)


def run_mode(
    stream,
    mode: str,
    seed: int,
    *,
    n_burn: int,
    alpha: float,
    n_recent_default: int,
    hard_frac: float = 0.2,
    cv_folds: int = 3,
) -> dict:
    spec = mode_spec(mode)
    mse_next: List[float] = []
    mse_next_hard: List[float] = []
    mse_next_easy: List[float] = []
    gate_on: List[int] = []
    drifts: List[float] = []
    lams: List[float] = []
    powers: List[float] = []
    families: List[str] = []
    beijing_on: List[int] = []
    hard_rows: List[dict] = []

    X0, y0 = stream[0]
    rfperm = fit_online_rfperm(X0, y0, seed=seed)

    for t in range(1, len(stream)):
        Xc, yc = stream[t]
        burn = t <= n_burn
        step = update_online_rfperm(
            rfperm, Xc, yc, burn_in=burn, alpha=alpha, ewma=True, fdr="alpha_investing"
        )
        rejected = bool(step["reject"]) and not burn
        gate_on.append(int(rejected))

        n_recent = spec["n_recent"] or n_recent_default
        windows = None
        po = None
        lam_used = 0.0
        drift = 0.0
        power_used = 0.0
        fam_used = "uniform"

        if mode == "uniform" or not rejected:
            w = np.ones(len(yc), float)
        else:
            windows = build_recent_ood_windows(
                stream, t, n_recent=n_recent, window_mode="recent_ood"
            )
            po = refit_po_on_windows(windows, seed=seed + t, blend_mu_gap=0.25)
            # Drift vs fixed f_ref (RFPerm yardstick) — not in-sample μ0.
            f_ref = rfperm.f_ref
            control_resid = np.abs(
                windows.y_recent - f_ref.predict(windows.X_recent)
            )
            ood_resid_ref = np.abs(yc - f_ref.predict(Xc))
            drift = drift_intensity(
                ood_resid_ref,
                control_resid=control_resid,
                p=float(step["p"]),
                T=float(step.get("T", 0.0)),
                alpha=alpha,
            )
            drifts.append(drift)
            beijing_on.append(int(is_beijing_class_drift(drift)))
            power_used = 0.0
            fam_used = "uniform"
            kind = spec.get("kind", "fixed")

            if kind == "fixed":
                if spec["adapt"]:
                    lam_used = adaptive_temper(
                        drift, lam_max=0.75, drift_gate=float(spec["drift_gate"])
                    )
                else:
                    lam_used = 0.0
                w = gated_obs_po_weights(
                    po,
                    reject=True,
                    mode=spec["obs"],  # type: ignore[arg-type]
                    soft=False,
                    p=float(step["p"]),
                    alpha=alpha,
                    temper=lam_used,
                    topk_frac=hard_frac,
                    boost_max=3.0,
                )
                fam_used = str(spec["obs"])
                power_used = 0.25 if spec["obs"] == "qrt" else 0.0

            elif kind == "cv_power":
                cap = None
                if spec.get("temper_cap"):
                    cap = adaptive_temper(
                        drift,
                        lam_max=0.75,
                        drift_gate=float(spec.get("drift_gate", BEIJING_DRIFT_GATE)),
                    )
                    if cap <= 0.0:
                        w = np.ones(len(yc), float)
                        lam_used = 0.0
                        power_used = 0.0
                        fam_used = "uniform"
                        lams.append(lam_used)
                        powers.append(power_used)
                        families.append(fam_used)
                        model = fit_rf(Xc, yc, w, seed + t)
                        if t + 1 < len(stream):
                            Xn, yn = stream[t + 1]
                            pred = model.predict(Xn)
                            err2 = (yn - pred) ** 2
                            mse_next.append(float(np.mean(err2)))
                            Xp, yp = stream[t - 1]
                            mu0 = fit_rf(Xp, yp, np.ones(len(yp)), seed + 99 + t)
                            po_next = np.abs(yn - mu0.predict(Xn))
                            hard_m = hard_subset_mask(po_next, frac=hard_frac)
                            mse_next_hard.append(float(np.mean(err2[hard_m])))
                            mse_next_easy.append(float(np.mean(err2[~hard_m])))
                        continue
                powers_grid = SOFT_POWERS if spec.get("powers") == "soft" else None
                sel_kw = dict(
                    n_folds=cv_folds,
                    seed=seed + t,
                    temper_cap=cap,
                    objective=spec.get("objective", "all"),
                    hard_frac=hard_frac,
                )
                if powers_grid is not None:
                    sel_kw["powers"] = powers_grid
                sel = cv_select_power(Xc, yc, po, **sel_kw)
                w = np.asarray(sel["weights"], float)
                lam_used = float(sel["temper"])
                power_used = float(sel["power"])
                fam_used = f"PO^{power_used:g}/{sel.get('objective','all')}"

            elif kind == "blend":
                lam_used = adaptive_temper(
                    drift, lam_max=0.75, drift_gate=float(spec["drift_gate"])
                )
                w = blend_hard_qrt_weights(
                    po,
                    lam=lam_used,
                    mix=float(spec.get("mix", 0.5)),
                    topk_frac=hard_frac,
                    boost_max=3.0,
                )
                fam_used = "blend_hard_qrt"
                power_used = 0.25

            elif kind == "dual":
                # mild reject → hard_support; beijing → soft CV power cap
                if is_beijing_class_drift(drift):
                    cap = adaptive_temper(
                        drift, lam_max=0.75, drift_gate=BEIJING_DRIFT_GATE
                    )
                    sel = cv_select_power(
                        Xc,
                        yc,
                        po,
                        powers=SOFT_POWERS,
                        n_folds=cv_folds,
                        seed=seed + t,
                        temper_cap=cap,
                        objective=spec.get("bj_objective", "all"),
                        hard_frac=hard_frac,
                    )
                    w = np.asarray(sel["weights"], float)
                    lam_used = float(sel["temper"])
                    power_used = float(sel["power"])
                    fam_used = f"dual_cv_PO^{power_used:g}"
                else:
                    lam_used = adaptive_temper(drift, lam_max=0.75, drift_gate=0.20)
                    w = gated_obs_po_weights(
                        po,
                        reject=True,
                        mode="hard_support",
                        soft=False,
                        p=float(step["p"]),
                        alpha=alpha,
                        temper=lam_used,
                        topk_frac=hard_frac,
                        boost_max=3.0,
                    )
                    fam_used = "dual_hard"
                    power_used = 0.0
            else:
                raise RuntimeError(kind)

            lams.append(lam_used)
            powers.append(power_used)
            families.append(fam_used)

        model = fit_rf(Xc, yc, w, seed + t)

        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            pred = model.predict(Xn)
            err2 = (yn - pred) ** 2
            mse_next.append(float(np.mean(err2)))

            if windows is not None:
                mu0 = fit_rf(
                    windows.X_recent,
                    windows.y_recent,
                    np.ones(len(windows.y_recent)),
                    seed + 99 + t,
                )
            else:
                Xp, yp = stream[t - 1]
                mu0 = fit_rf(Xp, yp, np.ones(len(yp)), seed + 99 + t)
            po_next = np.abs(yn - mu0.predict(Xn))
            hard_m = hard_subset_mask(po_next, frac=hard_frac)
            mse_next_hard.append(float(np.mean(err2[hard_m])))
            mse_next_easy.append(float(np.mean(err2[~hard_m])))

        if rejected and po is not None and windows is not None:
            Xr = np.vstack([windows.X_recent, windows.X_ood])
            yr = np.concatenate([windows.y_recent, windows.y_ood])
            mu_o = fit_rf(Xr, yr, np.ones(len(yr)), seed + 11 + t)
            truth = np.abs(yc - mu_o.predict(Xc))
            hard_rows.append(
                {
                    "po": hard_rank_metrics(po, truth),
                    "drift": drift,
                    "lam": lam_used,
                    "power": power_used if rejected else 0.0,
                    "family": fam_used if rejected else "uniform",
                }
            )

    def _avg(xs: List[float]) -> float:
        return float(np.mean(xs)) if xs else float("nan")

    def _mean_hard() -> dict:
        if not hard_rows:
            return {"spearman": float("nan"), "precision_at_k": float("nan"), "n": 0}
        sp = [h["po"]["spearman"] for h in hard_rows if h["po"]["spearman"] == h["po"]["spearman"]]
        pk = [
            h["po"]["precision_at_k"]
            for h in hard_rows
            if h["po"]["precision_at_k"] == h["po"]["precision_at_k"]
        ]
        return {
            "spearman": float(np.mean(sp)) if sp else float("nan"),
            "precision_at_k": float(np.mean(pk)) if pk else float("nan"),
            "n": len(hard_rows),
        }

    return {
        "mse_next": mse_next,
        "mse_next_hard": mse_next_hard,
        "mse_next_easy": mse_next_easy,
        "gate_on": gate_on,
        "mse_mean": _avg(mse_next),
        "mse_hard_mean": _avg(mse_next_hard),
        "mse_easy_mean": _avg(mse_next_easy),
        "hard_po": _mean_hard(),
        "n_reject": int(sum(gate_on)),
        "duty": float(np.mean(gate_on)) if gate_on else 0.0,
        "drift_mean": _avg(drifts),
        "lam_mean": _avg(lams),
        "power_mean": _avg(powers),
        "family_mode": (
            (lambda vals, counts: str(vals[int(np.argmax(counts))]))(
                *np.unique(np.asarray(families), return_counts=True)
            )
            if families
            else "—"
        ),
        "beijing_frac": _avg([float(x) for x in beijing_on]),
    }


def run_dataset(
    name: str,
    root: Path,
    *,
    seed: int,
    batch_size: int,
    n_batches: int,
    pca_d: int,
    n_burn: int,
    alpha: float,
    n_recent: int,
    cv_folds: int = 3,
) -> dict:
    X, y = load_xy(name, root)
    X = pca_fit(X, pca_d, seed)
    stream = make_stream(X, y, batch_size, n_batches)
    results = {}
    for mode in MODES:
        results[mode] = run_mode(
            stream,
            mode,
            seed,
            n_burn=n_burn,
            alpha=alpha,
            n_recent_default=n_recent,
            cv_folds=cv_folds,
        )
    annotate_results_with_sig(
        results,
        preferred_gate_modes=(
            "gated_hard_adapt",
            "dual",
            "dual_r2",
            "blend_25",
            "blend_50",
            "blend_75",
        ),
    )
    gate = None
    for m in (
        "gated_hard_adapt",
        "dual",
        "dual_r2",
        "uniform",
    ):
        if m in results and results[m].get("gate_on"):
            g = np.asarray(results[m]["gate_on"], bool)
            n = len(results[m]["mse_next"])
            gate = g[:n]
            break
    if gate is not None and gate.any():
        for m, r in results.items():
            mh = np.asarray(r["mse_next_hard"], float)[: len(gate)]
            me = np.asarray(r["mse_next_easy"], float)[: len(gate)]
            r["mse_hard_mean_sig"] = float(mh[gate].mean()) if len(mh) == len(gate) else float("nan")
            r["mse_easy_mean_sig"] = float(me[gate].mean()) if len(me) == len(gate) else float("nan")
    return {
        "dataset": name,
        "n_batches": len(stream),
        "results": results,
        "drift_mean": results["gated_hard_adapt"].get("drift_mean"),
        "beijing_class": is_beijing_class_drift(
            float(results["gated_hard_adapt"].get("drift_mean") or 0.0)
        ),
    }


def _bars(all_ds: dict, out: Path):
    packs = list(all_ds.keys())
    fig, axes = plt.subplots(1, 2, figsize=(max(11, 1.8 * len(packs)), 4.6))
    x = np.arange(len(packs))
    w = 0.11
    mid = (len(MODES) - 1) / 2.0

    for ax, key, title in [
        (axes[0], "mse_mean_sig", "sig-only next MSE (all rows / pack)"),
        (axes[1], "mse_hard_mean_sig", "sig-only next MSE (hard top-20%)"),
    ]:
        for i, m in enumerate(MODES):
            vals = []
            for p in packs:
                r = all_ds[p]["results"][m]
                vals.append(r.get(key, r.get("mse_mean")))
            ax.bar(x + (i - mid) * w, vals, w, label=m, color=COLORS[m])
        ax.set_xticks(x)
        ax.set_xticklabels(packs, rotation=18, ha="right")
        ax.set_ylabel("MSE (↓)")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].legend(ncol=2, fontsize=6)
    fig.tight_layout()
    path = out / "obs_po_v7_mse.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def report(all_ds: dict) -> str:
    def f(v, pct=False):
        if v is None or (isinstance(v, float) and v != v):
            return "—"
        return f"{100*v:+.1f}%" if pct else f"{v:.4g}"

    short = {m: short_name(m) for m in MODES}
    lines = [
        "# Observation-level PO hard-reweight v7 (blend mix + dual_r2)",
        "",
        "## What changed vs v6",
        "",
        "- Keep **`dual`** (v6 beijing packMSE winner) and **`hard_m`**.",
        "- **`dual_r2`**: same dual policy with `n_recent=2` for PO/μ0 windows.",
        "- **Blend mix scan**: `blend_25` / `blend_50` / `blend_75` (hard_support vs qrt).",
        "- Drop cv_hard_cap from primary table (did not transfer in v6).",
        "",
        "## Drift / dual diagnostics",
        "",
        "| dataset | drift | beijing? | dual mode | dual_r2 mode |",
        "|---|---:|:---:|---|---|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]
        bj = "yes" if blob.get("beijing_class") else "no"
        lines.append(
            f"| `{ds}` | {f(blob.get('drift_mean'))} | {bj} | "
            f"`{r['dual'].get('family_mode')}` | `{r['dual_r2'].get('family_mode')}` |"
        )

    hdr = " | ".join(short[m] for m in MODES)
    lines += [
        "",
        "## Sig-only hard top-20% next-MSE (↓)  ← primary",
        "",
        f"| dataset | {hdr} | best |",
        "|" + "---|---:" * len(MODES) + "|---|",
    ]
    wins_h = {m: 0 for m in MODES}
    for ds, blob in all_ds.items():
        r = blob["results"]

        def score_h(m, r=r):
            v = r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))
            return v if v == v else 1e99

        best = min(MODES, key=score_h)
        wins_h[best] += 1
        cells = [f(r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))) for m in MODES]
        lines.append(f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` |")
    lines += [
        "",
        "**Wins (hard):** " + ", ".join(f"`{short[k]}`={v}" for k, v in wins_h.items()),
        "",
        "## Sig-only pack next-MSE (↓)  ← beijing-conditional",
        "",
        f"| dataset | {hdr} | best | beijing? |",
        "|" + "---|---:" * len(MODES) + "|---|:---:|",
    ]
    wins = {m: 0 for m in MODES}
    wins_bj = {m: 0 for m in MODES}
    n_bj = 0
    for ds, blob in all_ds.items():
        r = blob["results"]
        bj = bool(blob.get("beijing_class"))

        def score(m, r=r):
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return v if v == v else 1e99

        best = min(MODES, key=score)
        wins[best] += 1
        if bj:
            n_bj += 1
            wins_bj[best] += 1
        cells = [f(r[m].get("mse_mean_sig", r[m]["mse_mean"])) for m in MODES]
        lines.append(
            f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` | {'yes' if bj else 'no'} |"
        )
    lines += [
        "",
        "**Wins (all):** " + ", ".join(f"`{short[k]}`={v}" for k, v in wins.items()),
        f"**Wins (beijing only, n={n_bj}):** "
        + ", ".join(f"`{short[k]}`={v}" for k, v in wins_bj.items()),
        "",
        "## Rel. pack MSE vs uniform",
        "",
        "| dataset | hard_m | dual | dual_r2 | b25 | b50 | b75 | drift |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]
        u = r["uniform"].get("mse_mean_sig", r["uniform"]["mse_mean"])

        def rel(m, r=r, u=u):
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return (v / u - 1.0) if u == u and u > 0 else float("nan")

        lines.append(
            f"| `{ds}` | {f(rel('gated_hard_adapt'), pct=True)} | "
            f"{f(rel('dual'), pct=True)} | {f(rel('dual_r2'), pct=True)} | "
            f"{f(rel('blend_25'), pct=True)} | {f(rel('blend_50'), pct=True)} | "
            f"{f(rel('blend_75'), pct=True)} | {f(blob.get('drift_mean'))} |"
        )

    lines += [
        "",
        "## Hard-rank",
        "",
        "| dataset | spearman | P@20% | n_reject |",
        "|---|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        h = blob["results"]["gated_hard_adapt"]["hard_po"]
        nrej = blob["results"]["gated_hard_adapt"]["n_reject"]
        lines.append(f"| `{ds}` | {f(h['spearman'])} | {f(h['precision_at_k'])} | {nrej} |")

    lines += [
        "",
        "### Takeaway",
        "",
        "- If `dual_r2` ≥ `dual` on beijing packMSE, prefer n_recent=2 in the recipe.",
        "- Pick blend mix by beijing packMSE / hard-subset tradeoff; lock best mix.",
        "- Default unified policy remains dual-family unless mix clearly dominates.",
        "",
        "See `docs/agod/AGOD_obs_po_weights.md`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=[
            "metro_interstate",
            "beijing_pm25",
            "stocks_AAPL",
            "stocks_MSFT",
            "stocks_IWM",
            "waymo_proxy",
        ],
    )
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--pca-d", type=int, default=16)
    ap.add_argument("--n-burn", type=int, default=5)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--n-recent", type=int, default=1)
    ap.add_argument("--cv-folds", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_obs_po_weights"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    all_ds: Dict[str, dict] = {}
    for name in args.datasets:
        print(f"=== {name} ===", flush=True)
        try:
            blob = run_dataset(
                name,
                args.root,
                seed=args.seed,
                batch_size=args.batch_size,
                n_batches=args.n_batches,
                pca_d=args.pca_d,
                n_burn=args.n_burn,
                alpha=args.alpha,
                n_recent=args.n_recent,
                cv_folds=args.cv_folds,
            )
        except Exception as e:
            print(f"[skip] {name}: {e}", flush=True)
            continue
        all_ds[name] = blob
        r = blob["results"]
        bj = "BJ" if blob.get("beijing_class") else "calm"
        print(
            "  drift={:.3f} ({}) dual={} dual_r2={} | sigPack ".format(
                float(blob.get("drift_mean") or float("nan")),
                bj,
                r["dual"].get("family_mode"),
                r["dual_r2"].get("family_mode"),
            )
            + " ".join(
                f"{short_name(m)}={r[m].get('mse_mean_sig', r[m]['mse_mean']):.4g}" for m in MODES
            )
            + " | hard "
            + " ".join(
                f"{short_name(m)}={r[m].get('mse_hard_mean_sig', r[m].get('mse_hard_mean', float('nan'))):.4g}"
                for m in MODES
            ),
            flush=True,
        )

    payload = {
        "version": 7,
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "cv_folds": args.cv_folds,
        "modes": list(MODES),
        "beijing_drift_gate": BEIJING_DRIFT_GATE,
        "datasets": all_ds,
        "note": "v7 blend mix scan + dual n_recent=2",
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_ds)
    (args.out / "OBS_PO_WEIGHTS_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_obs_po_weights.md").write_text(md, encoding="utf-8")
    plots = [_bars(all_ds, args.out)] if all_ds else []
    print(md)
    print("plots:", [str(p) for p in plots])


def short_name(m: str) -> str:
    return {
        "uniform": "uni",
        "gated_hard_adapt": "hard_m",
        "dual": "dual",
        "dual_r2": "dual_r2",
        "blend_25": "b25",
        "blend_50": "b50",
        "blend_75": "b75",
    }.get(m, m)


if __name__ == "__main__":
    main()
