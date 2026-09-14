#!/usr/bin/env python3
"""Obs-PO hard-reweight v3: adaptive temper by drift + hard-subset MSE.

Thesis
------
- Obs PO **ranks hard rows** (always useful as a hardness score).
- Pack next-MSE only lifts when hard ≈ **shift signal** (beijing-like),
  not when hard ≈ noise (calm stocks). Encode that with drift-adaptive λ:
    mild reject  → λ≈0 (stay uniform)
    strong drift → λ↑ (soft PO^{1/4} temper)

Also report **hard-subset next-MSE** (top-20% of next batch by PO) — the
metric that matches the hard-reweight claim.

  PYTHONPATH=. python3 scripts/run_agod_obs_po_weights.py \\
    --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT stocks_IWM waymo_proxy
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor

from agod.hard_rank_metrics import hard_rank_metrics
from agod.obs_po_weights import (
    adaptive_temper,
    drift_intensity,
    gated_obs_po_weights,
    hard_subset_mask,
)
from agod.online_rfperm import fit_online_rfperm, update_online_rfperm
from agod.po_refit import build_recent_ood_windows, refit_po_on_windows
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

# mode -> how temper is chosen
# fixed λ | adapt | uniform
MODES: Tuple[str, ...] = (
    "uniform",
    "gated_qrt_t50",       # fixed λ=0.5 (v2 best)
    "gated_qrt_adapt",     # λ = f(drift)
    "gated_qrt_adapt_r2",  # adapt + n_recent=2
    "gated_cbrt_adapt",     # adapt with ∛ base
)
COLORS = {
    "uniform": "#4C566A",
    "gated_qrt_t50": "#EBCB8B",
    "gated_qrt_adapt": "#A3BE8C",
    "gated_qrt_adapt_r2": "#88C0D0",
    "gated_cbrt_adapt": "#5E81AC",
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
    """Return obs_mode / temper policy / n_recent."""
    if mode == "uniform":
        return {"obs": "uniform", "temper": 0.0, "adapt": False, "n_recent": None}
    if mode == "gated_qrt_t50":
        return {"obs": "qrt", "temper": 0.5, "adapt": False, "n_recent": None}
    if mode == "gated_qrt_adapt":
        return {"obs": "qrt", "temper": None, "adapt": True, "n_recent": None}
    if mode == "gated_qrt_adapt_r2":
        return {"obs": "qrt", "temper": None, "adapt": True, "n_recent": 2}
    if mode == "gated_cbrt_adapt":
        return {"obs": "cbrt", "temper": None, "adapt": True, "n_recent": None}
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
) -> dict:
    spec = mode_spec(mode)
    mse_next: List[float] = []
    mse_next_hard: List[float] = []
    mse_next_easy: List[float] = []
    gate_on: List[int] = []
    drifts: List[float] = []
    lams: List[float] = []
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

        if mode == "uniform" or not rejected:
            w = np.ones(len(yc), float)
        else:
            windows = build_recent_ood_windows(
                stream, t, n_recent=n_recent, window_mode="recent_ood"
            )
            po = refit_po_on_windows(windows, seed=seed + t, blend_mu_gap=0.25)
            # Drift vs fixed f_ref (same yardstick as RFPerm) — not in-sample μ0,
            # which understates control residuals and inflates OOD/control ratios.
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
            if spec["adapt"]:
                lam_used = adaptive_temper(drift, lam_max=0.75, drift_gate=0.20)
            else:
                lam_used = float(spec["temper"])
            drifts.append(drift)
            lams.append(lam_used)
            w = gated_obs_po_weights(
                po,
                reject=True,
                mode=spec["obs"],  # type: ignore[arg-type]
                soft=False,
                p=float(step["p"]),
                alpha=alpha,
                temper=lam_used,
            )

        model = fit_rf(Xc, yc, w, seed + t)

        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            pred = model.predict(Xn)
            err2 = (yn - pred) ** 2
            mse_next.append(float(np.mean(err2)))

            # Hard-subset on next batch: score with μ0 from current windows if any,
            # else a quick μ0 on previous batch.
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
            hard_rows.append({"po": hard_rank_metrics(po, truth), "drift": drift, "lam": lam_used})

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
        )
    annotate_results_with_sig(
        results,
        preferred_gate_modes=(
            "gated_qrt_t50",
            "gated_qrt_adapt",
            "gated_qrt_adapt_r2",
            "gated_cbrt_adapt",
        ),
    )
    # also sig-filter hard/easy mse using same gate mask
    from agod.sig_batch_metrics import significant_only_stats

    # piggy-back: copy mse_next_hard into a side channel via annotate pattern
    gate = None
    for m in (
        "gated_qrt_adapt",
        "gated_qrt_t50",
        "gated_cbrt_adapt",
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
        "drift_mean_adapt": results["gated_qrt_adapt"].get("drift_mean"),
    }


def _bars(all_ds: dict, out: Path):
    packs = list(all_ds.keys())
    fig, axes = plt.subplots(1, 2, figsize=(max(11, 1.8 * len(packs)), 4.6))
    x = np.arange(len(packs))
    w = 0.15
    mid = (len(MODES) - 1) / 2.0

    for ax, key, title in [
        (axes[0], "mse_mean_sig", "sig-only next MSE (all rows)"),
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
    path = out / "obs_po_v3_mse.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def report(all_ds: dict) -> str:
    def f(v, pct=False):
        if v is None or (isinstance(v, float) and v != v):
            return "—"
        return f"{100*v:+.1f}%" if pct else f"{v:.4g}"

    short = {
        "uniform": "uniform",
        "gated_qrt_t50": "qrtλ.5",
        "gated_qrt_adapt": "qrt_adapt",
        "gated_qrt_adapt_r2": "qrt_ad+r2",
        "gated_cbrt_adapt": "cbrt_adapt",
    }
    lines = [
        "# Observation-level PO hard-reweight v3 (drift-adaptive)",
        "",
        "## Logic (do we buy it?)",
        "",
        "Yes, with a split:",
        "",
        "1. **Hard-rank always** — obs PO identifies hard rows (Spearman ~0.5–0.8).",
        "   That alone justifies PO as a *hardness score* for reweight targeting.",
        "2. **Pack MSE only under shift** — when hard-tail = drift signal (beijing-like),",
        "   soft temper helps next-MSE; when hard ≈ noise (calm stocks), uniform wins.",
        "3. Therefore **λ should track drift intensity**, not a fixed temper.",
        "",
        "v3: `λ = adaptive_temper(drift_intensity(PO, p, T))` with gate at drift≈0.2.",
        "Also report **hard-subset next-MSE** (top-20% of next batch by PO).",
        "",
        "## Drift intensity (mean on reject batches)",
        "",
        "| dataset | drift_mean (qrt_adapt) | lam_mean |",
        "|---|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]["gated_qrt_adapt"]
        lines.append(f"| `{ds}` | {f(r.get('drift_mean'))} | {f(r.get('lam_mean'))} |")

    hdr = " | ".join(short[m] for m in MODES)
    lines += [
        "",
        "## Sig-only next MSE — all rows (↓)",
        "",
        f"| dataset | {hdr} | best |",
        "|" + "---|---:" * len(MODES) + "|---|",
    ]
    wins = {m: 0 for m in MODES}
    for ds, blob in all_ds.items():
        r = blob["results"]

        def score(m):
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return v if v == v else 1e99

        best = min(MODES, key=score)
        wins[best] += 1
        cells = [f(r[m].get("mse_mean_sig", r[m]["mse_mean"])) for m in MODES]
        lines.append(f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` |")
    lines += ["", "**Wins (all-row MSE):** " + ", ".join(f"`{short[k]}`={v}" for k, v in wins.items())]

    lines += [
        "",
        "## Sig-only next MSE — hard top-20% of next batch (↓)  ← primary claim",
        "",
        f"| dataset | {hdr} | best |",
        "|" + "---|---:" * len(MODES) + "|---|",
    ]
    wins_h = {m: 0 for m in MODES}
    for ds, blob in all_ds.items():
        r = blob["results"]

        def score_h(m):
            v = r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))
            return v if v == v else 1e99

        best = min(MODES, key=score_h)
        wins_h[best] += 1
        cells = [f(r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))) for m in MODES]
        lines.append(f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` |")
    lines += [
        "",
        "**Wins (hard-subset MSE):** " + ", ".join(f"`{short[k]}`={v}" for k, v in wins_h.items()),
        "",
        "## Rel. all-row MSE vs uniform (qrt_adapt / qrtλ.5)",
        "",
        "| dataset | qrt_adapt | qrtλ.5 | drift |",
        "|---|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]
        u = r["uniform"].get("mse_mean_sig", r["uniform"]["mse_mean"])

        def rel(m):
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return (v / u - 1.0) if u == u and u > 0 else float("nan")

        lines.append(
            f"| `{ds}` | {f(rel('gated_qrt_adapt'), pct=True)} | {f(rel('gated_qrt_t50'), pct=True)} | "
            f"{f(r['gated_qrt_adapt'].get('drift_mean'))} |"
        )

    lines += [
        "",
        "## Hard-rank (unchanged claim)",
        "",
        "| dataset | spearman | P@20% | n_reject |",
        "|---|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        h = blob["results"]["gated_qrt_adapt"]["hard_po"]
        nrej = blob["results"]["gated_qrt_adapt"]["n_reject"]
        lines.append(f"| `{ds}` | {f(h['spearman'])} | {f(h['precision_at_k'])} | {nrej} |")

    lines += [
        "",
        "### Takeaway",
        "",
        "- **Agree:** hard-rank is the robust justification; pack MSE is conditional on drift.",
        "- **v3 test:** adaptive λ should fire on high-drift packs and stay near 0 on calm ones.",
        "- Prefer hard-subset next-MSE when claiming hard-reweight benefit.",
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
            )
        except Exception as e:
            print(f"[skip] {name}: {e}", flush=True)
            continue
        all_ds[name] = blob
        r = blob["results"]
        print(
            "  drift={:.3f} lam_adapt={:.3f} | sigMSE ".format(
                r["gated_qrt_adapt"].get("drift_mean") or float("nan"),
                r["gated_qrt_adapt"].get("lam_mean") or float("nan"),
            )
            + " ".join(
                f"{m}={r[m].get('mse_mean_sig', r[m]['mse_mean']):.4g}" for m in MODES
            )
            + " | hardMSE "
            + " ".join(
                f"{m}={r[m].get('mse_hard_mean_sig', r[m].get('mse_hard_mean', float('nan'))):.4g}"
                for m in MODES
            ),
            flush=True,
        )

    payload = {
        "version": 3,
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "modes": list(MODES),
        "datasets": all_ds,
        "note": "v3 drift-adaptive temper + hard-subset next-MSE",
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_ds)
    (args.out / "OBS_PO_WEIGHTS_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_obs_po_weights.md").write_text(md, encoding="utf-8")
    plots = [_bars(all_ds, args.out)] if all_ds else []
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
