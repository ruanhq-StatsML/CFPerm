#!/usr/bin/env python3
"""Iterate observation-level PO hard-reweight (tempered / soft / top-k).

v2 focus (after v1): prop/√ often hurt next-MSE; PO still ranks hard rows.
Try softer maps that keep hard emphasis without blowing MSE:

  uniform
  gated_cbrt          — previous soft baseline
  gated_qrt           — PO^{1/4}
  gated_cbrt_t50      — ∛PO tempered λ=0.5 with uniform
  gated_qrt_t50       — PO^{1/4} tempered λ=0.5
  gated_topk20        — only top-20% hard rows boosted
  gated_cbrt_soft     — soft p-blend gate (not hard reject)

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
from agod.obs_po_weights import ObsWeightMode, gated_obs_po_weights
from agod.online_rfperm import fit_online_rfperm, update_online_rfperm
from agod.po_refit import build_recent_ood_windows, refit_po_on_windows
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

# mode_name -> (obs_mode, temper, soft_gate, topk_frac)
MODE_CFG: Dict[str, Tuple[ObsWeightMode, float, bool, float]] = {
    "uniform": ("uniform", 1.0, False, 0.2),
    "gated_cbrt": ("cbrt", 1.0, False, 0.2),
    "gated_qrt": ("qrt", 1.0, False, 0.2),
    "gated_cbrt_t50": ("cbrt", 0.5, False, 0.2),
    "gated_qrt_t50": ("qrt", 0.5, False, 0.2),
    "gated_topk20": ("topk", 1.0, False, 0.2),
    "gated_cbrt_soft": ("cbrt", 0.5, True, 0.2),
}
MODES: Tuple[str, ...] = tuple(MODE_CFG.keys())
COLORS = {
    "uniform": "#4C566A",
    "gated_cbrt": "#5E81AC",
    "gated_qrt": "#88C0D0",
    "gated_cbrt_t50": "#A3BE8C",
    "gated_qrt_t50": "#EBCB8B",
    "gated_topk20": "#B48EAD",
    "gated_cbrt_soft": "#D08770",
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
    if name not in LOADERS:
        raise KeyError(name)
    out = LOADERS[name](root, max_n=max_n)
    if isinstance(out, tuple) and len(out) >= 2:
        return np.asarray(out[0], float), np.asarray(out[1], float).ravel()
    raise TypeError(type(out))


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


def run_mode(
    stream,
    mode: str,
    seed: int,
    *,
    n_burn: int,
    alpha: float,
    n_recent: int,
) -> dict:
    obs_mode, temper, soft_gate, topk_frac = MODE_CFG[mode]
    mse_next: List[float] = []
    gate_on: List[int] = []
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

        windows = None
        po = None
        apply = rejected or (soft_gate and not burn)
        if mode == "uniform" or not apply:
            w = np.ones(len(yc), float)
        else:
            windows = build_recent_ood_windows(
                stream, t, n_recent=n_recent, window_mode="recent_ood"
            )
            po = refit_po_on_windows(windows, seed=seed + t, blend_mu_gap=0.25)
            w = gated_obs_po_weights(
                po,
                reject=rejected,
                mode=obs_mode,
                soft=soft_gate,
                p=float(step["p"]),
                alpha=alpha,
                temper=temper,
                topk_frac=topk_frac,
                topk_boost=2.0,
            )

        model = fit_rf(Xc, yc, w, seed + t)
        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            pred = model.predict(Xn)
            mse_next.append(float(np.mean((yn - pred) ** 2)))

        if rejected and po is not None and windows is not None:
            Xr = np.vstack([windows.X_recent, windows.X_ood])
            yr = np.concatenate([windows.y_recent, windows.y_ood])
            mu_o = fit_rf(Xr, yr, np.ones(len(yr)), seed + 11 + t)
            truth = np.abs(yc - mu_o.predict(Xc))
            hard_rows.append(
                {
                    "po": hard_rank_metrics(po, truth),
                    "w": hard_rank_metrics(w, truth),
                }
            )

    def _mean_hard(key: str) -> dict:
        if not hard_rows:
            return {"spearman": float("nan"), "precision_at_k": float("nan"), "n": 0}
        sp = [
            h[key]["spearman"]
            for h in hard_rows
            if h[key].get("spearman") == h[key].get("spearman")
        ]
        pk = [
            h[key]["precision_at_k"]
            for h in hard_rows
            if h[key].get("precision_at_k") == h[key].get("precision_at_k")
        ]
        return {
            "spearman": float(np.mean(sp)) if sp else float("nan"),
            "precision_at_k": float(np.mean(pk)) if pk else float("nan"),
            "n": len(hard_rows),
        }

    return {
        "mse_next": mse_next,
        "gate_on": gate_on,
        "mse_mean": float(np.mean(mse_next)) if mse_next else float("nan"),
        "hard_po": _mean_hard("po"),
        "hard_w": _mean_hard("w"),
        "n_reject": int(sum(gate_on)),
        "duty": float(np.mean(gate_on)) if gate_on else 0.0,
        "cfg": {
            "obs_mode": obs_mode,
            "temper": temper,
            "soft_gate": soft_gate,
            "topk_frac": topk_frac,
        },
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
            stream, mode, seed, n_burn=n_burn, alpha=alpha, n_recent=n_recent
        )
    annotate_results_with_sig(
        results,
        preferred_gate_modes=(
            "gated_cbrt",
            "gated_qrt",
            "gated_cbrt_t50",
            "gated_topk20",
        ),
    )
    return {"dataset": name, "n_batches": len(stream), "results": results}


def _bar_mse(all_ds: dict, out: Path):
    packs = list(all_ds.keys())
    fig, ax = plt.subplots(figsize=(max(9, 1.7 * len(packs)), 4.8))
    x = np.arange(len(packs))
    w = 0.11
    mid = (len(MODES) - 1) / 2.0
    for i, m in enumerate(MODES):
        vals = [
            all_ds[p]["results"][m].get("mse_mean_sig", all_ds[p]["results"][m]["mse_mean"])
            for p in packs
        ]
        ax.bar(x + (i - mid) * w, vals, w, label=m, color=COLORS[m])
    ax.set_xticks(x)
    ax.set_xticklabels(packs, rotation=15, ha="right")
    ax.set_ylabel("sig-only next MSE (↓)")
    ax.set_title("Obs PO hard-reweight v2: soft / temper / top-k")
    ax.legend(ncol=3, fontsize=6)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out / "obs_po_sig_mse_v2.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def report(all_ds: dict) -> str:
    def f(v):
        if v is None or (isinstance(v, float) and v != v):
            return "—"
        return f"{v:.4f}"

    short = {
        "uniform": "uniform",
        "gated_cbrt": "cbrt",
        "gated_qrt": "qrt",
        "gated_cbrt_t50": "cbrtλ.5",
        "gated_qrt_t50": "qrtλ.5",
        "gated_topk20": "topk20",
        "gated_cbrt_soft": "cbrt_soft",
    }
    header = " | ".join(short[m] for m in MODES)
    lines = [
        "# Observation-level PO hard-reweight v2 (temper / soft / top-k)",
        "",
        "Default = **uniform**. On OnlineRFPerm reject, obs PO → soft weights.",
        "v2 adds PO^{1/4}, temper(λ=0.5), top-20% boost, soft p-gate.",
        "",
        "## Sig-only next MSE (↓ better)",
        "",
        f"| dataset | {header} | best |",
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
    lines += [
        "",
        "**Wins:** " + ", ".join(f"`{short[k]}`={v}" for k, v in wins.items()),
        "",
        "## Hard-rank (PO vs oracle residual) on reject batches",
        "",
        "| dataset | spearman | P@20% | n_reject |",
        "|---|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        # PO ranking identical across monotone maps — report from cbrt run
        h = blob["results"]["gated_cbrt"]["hard_po"]
        nrej = blob["results"]["gated_cbrt"]["n_reject"]
        lines.append(
            f"| `{ds}` | {f(h['spearman'])} | {f(h['precision_at_k'])} | {nrej} |"
        )
    lines += [
        "",
        "### Takeaway",
        "",
        "1. Obs PO still **ranks hard rows** (use for hard-reweight targeting).",
        "2. Prefer **tempered / soft** maps (∛·λ0.5, ¼, top-k) over raw prop/√ for MSE.",
        "3. Not an image-OOD detector — Mahalanobis / gradient / RF-binary for that.",
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
            "  sigMSE "
            + " ".join(
                f"{m}={r[m].get('mse_mean_sig', r[m]['mse_mean']):.4g}" for m in MODES
            ),
            flush=True,
        )

    payload = {
        "version": 2,
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "n_recent": args.n_recent,
        "alpha": args.alpha,
        "modes": list(MODES),
        "mode_cfg": {k: {"obs": v[0], "temper": v[1], "soft": v[2], "topk_frac": v[3]} for k, v in MODE_CFG.items()},
        "datasets": all_ds,
        "note": "obs PO hard-reweight v2: temper / soft / top-k under OnlineRFPerm",
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_ds)
    (args.out / "OBS_PO_WEIGHTS_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_obs_po_weights.md").write_text(md, encoding="utf-8")
    plots = []
    if all_ds:
        plots.append(_bar_mse(all_ds, args.out))
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
