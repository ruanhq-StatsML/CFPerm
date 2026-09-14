#!/usr/bin/env python3
"""Compare PO power weights: PO^{1/2} vs PO^{1/3} (gated T0/T1 re-fit).

Default remains uniform. On OnlineRFPerm reject:
  gated_sqrt : w ∝ PO^{1/2}
  gated_cbrt : w ∝ PO^{1/3}   ← softer adaptation

  PYTHONPATH=. python3 scripts/run_agod_po_cbrt_compare.py \\
    --datasets metro_interstate beijing_pm25 stocks_AAPL waymo_proxy
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error

from agod.online_rfperm import fit_online_rfperm, update_online_rfperm
from agod.po_iptw import dre_weights, instance_po_risk, po_iptw_weights
from agod.po_refit import current_batch_cbrt_po_weights, current_batch_sqrt_po_weights
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

MODES = ("uniform", "sqrt", "cbrt", "gated_sqrt", "gated_cbrt", "dre")
COLORS = {
    "uniform": "#4C566A",
    "sqrt": "#A3BE8C",
    "cbrt": "#EBCB8B",
    "gated_sqrt": "#88C0D0",
    "gated_cbrt": "#5E81AC",
    "dre": "#BF616A",
}


def ensure_loaders() -> None:
    for t in ("MSFT", "IWM", "AAPL"):
        key = f"stocks_{t}"
        if key not in LOADERS:
            LOADERS[key] = (lambda tk: (lambda root, max_n=20000: load_stocks(root, tk, max_n)))(t)


def pca_fit(X, d, seed):
    if X.shape[1] <= d:
        return X.astype(np.float32)
    return PCA(n_components=d, random_state=seed).fit_transform(X).astype(np.float32)


def fit_rf(X, y, w, seed):
    m = RandomForestRegressor(
        n_estimators=40, max_depth=8, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    m.fit(X, y, sample_weight=w)
    return m


def batch_po(X0, y0, X1, y1, seed):
    if len(X0) < 16 or len(X1) < 16:
        return 0.0
    rf = RandomForestRegressor(
        n_estimators=20, max_depth=4, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    rf.fit(X0, y0)
    e0 = float(np.mean(np.abs(y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(y1 - rf.predict(X1))))
    return max(e1 - e0, 0.0)


def make_stream(X, y, bs, n_batches):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, len(X), bs)]


def run_mode(stream, mode: str, seed: int, *, n_burn: int, alpha: float, n_control: int) -> dict:
    mse_next, mae_next, gate_on, posthoc = [], [], [], []
    X0, y0 = stream[0]
    probe = fit_rf(X0, y0, np.ones(len(y0)), seed)
    rfperm = fit_online_rfperm(X0, y0, seed=seed)

    for t in range(1, len(stream)):
        Xp, yp = stream[t - 1]
        Xc, yc = stream[t]
        burn = t <= n_burn
        step = update_online_rfperm(
            rfperm, Xc, yc, burn_in=burn, alpha=alpha, ewma=True, fdr="alpha_investing"
        )
        rejected = bool(step["reject"])
        gate_on.append(int(rejected))
        used = 0

        if mode == "dre":
            w = dre_weights(Xp, Xc, seed=seed + t)
        elif mode == "uniform":
            w = np.ones(len(yc))
        elif mode in ("sqrt", "cbrt"):
            po_b = batch_po(Xp, yp, Xc, yc, seed + t)
            pred = probe.predict(Xc)
            po_row = instance_po_risk(yc, pred, batch_po=po_b, mix=0.5)
            w = po_iptw_weights(po_row, mode=mode)  # type: ignore[arg-type]
        elif mode == "gated_sqrt":
            if rejected:
                w = current_batch_sqrt_po_weights(
                    stream, t, seed=seed + 7 * t, n_control=n_control
                )
                used = 1
            else:
                w = np.ones(len(yc))
        elif mode == "gated_cbrt":
            if rejected:
                w = current_batch_cbrt_po_weights(
                    stream, t, seed=seed + 7 * t, n_control=n_control
                )
                used = 1
            else:
                w = np.ones(len(yc))
        else:
            raise ValueError(mode)
        posthoc.append(used)

        model = fit_rf(Xc, yc, w, seed + 17 * t)
        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            pn = model.predict(Xn)
            mse_next.append(float(mean_squared_error(yn, pn)))
            mae_next.append(float(mean_absolute_error(yn, pn)))
        probe = fit_rf(Xc, yc, w, seed + 31 * t)

    mse = np.asarray(mse_next, float)
    mae = np.asarray(mae_next, float)
    return {
        "mode": mode,
        "mse_next": mse_next,
        "mae_next": mae_next,
        "mse_mean": float(mse.mean()) if len(mse) else float("nan"),
        "mse_std": float(mse.std()) if len(mse) else float("nan"),
        "mae_mean": float(mae.mean()) if len(mae) else float("nan"),
        "cum_mse": float(mse.sum()) if len(mse) else float("nan"),
        "gate_duty": float(np.mean(gate_on)) if gate_on else 0.0,
        "posthoc_duty": float(np.mean(posthoc)) if posthoc else 0.0,
        "gate_on": gate_on,
        "posthoc_on": posthoc,
    }


def regrets(results: Dict[str, dict]) -> Dict[str, dict]:
    base_u = np.asarray(results["uniform"]["mse_next"], float)
    base_d = np.asarray(results["dre"]["mse_next"], float)
    out = {}
    for mode, r in results.items():
        m = np.asarray(r["mse_next"], float)
        ex_u = m - base_u[: len(m)]
        ex_d = m - base_d[: len(m)]
        out[mode] = {
            "cum_regret_vs_uniform": float(ex_u.sum()),
            "cum_regret_vs_dre": float(ex_d.sum()),
            "rel_mse_vs_uniform": float(m.mean() / (base_u.mean() + 1e-12)),
            "win_rate_vs_dre": float(np.mean(ex_d < 0)),
        }
    return out


def plot_all(all_res: dict, out_dir: Path) -> List[Path]:
    paths = []
    ds = list(all_res.keys())
    fig, ax = plt.subplots(figsize=(max(8, 1.7 * len(ds)), 4.6))
    x = np.arange(len(ds))
    w = 0.13
    for i, mode in enumerate(MODES):
        vals = [
            all_res[d]["results"][mode].get(
                "rel_mse_vs_uniform_sig",
                all_res[d]["regret"][mode]["rel_mse_vs_uniform"],
            )
            for d in ds
        ]
        ax.bar(x + (i - 2.5) * w, vals, w, label=mode, color=COLORS[mode])
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("Next-MSE / uniform (significant only)")
    ax.set_title("PO^{1/2} vs PO^{1/3} — significant batches only")
    ax.legend(ncol=3, fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "mse_rel_cbrt_bars.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)
    return paths


def report(all_res: dict) -> str:
    lines = [
        "# PO^{1/3} vs PO^{1/2} adaptation (OnlineRFPerm-gated)",
        "",
        "Default = **uniform**. On significant reject, T0/T1 re-fit then:",
        "`gated_sqrt`: w∝PO^{1/2}; `gated_cbrt`: w∝PO^{1/3} (softer).",
        "",
        "Primary metric = next-MSE on **significant (reject) batches only** —",
        "non-reject steps stay uniform and are excluded.",
        "",
        "| dataset | n_sig | unif | always-√ | always-∛ | gated-√ | **gated-∛** | dre | best |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    wins = {m: 0 for m in MODES}
    cbrt_beats_sqrt = 0
    n = 0
    for d, blob in all_res.items():
        res = blob["results"]
        means = {m: res[m].get("mse_mean_sig", res[m]["mse_mean"]) for m in MODES}
        best = min(means, key=means.get)
        wins[best] += 1
        n += 1
        if means["gated_cbrt"] <= means["gated_sqrt"]:
            cbrt_beats_sqrt += 1
        n_sig = res["gated_cbrt"].get("n_significant", "?")
        lines.append(
            f"| `{d}` | {n_sig} | {means['uniform']:.4g} | {means['sqrt']:.4g} | "
            f"{means['cbrt']:.4g} | {means['gated_sqrt']:.4g} | "
            f"**{means['gated_cbrt']:.4g}** | {means['dre']:.4g} | `{best}` |"
        )
    lines += [
        "",
        f"**Wins (sig-only):** " + ", ".join(f"`{m}`={wins[m]}" for m in MODES),
        f"**gated-∛ ≤ gated-√ (sig-only):** `{cbrt_beats_sqrt}/{n}`",
        "",
        "```python",
        "w = po ** (1/3)   # softer than sqrt; closer to uniform",
        "metric = mean(mse_next[reject])  # drop non-significant",
        "```",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ensure_loaders()
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=[
            "metro_interstate",
            "beijing_pm25",
            "stocks_AAPL",
            "waymo_proxy",
            "stocks_MSFT",
            "stocks_IWM",
        ],
    )
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--n-burn", type=int, default=5)
    ap.add_argument("--n-control", type=int, default=1)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--max-n", type=int, default=8000)
    ap.add_argument("--pca-d", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_po_cbrt"))
    args = ap.parse_args()

    need = args.batch_size * args.n_batches
    all_res: Dict[str, dict] = {}
    skipped: List[str] = []

    for name in args.datasets:
        loader = LOADERS.get(name)
        if loader is None:
            print(f"[skip] unknown {name}")
            skipped.append(name)
            continue
        packed = loader(args.root, max_n=max(need, args.max_n))
        if packed is None or len(packed[0]) < need:
            print(f"[skip] {name}")
            skipped.append(name)
            continue
        X, y, meta = packed
        X = pca_fit(X[:need], args.pca_d, args.seed)
        y = y[:need]
        stream = make_stream(X, y, args.batch_size, args.n_batches)
        print(f"=== {name} n={len(X)} d={X.shape[1]} batches={len(stream)} ===", flush=True)
        results = {}
        for mode in MODES:
            print(f"  [{mode}]", flush=True)
            results[mode] = run_mode(
                stream,
                mode,
                args.seed,
                n_burn=args.n_burn,
                alpha=args.alpha,
                n_control=args.n_control,
            )
            print(
                f"    mse={results[mode]['mse_mean']:.6g} "
                f"posthoc={results[mode]['posthoc_duty']:.2f}"
            )
        annotate_results_with_sig(results, preferred_gate_modes=("gated_cbrt", "gated_sqrt"))
        all_res[name] = {"results": results, "regret": regrets(results), "meta": meta}

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "modes": list(MODES),
        "skipped": skipped,
        "datasets": all_res,
        "note": (
            "cbrt = PO**(1/3); softer than sqrt; default still uniform. "
            "Primary compare uses mse_mean_sig (reject batches only)."
        ),
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_res)
    (args.out / "PO_CBRT_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_po_cbrt.md").write_text(md, encoding="utf-8")
    plots = plot_all(all_res, args.out) if all_res else []
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
