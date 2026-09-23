#!/usr/bin/env python3
"""Streaming PO-OOD vs DRE on stream packs — MSE + regret + plots.

  PYTHONPATH=. python3 scripts/run_agod_po_ood_stream_metrics.py \\
    --batch-size 100 --n-batches 40
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

from agod.po_iptw import dre_weights, instance_po_risk, po_iptw_weights
from agod.stream_packs import LOADERS

MODES = ("uniform", "prop", "sqrt", "inv", "dre")
COLORS = {
    "uniform": "#4C566A",
    "prop": "#D08770",
    "sqrt": "#A3BE8C",
    "inv": "#B48EAD",
    "dre": "#BF616A",
}


def pca_fit(X: np.ndarray, d: int, seed: int) -> np.ndarray:
    if X.shape[1] <= d:
        return X.astype(np.float32)
    return PCA(n_components=d, random_state=seed).fit_transform(X).astype(np.float32)


def batch_po(X0, y0, X1, y1, seed: int) -> float:
    if len(X0) < 16 or len(X1) < 16:
        return 0.0
    rf = RandomForestRegressor(
        n_estimators=20, max_depth=4, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    rf.fit(X0, y0)
    e0 = float(np.mean(np.abs(y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(y1 - rf.predict(X1))))
    return max(e1 - e0, 0.0)


def fit_rf(X, y, w, seed: int) -> RandomForestRegressor:
    m = RandomForestRegressor(
        n_estimators=40, max_depth=8, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    m.fit(X, y, sample_weight=w)
    return m


def make_stream(X, y, bs: int, n_batches: int):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, len(X), bs)]


def dataset_properties(X, y, meta: dict, bs: int) -> dict:
    n_b = len(X) // bs
    means = [X[i * bs : (i + 1) * bs].mean(0) for i in range(n_b)]
    drifts = [float(np.linalg.norm(means[i] - means[i - 1])) for i in range(1, len(means))]
    y_m = [float(y[i * bs : (i + 1) * bs].mean()) for i in range(n_b)]
    y_d = [abs(y_m[i] - y_m[i - 1]) for i in range(1, len(y_m))]
    x_std = float(np.std(X) + 1e-8)
    return {
        **meta,
        "n": int(len(X)),
        "d": int(X.shape[1]),
        "batch_size": bs,
        "n_batches": n_b,
        "y_mean": float(np.mean(y)),
        "y_std": float(np.std(y)),
        "x_mean_drift_avg": float(np.mean(drifts)) if drifts else 0.0,
        "x_mean_drift_p90": float(np.percentile(drifts, 90)) if drifts else 0.0,
        "y_mean_drift_avg": float(np.mean(y_d)) if y_d else 0.0,
        "gradual_shift_score": float(np.mean(drifts) / x_std) if drifts else 0.0,
    }


def run_mode(stream, mode: str, seed: int) -> dict:
    mse_cur: List[float] = []
    mse_next: List[float] = []
    mae_next: List[float] = []
    batch_pos: List[float] = []
    X0, y0 = stream[0]
    probe = fit_rf(X0, y0, np.ones(len(y0)), seed)

    for t in range(1, len(stream)):
        Xp, yp = stream[t - 1]
        Xc, yc = stream[t]
        po_b = batch_po(Xp, yp, Xc, yc, seed + t)
        batch_pos.append(po_b)
        if mode == "dre":
            w = dre_weights(Xp, Xc, seed=seed + t)
        else:
            pred = probe.predict(Xc)
            po_row = instance_po_risk(yc, pred, batch_po=po_b, mix=0.5)
            w = po_iptw_weights(po_row, mode=mode)  # type: ignore[arg-type]
        model = fit_rf(Xc, yc, w, seed + 17 * t)
        mse_cur.append(float(mean_squared_error(yc, model.predict(Xc))))
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
        "batch_po": batch_pos,
        "mse_cur": mse_cur,
        "mse_next": mse_next,
        "mae_next": mae_next,
        "mse_mean": float(mse.mean()) if len(mse) else float("nan"),
        "mse_std": float(mse.std()) if len(mse) else float("nan"),
        "mae_mean": float(mae.mean()) if len(mae) else float("nan"),
        "cum_mse": float(mse.sum()) if len(mse) else float("nan"),
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
            "avg_regret_vs_uniform": float(ex_u.mean()),
            "cum_regret_vs_dre": float(ex_d.sum()),
            "avg_regret_vs_dre": float(ex_d.mean()),
            "win_rate_vs_uniform": float(np.mean(ex_u < 0)),
            "win_rate_vs_dre": float(np.mean(ex_d < 0)),
            "cum_excess_traj_vs_uniform": [float(v) for v in np.cumsum(ex_u)],
            "cum_excess_traj_vs_dre": [float(v) for v in np.cumsum(ex_d)],
        }
    return out


def plot_all(all_res: dict, out_dir: Path) -> List[Path]:
    paths: List[Path] = []
    ds_names = list(all_res.keys())
    if not ds_names:
        return paths

    fig, ax = plt.subplots(figsize=(max(8, 1.5 * len(ds_names)), 4.5))
    x = np.arange(len(ds_names))
    width = 0.15
    for i, mode in enumerate(MODES):
        vals = [all_res[d]["results"][mode]["mse_mean"] for d in ds_names]
        ax.bar(x + (i - 2) * width, vals, width, label=mode, color=COLORS[mode])
    ax.set_xticks(x)
    ax.set_xticklabels(ds_names, rotation=20, ha="right")
    ax.set_ylabel("Next-batch MSE (mean)")
    ax.set_title("PO-OOD vs DRE — next-batch MSE by dataset")
    ax.legend(ncol=5, fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "mse_mean_bars.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for d in ds_names:
        traj = all_res[d]["regret"]["sqrt"]["cum_excess_traj_vs_dre"]
        ax.plot(traj, label=d, linewidth=1.8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Stream step t")
    ax.set_ylabel("Cum excess MSE vs DRE (↓ better for √PO)")
    ax.set_title("√PO cumulative regret vs DRE")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = out_dir / "cum_regret_sqrt_vs_dre.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    for d in ds_names:
        fig, ax = plt.subplots(figsize=(8, 4))
        for mode in MODES:
            ax.plot(
                all_res[d]["results"][mode]["mse_next"],
                label=mode,
                color=COLORS[mode],
                alpha=0.9,
            )
        ax.set_title(f"{d}: next-batch MSE trajectory")
        ax.set_xlabel("Stream step t")
        ax.set_ylabel("MSE")
        ax.legend(fontsize=8, ncol=5)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p = out_dir / f"traj_mse_{d}.png"
        fig.savefig(p, dpi=130)
        plt.close(fig)
        paths.append(p)

    fig, ax = plt.subplots(figsize=(10, 0.7 + 0.45 * len(ds_names)))
    ax.axis("off")
    rows = []
    for d in ds_names:
        r = all_res[d]["regret"]
        rows.append(
            [
                d,
                f"{all_res[d]['results']['sqrt']['mse_mean']:.4g}",
                f"{all_res[d]['results']['dre']['mse_mean']:.4g}",
                f"{r['sqrt']['cum_regret_vs_dre']:.4g}",
                f"{r['sqrt']['win_rate_vs_dre']:.2f}",
                f"{r['sqrt']['cum_regret_vs_uniform']:.4g}",
            ]
        )
    table = ax.table(
        cellText=rows,
        colLabels=["dataset", "√PO MSE", "DRE MSE", "cumR(√−DRE)", "win% vs DRE", "cumR vs unif"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.35)
    ax.set_title("Regret / MSE summary (√PO vs DRE)", pad=18)
    fig.tight_layout()
    p = out_dir / "regret_summary_table.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for d in ds_names:
        prop = all_res[d]["properties"]
        adv = all_res[d]["results"]["dre"]["mse_mean"] - all_res[d]["results"]["sqrt"]["mse_mean"]
        ax.scatter(prop["gradual_shift_score"], adv, s=70)
        ax.annotate(
            d,
            (prop["gradual_shift_score"], adv),
            fontsize=8,
            xytext=(4, 4),
            textcoords="offset points",
        )
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("Gradual-shift score (mean Δμ_x / std_x)")
    ax.set_ylabel("MSE advantage √PO − DRE (↑ √PO better)")
    ax.set_title("√PO advantage vs shift intensity")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = out_dir / "shift_vs_po_advantage.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)
    return paths


def markdown_report(all_res: dict) -> str:
    lines = [
        "# Stream packs: PO-OOD vs DRE — MSE & regret",
        "",
        "Real streams + Waymo kinematics proxy. `batch_size=100`.",
        "Metrics: next-batch **MSE/MAE**, cumulative **regret** vs uniform & DRE.",
        "",
        "| dataset | n | d | shift | √PO MSE | DRE MSE | unif MSE | cumR(√−DRE) | win% vs DRE | best |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    beats = 0
    n = 0
    for d, blob in all_res.items():
        res, reg, prop = blob["results"], blob["regret"], blob["properties"]
        means = {m: res[m]["mse_mean"] for m in MODES}
        best = min(means, key=means.get)
        if means["sqrt"] < means["dre"]:
            beats += 1
        n += 1
        lines.append(
            f"| `{d}` | {prop['n']} | {prop['d']} | {prop['gradual_shift_score']:.4f} | "
            f"{means['sqrt']:.4g} | {means['dre']:.4g} | {means['uniform']:.4g} | "
            f"{reg['sqrt']['cum_regret_vs_dre']:.4g} | {reg['sqrt']['win_rate_vs_dre']:.2f} | `{best}` |"
        )
    lines += [
        "",
        f"**√PO vs DRE head-to-head:** `{beats}/{n}`",
        "",
        "```python",
        "w_i = np.sqrt(PO-risk(X_i, Y_i, T_i=1))",
        "```",
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
            "stocks_SPY",
            "stocks_QQQ",
            "waymo_proxy",
            "affec",
        ],
    )
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--max-n", type=int, default=8000)
    ap.add_argument("--pca-d", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_po_ood_metrics"))
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
        try:
            packed = loader(args.root, max_n=max(need, args.max_n))
        except Exception as e:  # noqa: BLE001
            print(f"[skip] {name}: {e}")
            skipped.append(name)
            continue
        if packed is None:
            print(f"[skip] {name}: missing")
            skipped.append(name)
            continue
        X, y, meta = packed
        if len(X) < need:
            print(f"[skip] {name}: need {need}, got {len(X)}")
            skipped.append(name)
            continue
        X = pca_fit(X[:need], args.pca_d, args.seed)
        y = y[:need]
        props = dataset_properties(X, y, meta, args.batch_size)
        stream = make_stream(X, y, args.batch_size, args.n_batches)
        print(
            f"=== {name} n={len(X)} d={X.shape[1]} batches={len(stream)} "
            f"shift={props['gradual_shift_score']:.4f} ===",
            flush=True,
        )
        results = {}
        for mode in MODES:
            print(f"  [{mode}]", flush=True)
            results[mode] = run_mode(stream, mode, args.seed)
            print(
                f"    mse={results[mode]['mse_mean']:.6g} "
                f"mae={results[mode]['mae_mean']:.6g} "
                f"cum={results[mode]['cum_mse']:.6g}"
            )
        all_res[name] = {
            "results": results,
            "regret": regrets(results),
            "properties": props,
        }

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "modes": list(MODES),
        "skipped": skipped,
        "datasets": all_res,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = markdown_report(all_res)
    (args.out / "PO_OOD_METRICS_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_po_ood_stream_metrics.md").write_text(md, encoding="utf-8")
    plots = plot_all(all_res, args.out)
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
