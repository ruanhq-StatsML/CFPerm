#!/usr/bin/env python3
"""Sample-chunk adjacent boards across datasets (business FS dashboard).

Nothing weird: sort by time (or index), cut every N rows (1000 / 2000),
adjacent chunks t→t+1 as train→test, same SelectKBest→HGB smoke, joint board.

Datasets:
  - diffusiondb       : prompt TF-IDF → image_nsfw
  - tencent_gr        : edge feature_grid → y_convert (by e_last_ts)
  - waymo_proxy       : tabular X → y (by row order)
  - metro_interstate  : traffic volume (hourly, sorted)
  - beijing_pm25      : PM2.5 (hourly, sorted)

  PYTHONPATH=. python3 scripts/run_sample_chunk_adjacent_board.py \\
    --chunk-sizes 1000,2000 --out results/sample_chunk_adjacent_board
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from agod.stream_packs import load_beijing_pm25, load_metro_interstate
from scripts.run_diffusiondb_temporal_fsds import (
    build_light_token_X,
    load_subset as load_diffusion_subset,
)

ROOT = Path(__file__).resolve().parents[1]


def chunk_by_n(n_rows: int, chunk_size: int) -> np.ndarray:
    """Integer chunk id for each row after sort (every ``chunk_size`` samples)."""
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0")
    return (np.arange(n_rows) // int(chunk_size)).astype(int)


def adjacent_chunk_pairs(T: np.ndarray) -> List[Tuple[int, int]]:
    occ = sorted(int(t) for t in np.unique(T))
    return list(zip(occ[:-1], occ[1:]))


def fs_adjacent(
    X: np.ndarray,
    y: np.ndarray,
    T: np.ndarray,
    feat_names: List[str],
    *,
    select_k: int,
    y_quantile: float,
    seed: int,
    max_pairs: int,
    min_n: int = 80,
) -> List[Dict[str, Any]]:
    pairs = adjacent_chunk_pairs(T)
    if max_pairs > 0 and len(pairs) > max_pairs:
        idx = np.linspace(0, len(pairs) - 1, num=max_pairs, dtype=int)
        pairs = [pairs[i] for i in idx]
    rows = []
    for t0, t1 in pairs:
        m0, m1 = T == t0, T == t1
        if int(m0.sum()) < min_n or int(m1.sum()) < min_n:
            continue
        y0, y1 = y[m0], y[m1]
        # binary labels stay binary; continuous → early-quantile binary
        uniq0 = set(np.unique(y0[~np.isnan(y0)]).tolist()) if len(y0) else set()
        if uniq0 and uniq0 <= {0.0, 1.0}:
            yb0, yb1 = y0.astype(int), y1.astype(int)
            thr = 0.5
        else:
            thr = float(np.quantile(y0, y_quantile))
            yb0 = (y0 >= thr).astype(int)
            yb1 = (y1 >= thr).astype(int)
        if len(np.unique(yb0)) < 2 or len(np.unique(yb1)) < 2:
            continue
        k = min(select_k, X.shape[1], max(1, int(m0.sum()) - 1))
        try:
            pre = Pipeline(
                [("sc", StandardScaler()), ("var", VarianceThreshold(1e-10))]
            )
            Xv = pre.fit_transform(X[m0], yb0)
            k_eff = min(k, Xv.shape[1])
            sel = SelectKBest(f_classif, k=k_eff)
            Xt = sel.fit_transform(Xv, yb0)
            Xte = sel.transform(pre.transform(X[m1]))
        except Exception:
            continue
        hgb = HistGradientBoostingClassifier(
            max_depth=3, learning_rate=0.1, max_iter=60, random_state=seed
        )
        hgb.fit(Xt, yb0)
        proba = hgb.predict_proba(Xte)[:, 1]
        auc = float(roc_auc_score(yb1, proba))
        var_mask = pre.named_steps["var"].get_support()
        cols_var = [c for c, m in zip(feat_names, var_mask) if m]
        ranking = (
            pd.DataFrame({"feature": cols_var, "f_score": sel.scores_})
            .sort_values("f_score", ascending=False)
            .head(5)
        )
        # cmean top on continuous / binary X
        dmu = X[m1].mean(0) - X[m0].mean(0)
        top_cmean = [feat_names[i] for i in np.argsort(-np.abs(dmu))[:5]]
        rows.append(
            {
                "t0": int(t0),
                "t1": int(t1),
                "n0": int(m0.sum()),
                "n1": int(m1.sum()),
                "y_mean_0": float(np.mean(y0)),
                "y_mean_1": float(np.mean(y1)),
                "delta_Y": float(np.mean(y1) - np.mean(y0)),
                "y_threshold": thr,
                "hgb_auc": auc,
                "top_fsds": ranking["feature"].tolist(),
                "top_cmean": top_cmean,
            }
        )
    return rows


def load_diffusion(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    meta = ROOT / "data/diffusiondb/metadata.parquet"
    df = load_diffusion_subset(meta, n_sample=n_sample, seed=seed)
    df = df.sort_values("timestamp").reset_index(drop=True)
    X, names, _ = build_light_token_X(
        df["prompt_clean"].tolist(), method="tfidf", max_features=256, seed=seed
    )
    y = df["image_nsfw"].to_numpy(float)
    return X, y, names, {"dataset": "diffusiondb", "y": "image_nsfw", "sort": "timestamp", "n": len(df)}


def load_tencent(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    path = ROOT / "results/tencent_gr_tabular_ui/feature_grid.parquet"
    df = pd.read_parquet(path)
    # drop ids; keep numeric feats
    y_col = "y_convert" if "y_convert" in df.columns else "e_ctr"
    ts_col = "e_last_ts" if "e_last_ts" in df.columns else None
    if ts_col:
        df = df.sort_values(ts_col).reset_index(drop=True)
    if n_sample < len(df):
        # evenly subsample along time then re-sort
        idx = np.linspace(0, len(df) - 1, num=n_sample, dtype=int)
        df = df.iloc[idx].reset_index(drop=True)
    drop = {
        "user_id",
        "item_id",
        y_col,
        "e_last_ts",
        "last_ts",
        # direct convert leakage for adjacent FS board
        "e_n_cnv",
        "u_n_cnv",
        "u_has_convert",
        "u_log1p_n_cnv",
        "u_cvr",
        "u_ctcvr",
        "u_user_ctcvr_rank",
        "i_n_cnv",
        "i_n_as_convert_terminal",
        "i_log1p_n_cnv",
        "i_item_cnv_rank",
        "i_cvr",
        "i_ctcvr",
    }
    feat_cols = [
        c
        for c in df.columns
        if c not in drop and pd.api.types.is_numeric_dtype(df[c])
    ]
    X = df[feat_cols].to_numpy(dtype=float)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    y = df[y_col].to_numpy(float)
    return X, y, feat_cols, {
        "dataset": "tencent_gr_edges",
        "y": y_col,
        "sort": ts_col or "row",
        "n": len(df),
        "d": len(feat_cols),
    }


def load_waymo(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    blob = np.load(ROOT / "data/stream_packs/waymo_proxy/waymo_proxy_xy.npz")
    X = np.asarray(blob["X"], dtype=float)
    y = np.asarray(blob["y"], dtype=float).ravel()
    if n_sample < len(X):
        X, y = X[:n_sample], y[:n_sample]
    names = [f"x{j}" for j in range(X.shape[1])]
    return X, y, names, {
        "dataset": "waymo_proxy",
        "y": "proxy_y",
        "sort": "row_index",
        "n": len(y),
        "d": X.shape[1],
    }


def load_metro(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    del seed  # time-ordered pack; seed unused
    X, y, meta = load_metro_interstate(ROOT, max_n=n_sample)
    names = [f"m{j}" for j in range(X.shape[1])]
    return X.astype(float), y.astype(float), names, {
        "dataset": "metro_interstate",
        "y": meta.get("target", "traffic_volume"),
        "sort": "date_time",
        "n": int(meta["n"]),
        "d": int(meta["d"]),
    }


def load_beijing(n_sample: int, seed: int) -> Tuple[np.ndarray, np.ndarray, List[str], Dict]:
    del seed
    X, y, meta = load_beijing_pm25(ROOT, max_n=n_sample)
    names = [f"p{j}" for j in range(X.shape[1])]
    return X.astype(float), y.astype(float), names, {
        "dataset": "beijing_pm25",
        "y": meta.get("target", "pm2.5"),
        "sort": "stamp",
        "n": int(meta["n"]),
        "d": int(meta["d"]),
    }


LOADERS = {
    "diffusiondb": load_diffusion,
    "tencent_gr": load_tencent,
    "waymo_proxy": load_waymo,
    "metro_interstate": load_metro,
    "beijing_pm25": load_beijing,
}


def plot_dataset_board(name: str, by_chunk: Dict[int, List[Dict]], out_png: Path) -> None:
    sizes = sorted(by_chunk.keys())
    fig, axes = plt.subplots(1, len(sizes), figsize=(4.0 * len(sizes), 3.4), squeeze=False)
    for ax, cs in zip(axes[0], sizes):
        rows = by_chunk[cs]
        if not rows:
            ax.set_title(f"chunk={cs} (empty)")
            ax.axis("off")
            continue
        xs = np.arange(len(rows))
        dys = [r["delta_Y"] for r in rows]
        aucs = [r["hgb_auc"] for r in rows]
        ax.bar(xs - 0.15, dys, width=0.3, color="#4C78A8", label="ΔȲ")
        ax2 = ax.twinx()
        ax2.plot(xs + 0.15, aucs, "o-", color="#F58518", ms=4, label="HGB AUC")
        ax.axhline(0, color="gray", lw=0.7)
        ax.set_title(f"every {cs} samples")
        ax.set_xlabel("adjacent chunk idx")
        ax.set_ylabel("ΔȲ")
        ax2.set_ylabel("AUC")
    fig.suptitle(f"{name}: sample-chunk adjacent board", y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--datasets",
        type=str,
        default="diffusiondb,tencent_gr,waymo_proxy,metro_interstate,beijing_pm25",
    )
    ap.add_argument("--chunk-sizes", type=str, default="1000,2000")
    ap.add_argument("--n-sample", type=int, default=8000, help="cap rows per dataset")
    ap.add_argument("--select-k", type=int, default=20)
    ap.add_argument("--y-quantile", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-pairs", type=int, default=12)
    ap.add_argument(
        "--out", type=Path, default=ROOT / "results/sample_chunk_adjacent_board"
    )
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    chunk_sizes = [int(x) for x in args.chunk_sizes.split(",") if x.strip()]
    datasets = [x.strip() for x in args.datasets.split(",") if x.strip()]
    t0 = time.time()
    report: Dict[str, Any] = {
        "note": (
            "Sort → every N samples as a window → adjacent chunk FS board. "
            "Business grain = sample count (1000/2000); joint viz only."
        ),
        "chunk_sizes": chunk_sizes,
        "datasets": {},
    }

    md = [
        "# Sample-chunk adjacent boards (multi-dataset)",
        "",
        report["note"],
        "",
    ]

    for ds in datasets:
        if ds not in LOADERS:
            print("skip unknown", ds)
            continue
        print("load", ds, flush=True)
        try:
            X, y, names, meta = LOADERS[ds](args.n_sample, args.seed)
        except Exception as e:
            report["datasets"][ds] = {"ok": False, "error": str(e)}
            continue
        by_chunk: Dict[int, List[Dict]] = {}
        for cs in chunk_sizes:
            if cs >= len(y):
                by_chunk[cs] = []
                continue
            T = chunk_by_n(len(y), cs)
            rows = fs_adjacent(
                X,
                y,
                T,
                names,
                select_k=args.select_k,
                y_quantile=args.y_quantile,
                seed=args.seed,
                max_pairs=args.max_pairs,
            )
            by_chunk[cs] = rows
            print(f"  chunk={cs} pairs_reported={len(rows)}", flush=True)

        png = args.out / f"{ds}_board.png"
        plot_dataset_board(ds, by_chunk, png)
        report["datasets"][ds] = {
            "ok": True,
            "meta": meta,
            "by_chunk": {
                str(cs): {
                    "n_pairs": len(rows),
                    "mean_delta_Y": float(np.mean([r["delta_Y"] for r in rows]))
                    if rows
                    else None,
                    "mean_auc": float(np.mean([r["hgb_auc"] for r in rows]))
                    if rows
                    else None,
                    "rows": rows,
                }
                for cs, rows in by_chunk.items()
            },
            "board_png": str(png.name),
        }
        md += [f"## {ds}", f"- meta: `{meta}`", f"- ![{ds}]({png.name})", ""]
        for cs, rows in by_chunk.items():
            md += [
                f"### every {cs} samples",
                f"- pairs={len(rows)}",
                "",
                "| t0→t1 | n0/n1 | ΔȲ | AUC | top_fsds |",
                "|---|---|---:|---:|---|",
            ]
            for r in rows[:8]:
                md.append(
                    f"| {r['t0']}→{r['t1']} | {r['n0']}/{r['n1']} | "
                    f"{r['delta_Y']:.4f} | {r['hgb_auc']:.3f} | "
                    f"{', '.join(r['top_fsds'][:3])} |"
                )
            md.append("")

    report["sec"] = float(time.time() - t0)
    (args.out / "summary.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n"
    )
    (args.out / "SAMPLE_CHUNK_BOARD.md").write_text("\n".join(md) + "\n")
    print(
        json.dumps(
            {
                "datasets": list(report["datasets"].keys()),
                "chunk_sizes": chunk_sizes,
                "out": str(args.out),
                "sec": report["sec"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
