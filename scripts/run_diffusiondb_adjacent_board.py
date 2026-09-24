#!/usr/bin/env python3
"""DiffusionDB adjacent-batch rolling board (business time cuts).

Nothing fancy: pick a calendar width (5 / 10 / 20 / … min), groupby window,
for each adjacent pair (T=t → T=t+1) run the same FS smoke (cmean / VIMP / FSDS),
and dump a joint dashboard table + small PNG.

  PYTHONPATH=. python3 scripts/run_diffusiondb_adjacent_board.py \\
    --n-sample 4000 --widths-min 5,10,20,60 \\
    --out results/diffusiondb_adjacent_board
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.run_diffusiondb_temporal_fsds import (
    build_light_token_X,
    concat_hyperparams,
    load_subset,
    mean_shift_by_feature,
    rf_domain_vimp,
    run_fsds_y,
)

ROOT = Path(__file__).resolve().parents[1]


def _assign_width_minutes(df: pd.DataFrame, width_min: float) -> pd.DataFrame:
    out = df.sort_values("timestamp").reset_index(drop=True).copy()
    ts = pd.to_datetime(out["timestamp"], utc=True)
    t0 = ts.min()
    minutes = (ts - t0).dt.total_seconds() / 60.0
    out["T"] = np.floor(minutes / float(width_min)).astype(int)
    return out


def _adjacent_pairs(T: np.ndarray) -> List[tuple[int, int]]:
    occupied = sorted(int(t) for t in np.unique(T))
    return [(a, b) for a, b in zip(occupied[:-1], occupied[1:]) if b == a + 1]


def run_width_board(
    df: pd.DataFrame,
    *,
    width_min: float,
    max_features: int,
    select_k: int,
    y_quantile: float,
    seed: int,
    concat_hp: bool,
    max_pairs: int,
) -> Dict:
    df_w = _assign_width_minutes(df, width_min)
    X, names, _ = build_light_token_X(
        df_w["prompt_clean"].tolist(),
        method="tfidf",
        max_features=max_features,
        seed=seed,
    )
    if concat_hp:
        X, names = concat_hyperparams(X, names, df_w)
    T = df_w["T"].to_numpy(int)
    Y = df_w["image_nsfw"].to_numpy(float)
    pairs = _adjacent_pairs(T)
    if max_pairs > 0:
        pairs = pairs[:max_pairs]

    rows = []
    for t0, t1 in pairs:
        m0, m1 = T == t0, T == t1
        if m0.sum() < 40 or m1.sum() < 40:
            continue
        y0, y1 = Y[m0], Y[m1]
        delta = float(y1.mean() - y0.mean())
        thr = float(np.quantile(y0, y_quantile))
        yb0 = (y0 >= thr).astype(int)
        yb1 = (y1 >= thr).astype(int)
        if len(np.unique(yb0)) < 2 or len(np.unique(yb1)) < 2:
            continue
        # cheap VIMP on the two-window subset only
        Xs = np.vstack([X[m0], X[m1]])
        Ts = np.concatenate([np.zeros(m0.sum()), np.ones(m1.sum())]).astype(int)
        vimp = rf_domain_vimp(Xs, Ts, seed=seed, n_trees=40)
        top_cov = [names[i] for i in np.argsort(-vimp)[:5]]
        cmean = mean_shift_by_feature(X[m0], X[m1], names)
        top_cmean = cmean.head(5)["feature"].tolist()
        res = run_fsds_y(
            X[m0], yb0, X[m1], yb1, names, select_k=select_k, seed=seed
        )
        auc = None
        top_fs = []
        if res.get("ok") and res.get("models", {}).get("hgb"):
            auc = float(res["models"]["hgb"]["auc"])
            top_fs = [r["feature"] for r in res["ranking"].head(5).to_dict("records")]
        rows.append(
            {
                "t0": int(t0),
                "t1": int(t1),
                "n0": int(m0.sum()),
                "n1": int(m1.sum()),
                "y_mean_0": float(y0.mean()),
                "y_mean_1": float(y1.mean()),
                "delta_Y": delta,
                "hgb_auc": auc,
                "top_cov": top_cov,
                "top_cmean": top_cmean,
                "top_fsds": top_fs,
            }
        )
    return {
        "width_min": float(width_min),
        "n_bins": int(df_w["T"].nunique()),
        "n_adjacent_pairs": len(pairs),
        "n_reported": len(rows),
        "rows": rows,
    }


def _plot_board(boards: List[Dict], out_png: Path) -> None:
    fig, axes = plt.subplots(1, len(boards), figsize=(4.2 * len(boards), 3.6), squeeze=False)
    for ax, b in zip(axes[0], boards):
        rows = b["rows"]
        if not rows:
            ax.set_title(f"{b['width_min']:.0f}min (empty)")
            ax.axis("off")
            continue
        xs = np.arange(len(rows))
        dys = [r["delta_Y"] for r in rows]
        aucs = [r["hgb_auc"] if r["hgb_auc"] is not None else np.nan for r in rows]
        ax.bar(xs - 0.15, dys, width=0.3, color="#4C78A8", label="ΔȲ")
        ax2 = ax.twinx()
        ax2.plot(xs + 0.15, aucs, "o-", color="#F58518", label="HGB AUC")
        ax.axhline(0, color="gray", lw=0.8)
        ax.set_title(f"width={b['width_min']:.0f} min")
        ax.set_xlabel("adjacent pair idx")
        ax.set_ylabel("ΔȲ")
        ax2.set_ylabel("AUC")
        ax.set_ylim(min(-0.05, min(dys) - 0.01), max(0.05, max(dys) + 0.01))
    fig.suptitle("DiffusionDB adjacent-batch board (business time cuts)", y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", type=Path, default=ROOT / "data/diffusiondb/metadata.parquet")
    ap.add_argument("--n-sample", type=int, default=4000)
    ap.add_argument("--widths-min", type=str, default="5,10,20,60")
    ap.add_argument("--max-features", type=int, default=256)
    ap.add_argument("--select-k", type=int, default=20)
    ap.add_argument("--y-quantile", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--concat-hyperparams", action="store_true")
    ap.add_argument("--max-pairs", type=int, default=12, help="cap pairs per width")
    ap.add_argument("--out", type=Path, default=ROOT / "results/diffusiondb_adjacent_board")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    widths = [float(x) for x in args.widths_min.split(",") if x.strip()]
    t0 = time.time()
    df = load_subset(args.meta, n_sample=args.n_sample, seed=args.seed)
    boards = []
    for w in widths:
        print(f"board width={w} min", flush=True)
        boards.append(
            run_width_board(
                df,
                width_min=w,
                max_features=args.max_features,
                select_k=args.select_k,
                y_quantile=args.y_quantile,
                seed=args.seed,
                concat_hp=args.concat_hyperparams,
                max_pairs=args.max_pairs,
            )
        )

    _plot_board(boards, args.out / "adjacent_board.png")
    blob = {
        "note": (
            "Business time cuts → adjacent batch FS → joint board. "
            "Not a new estimator; visualization of ΔȲ + late AUC over rolling pairs."
        ),
        "n_sample": args.n_sample,
        "concat_hyperparams": bool(args.concat_hyperparams),
        "widths_min": widths,
        "boards": boards,
        "sec": float(time.time() - t0),
    }
    (args.out / "summary.json").write_text(json.dumps(blob, indent=2, default=str) + "\n")

    lines = [
        "# DiffusionDB adjacent-batch board",
        "",
        blob["note"],
        "",
        f"![board](adjacent_board.png)",
        "",
    ]
    for b in boards:
        lines += [
            f"## width = {b['width_min']:.0f} min",
            f"- bins={b['n_bins']}, adjacent_pairs={b['n_adjacent_pairs']}, reported={b['n_reported']}",
            "",
            "| t0→t1 | n0/n1 | ΔȲ | HGB AUC | top_fsds |",
            "|---|---|---:|---:|---|",
        ]
        for r in b["rows"][:10]:
            auc = f"{r['hgb_auc']:.3f}" if r["hgb_auc"] is not None else "—"
            lines.append(
                f"| {r['t0']}→{r['t1']} | {r['n0']}/{r['n1']} | {r['delta_Y']:.4f} | "
                f"{auc} | {', '.join(r['top_fsds'][:3])} |"
            )
        lines.append("")
    (args.out / "ADJACENT_BOARD.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"widths": widths, "out": str(args.out), "sec": blob["sec"]}, indent=2))


if __name__ == "__main__":
    main()
