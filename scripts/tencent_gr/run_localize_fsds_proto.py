#!/usr/bin/env python3
"""TencentGR prototype: Localization → FSDS (nothing else).

Pipeline (locked scope)
-----------------------
1. **Localization** — node → ``(user, item)`` 图谱特征网格 + path/covisit 下钻
   (no NetworkX / ChannelAttribution / graph algos)
2. **FSDS** — variance → MI/F select → HistGradientBoosting / LogReg on
   ``y_convert`` at the edge row

Reuse an existing ``feature_grid.parquet`` if present; otherwise build it.

  PYTHONPATH=. python3 scripts/tencent_gr/run_localize_fsds_proto.py \\
    --root data/tencent_subset --max-users 3000 --select-k 24
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import (
    SelectKBest,
    VarianceThreshold,
    f_classif,
    mutual_info_classif,
)
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]

# Edge-level label is ``y_convert = (e_n_cnv > 0)`` → drop direct leaks.
LEAK_COLS = {
    "e_n_cnv",
    "y_convert",
    "user_id",
    "item_id",
    "e_last_ts",
}


def ensure_localize_grid(
    *,
    root: Path,
    grid_dir: Path,
    max_users: int,
    co_window: int,
    rebuild: bool,
) -> Tuple[pd.DataFrame, Path]:
    grid_path = grid_dir / "feature_grid.parquet"
    if grid_path.is_file() and not rebuild:
        print(f"reuse localize grid {grid_path}", flush=True)
        return pd.read_parquet(grid_path), grid_dir

    print("build localize grid (图谱特征 → (user,item)) ...", flush=True)
    from scripts.tencent_gr import run_tabular_user_item_feats as tab

    users = tab.iter_users(root / "seq", max_users)
    user_df, item_df, edge_df, convert_df, covisit_df, meta = tab.accumulate(
        users, co_window=co_window, top_covisit=10
    )
    user_df, item_df = tab.enrich_spectrum_cols(user_df, item_df)
    grid = tab.build_feature_grid(user_df, item_df, edge_df)
    grid_dir.mkdir(parents=True, exist_ok=True)
    user_df.to_parquet(grid_dir / "user_features.parquet", index=False)
    item_df.to_parquet(grid_dir / "item_features.parquet", index=False)
    edge_df.to_parquet(grid_dir / "edge_ui_features.parquet", index=False)
    convert_df.to_csv(grid_dir / "localize_convert_path.csv", index=False)
    covisit_df.to_parquet(grid_dir / "localize_item_covisit.parquet", index=False)
    grid.to_parquet(grid_path, index=False)
    feat_cols = [c for c in grid.columns if c not in LEAK_COLS]
    meta.update(
        {
            "n_grid_rows": int(len(grid)),
            "n_feature_cols": len(feat_cols),
            "feature_cols": feat_cols,
            "label": "y_convert",
        }
    )
    (grid_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    (grid_dir / "feature_cols.json").write_text(
        json.dumps({"label": "y_convert", "features": feat_cols}, indent=2)
    )
    return grid, grid_dir


def feature_matrix(grid: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    y = grid["y_convert"].astype(int).to_numpy()
    cols = [
        c
        for c in grid.columns
        if c not in LEAK_COLS and pd.api.types.is_numeric_dtype(grid[c])
    ]
    X = grid[cols].astype(np.float64).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return X.to_numpy(np.float64), y, cols


def select_and_train(
    X: np.ndarray,
    y: np.ndarray,
    names: List[str],
    *,
    select_k: int,
    score: str,
    seed: int,
    test_size: float,
) -> dict:
    Xtr, Xte, ytr, yte = train_test_split(
        X, y, test_size=test_size, random_state=seed, stratify=y
    )
    # variance filter then SelectKBest
    vt = VarianceThreshold(1e-8)
    Xtr_v = vt.fit_transform(Xtr)
    Xte_v = vt.transform(Xte)
    kept = [names[i] for i, on in enumerate(vt.get_support()) if on]
    k = min(select_k, len(kept), Xtr_v.shape[1])
    scorer = mutual_info_classif if score == "mi" else f_classif
    skb = SelectKBest(scorer, k=k)
    Xtr_s = skb.fit_transform(Xtr_v, ytr)
    Xte_s = skb.transform(Xte_v)
    sel_idx = skb.get_support(indices=True)
    selected = [
        {"name": kept[i], "score": float(skb.scores_[i])} for i in sel_idx
    ]
    selected.sort(key=lambda d: -d["score"])

    metrics = {}
    # HGB
    t0 = time.time()
    hgb = HistGradientBoostingClassifier(
        max_depth=6, learning_rate=0.08, max_iter=120, random_state=seed
    )
    hgb.fit(Xtr_s, ytr)
    proba = hgb.predict_proba(Xte_s)[:, 1]
    pred = (proba >= 0.5).astype(int)
    metrics["hgb"] = {
        "auc": float(roc_auc_score(yte, proba)),
        "ap": float(average_precision_score(yte, proba)),
        "acc": float(accuracy_score(yte, pred)),
        "sec": float(time.time() - t0),
        "n_selected": int(k),
    }
    # LogReg
    t0 = time.time()
    pipe = Pipeline(
        [
            ("sc", StandardScaler()),
            (
                "lr",
                LogisticRegression(
                    max_iter=400, class_weight="balanced", random_state=seed
                ),
            ),
        ]
    )
    pipe.fit(Xtr_s, ytr)
    proba = pipe.predict_proba(Xte_s)[:, 1]
    pred = (proba >= 0.5).astype(int)
    metrics["logreg"] = {
        "auc": float(roc_auc_score(yte, proba)),
        "ap": float(average_precision_score(yte, proba)),
        "acc": float(accuracy_score(yte, pred)),
        "sec": float(time.time() - t0),
        "n_selected": int(k),
    }
    return {
        "n_train": int(len(ytr)),
        "n_test": int(len(yte)),
        "pos_rate_train": float(ytr.mean()),
        "pos_rate_test": float(yte.mean()),
        "n_feat_in": len(names),
        "n_feat_after_var": len(kept),
        "selected": selected,
        "metrics": metrics,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument(
        "--grid-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_tabular_ui",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_localize_fsds",
    )
    ap.add_argument("--max-users", type=int, default=3000)
    ap.add_argument("--co-window", type=int, default=8)
    ap.add_argument("--rebuild-grid", action="store_true")
    ap.add_argument("--select-k", type=int, default=24)
    ap.add_argument("--score", choices=("f", "mi"), default="f")
    ap.add_argument("--test-size", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--max-rows",
        type=int,
        default=120000,
        help="Subsample grid rows for a fast prototype (stratified)",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    grid, grid_dir = ensure_localize_grid(
        root=args.root,
        grid_dir=args.grid_dir,
        max_users=args.max_users,
        co_window=args.co_window,
        rebuild=args.rebuild_grid,
    )
    print(
        f"grid={grid.shape} pos_rate={float(grid['y_convert'].mean()):.4f}",
        flush=True,
    )

    # stratified subsample for prototype speed
    if args.max_rows and len(grid) > args.max_rows:
        pos = grid[grid["y_convert"] == 1]
        neg = grid[grid["y_convert"] == 0]
        n_pos = min(len(pos), max(500, int(args.max_rows * 0.15)))
        n_neg = min(len(neg), args.max_rows - n_pos)
        grid = pd.concat(
            [
                pos.sample(n_pos, random_state=args.seed),
                neg.sample(n_neg, random_state=args.seed),
            ],
            ignore_index=True,
        )
        print(f"subsample → {grid.shape} pos_rate={float(grid['y_convert'].mean()):.4f}", flush=True)

    X, y, names = feature_matrix(grid)
    print(f"FSDS on {X.shape[1]} localize features, select_k={args.select_k}", flush=True)
    blob = select_and_train(
        X,
        y,
        names,
        select_k=args.select_k,
        score=args.score,
        seed=args.seed,
        test_size=args.test_size,
    )

    summary = {
        "pipeline": ["localization_(user,item)_图谱特征", "FSDS_select_train"],
        "dataset": str(args.root),
        "grid_dir": str(grid_dir),
        "co_window": args.co_window,
        "max_users": args.max_users,
        "select_k": args.select_k,
        "score": args.score,
        "leak_dropped": sorted(LEAK_COLS),
        "feature_names_in": names,
        **{k: blob[k] for k in blob if k != "selected"},
        "top_selected": blob["selected"][:20],
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    (args.out_dir / "selected_features.json").write_text(
        json.dumps(blob["selected"], indent=2)
    )

    md = [
        "# TencentGR prototype: Localization → FSDS",
        "",
        "Scope locked: **no network graph algos**. Two steps only.",
        "",
        "1. **Localization** — node → `(user, item)` 图谱特征网格 + path/covisit 下钻",
        "2. **FSDS** — variance → SelectKBest → HGB / LogReg on `y_convert`",
        "",
        f"- grid source: `{grid_dir}`",
        f"- rows used: **{blob['n_train'] + blob['n_test']}** "
        f"(train {blob['n_train']} / test {blob['n_test']})",
        f"- pos rate test: **{blob['pos_rate_test']:.4f}**",
        f"- features in → after var → selected: "
        f"**{blob['n_feat_in']} → {blob['n_feat_after_var']} → {args.select_k}**",
        f"- leak dropped: `{sorted(LEAK_COLS - {'y_convert', 'user_id', 'item_id', 'e_last_ts'})}` "
        f"(plus id/ts)",
        "",
        "## Metrics",
        "```json",
        json.dumps(blob["metrics"], indent=2),
        "```",
        "",
        "## Top selected (localization features)",
        "```json",
        json.dumps(blob["selected"][:15], indent=2),
        "```",
        "",
        "## Read",
        "- first/last/linear **share_*** columns = 三段启发式 map-back 到 item node",
        "- `u_*` / `i_*` / `e_*` = 图谱聚合列挂在 `(user, item)` edge 上",
        "- FSDS 只在这些列上选 + 训；不再接图算法模块",
        "",
    ]
    (args.out_dir / "LOCALIZE_FSDS_REPORT.md").write_text("\n".join(md))
    print(json.dumps(blob["metrics"], indent=2))
    print("top5:", [s["name"] for s in blob["selected"][:5]])
    print("wrote", args.out_dir)


if __name__ == "__main__":
    main()
