#!/usr/bin/env python3
"""Prototype: two time windows → localization subset → FSDS (no leakage).

Pipeline
--------
1. Propose W1 (early) / W2 (late) with **≥ ``gap_days`` (default 30)** between them.
2. ``feature_engineer(t_start, t_end)`` independently on each window
   (图谱特征 + 三段启发式 localization; no network).
3. Localization subset = top-k items by **W1** ``share_linear`` only.
4. FSDS: variance → SelectKBest → HGB/LogReg on W1 edges in that subset;
   evaluate on **W2** edges (temporal holdout). W2 never enters selection / share fitting.

  PYTHONPATH=. python3 scripts/tencent_gr/run_localize_fsds_time_proto.py \\
    --root data/tencent_subset --max-users 4000 --gap-days 30
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from time_window_feats import (  # noqa: E402
    DAY,
    feature_engineer,
    fsds_feature_columns,
    propose_two_windows,
    scan_time_range,
)

ROOT = Path(__file__).resolve().parents[2]


def _xy(grid: pd.DataFrame, cols: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    X = grid[cols].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(np.float64)
    y = grid["y_convert"].to_numpy(int)
    return X, y


def run_fsds(
    grid_tr: pd.DataFrame,
    grid_te: pd.DataFrame,
    cols: List[str],
    *,
    select_k: int,
    seed: int,
) -> Dict:
    Xtr, ytr = _xy(grid_tr, cols)
    Xte, yte = _xy(grid_te, cols)
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        return {"ok": False, "reason": "need both classes in train and test"}

    k = min(select_k, Xtr.shape[1], max(1, Xtr.shape[0] - 1))
    pipe_prep = Pipeline(
        [
            ("sc", StandardScaler(with_mean=True, with_std=True)),
            ("var", VarianceThreshold(1e-8)),
            ("sel", SelectKBest(f_classif, k=k)),
        ]
    )
    t0 = time.time()
    Xt = pipe_prep.fit_transform(Xtr, ytr)
    Xv = pipe_prep.transform(Xte)
    # recover selected names
    var_mask = pipe_prep.named_steps["var"].get_support()
    cols_var = [c for c, m in zip(cols, var_mask) if m]
    sel_mask = pipe_prep.named_steps["sel"].get_support()
    selected = [c for c, m in zip(cols_var, sel_mask) if m]

    out = {
        "ok": True,
        "n_train": int(len(ytr)),
        "n_test": int(len(yte)),
        "pos_train": float(ytr.mean()),
        "pos_test": float(yte.mean()),
        "n_in": len(cols),
        "n_selected": len(selected),
        "selected": selected,
        "sec_select": float(time.time() - t0),
        "models": {},
    }

    # HGB
    t0 = time.time()
    hgb = HistGradientBoostingClassifier(
        max_depth=6,
        learning_rate=0.08,
        max_iter=150,
        random_state=seed,
        class_weight="balanced",
    )
    hgb.fit(Xt, ytr)
    ph = hgb.predict_proba(Xv)[:, 1]
    out["models"]["hgb"] = {
        "auc": float(roc_auc_score(yte, ph)),
        "ap": float(average_precision_score(yte, ph)),
        "sec": float(time.time() - t0),
    }

    # LogReg (already standardized in pipe_prep)
    t0 = time.time()
    lr = LogisticRegression(max_iter=500, C=0.5, class_weight="balanced", random_state=seed)
    lr.fit(Xt, ytr)
    pl = lr.predict_proba(Xv)[:, 1]
    out["models"]["logreg"] = {
        "auc": float(roc_auc_score(yte, pl)),
        "ap": float(average_precision_score(yte, pl)),
        "sec": float(time.time() - t0),
    }
    return out


def shares_from_convert_df(convert_df: pd.DataFrame) -> pd.DataFrame:
    """Recompute first/last/linear shares from localization rows (train users only)."""
    if convert_df is None or len(convert_df) == 0:
        return pd.DataFrame(
            columns=[
                "item_id",
                "share_first",
                "share_last",
                "share_linear",
                "credit_first",
                "credit_last",
                "credit_linear",
                "n_as_convert_terminal",
            ]
        )
    first = convert_df[convert_df["touch"] == "first"].groupby("path_item").size()
    last = convert_df[convert_df["touch"] == "last"].groupby("path_item").size()
    lin = convert_df.groupby("path_item")["linear_credit"].sum()
    term = convert_df[convert_df["is_convert_item"] == 1].groupby("path_item").size()
    items = sorted(set(first.index) | set(last.index) | set(lin.index) | set(term.index))
    sf, sl, sn = float(first.sum() or 1.0), float(last.sum() or 1.0), float(lin.sum() or 1.0)
    rows = []
    for iid in items:
        cf = float(first.get(iid, 0.0))
        cl = float(last.get(iid, 0.0))
        cn = float(lin.get(iid, 0.0))
        rows.append(
            {
                "item_id": int(iid),
                "credit_first": cf,
                "credit_last": cl,
                "credit_linear": cn,
                "share_first": cf / sf,
                "share_last": cl / sl,
                "share_linear": cn / sn,
                "n_as_convert_terminal": int(term.get(iid, 0)),
            }
        )
    return pd.DataFrame(rows)


def attach_train_shares(grid: pd.DataFrame, share_df: pd.DataFrame) -> pd.DataFrame:
    """Overwrite localization share columns with train-only aggregates."""
    g = grid.copy()
    drop_cols = [
        c
        for c in g.columns
        if c.startswith("i_share_")
        or c.startswith("i_credit_")
        or c == "i_n_as_convert_terminal"
        or c == "i_item_credit_rank"
    ]
    g = g.drop(columns=drop_cols, errors="ignore")
    if share_df is None or len(share_df) == 0:
        for c in (
            "i_share_first",
            "i_share_last",
            "i_share_linear",
            "i_credit_first",
            "i_credit_last",
            "i_credit_linear",
        ):
            g[c] = 0.0
        return g
    s = share_df.rename(
        columns={
            "share_first": "i_share_first",
            "share_last": "i_share_last",
            "share_linear": "i_share_linear",
            "credit_first": "i_credit_first",
            "credit_last": "i_credit_last",
            "credit_linear": "i_credit_linear",
        }
    )
    keep = [
        "item_id",
        "i_share_first",
        "i_share_last",
        "i_share_linear",
        "i_credit_first",
        "i_credit_last",
        "i_credit_linear",
    ]
    g = g.merge(s[keep], on="item_id", how="left")
    for c in keep[1:]:
        g[c] = g[c].fillna(0.0)
    g["i_item_credit_rank"] = g["i_share_linear"].rank(method="average", ascending=False)
    return g


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=4000)
    ap.add_argument("--window-days", type=int, default=45)
    ap.add_argument("--gap-days", type=int, default=30, help="min gap between W1 end and W2 start")
    ap.add_argument("--localize-k", type=int, default=200)
    ap.add_argument("--select-k", type=int, default=32)
    ap.add_argument("--co-window", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--test-size", type=float, default=0.25)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_localize_fsds_time",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("scan timeline ...", flush=True)
    t_min, t_max, n_scan = scan_time_range(args.root, max_users=args.max_users)
    w1, w2 = propose_two_windows(
        t_min, t_max, window_days=args.window_days, gap_days=args.gap_days
    )
    gap_d = (w2.t_start - w1.t_end) / DAY
    print(
        f"span={ (t_max-t_min)/DAY:.1f}d users_scan_cap={args.max_users} "
        f"W1=[{w1.t_start},{w1.t_end}) {w1.n_days:.1f}d  "
        f"gap={gap_d:.1f}d  "
        f"W2=[{w2.t_start},{w2.t_end}) {w2.n_days:.1f}d",
        flush=True,
    )
    if gap_d < args.gap_days - 1e-6:
        raise SystemExit(f"gap {gap_d:.2f}d < required {args.gap_days}d")

    print("FE W1 ...", flush=True)
    pack1 = feature_engineer(
        args.root,
        w1.t_start,
        w1.t_end,
        max_users=args.max_users,
        co_window=args.co_window,
        window_name=w1.name,
        terminal_action=1,  # TencentGR-1M: click is the labeled success
    )
    print("FE W2 ...", flush=True)
    pack2 = feature_engineer(
        args.root,
        w2.t_start,
        w2.t_end,
        max_users=args.max_users,
        co_window=args.co_window,
        window_name=w2.name,
        terminal_action=1,
    )
    print(
        f"W1 edges={pack1['meta']['n_edges']} pos={pack1['meta']['pos_rate']:.4f} | "
        f"W2 edges={pack2['meta']['n_edges']} pos={pack2['meta']['pos_rate']:.4f}",
        flush=True,
    )

    # --- no-leak localization: user split first, shares from train users only ---
    g1 = pack1["grid"]
    g2 = pack2["grid"]
    if len(g1) < 80 or g1["y_convert"].nunique() < 2:
        raise SystemExit("W1 grid too small / single-class")

    gss = GroupShuffleSplit(n_splits=1, test_size=args.test_size, random_state=args.seed)
    tr_idx, va_idx = next(gss.split(g1, g1["y_convert"], g1["user_id"]))
    train_users = set(g1.iloc[tr_idx]["user_id"].tolist())
    conv_tr = pack1["convert_df"][pack1["convert_df"]["user_id"].isin(train_users)]
    share_tr = shares_from_convert_df(conv_tr)
    loc_items = (
        share_tr.sort_values("share_linear", ascending=False)
        .head(args.localize_k)["item_id"]
        .astype(int)
        .tolist()
    )
    loc_set = set(loc_items)
    print(f"localize subset from W1-train users only: k={len(loc_items)}", flush=True)

    g1 = attach_train_shares(g1, share_tr)
    g2 = attach_train_shares(g2, share_tr)  # freeze W1-train shares onto W2
    g1_loc = g1[g1["item_id"].isin(loc_set)].copy()
    g2_loc = g2[g2["item_id"].isin(loc_set)].copy()
    g1_tr = g1_loc[g1_loc["user_id"].isin(train_users)]
    g1_va = g1_loc[~g1_loc["user_id"].isin(train_users)]
    print(
        f"localized edges W1-train={len(g1_tr)} W1-holdout={len(g1_va)} W2={len(g2_loc)}",
        flush=True,
    )
    if len(g1_tr) < 40 or g1_tr["y_convert"].sum() < 5 or g1_tr["y_convert"].nunique() < 2:
        raise SystemExit(
            f"W1 train localized grid too small / few positives "
            f"(n={len(g1_tr)} pos={int(g1_tr['y_convert'].sum())})"
        )

    cols = fsds_feature_columns(g1_tr)
    cols = [c for c in cols if c in g2_loc.columns]
    print(f"FSDS candidate cols={len(cols)}", flush=True)

    res_w1 = run_fsds(g1_tr, g1_va, cols, select_k=args.select_k, seed=args.seed)
    res_w2 = run_fsds(g1_tr, g2_loc, cols, select_k=args.select_k, seed=args.seed)

    # save tables
    pack1["item_df"].to_parquet(args.out_dir / "W1_item_features.parquet", index=False)
    pack2["item_df"].to_parquet(args.out_dir / "W2_item_features.parquet", index=False)
    g1_loc = pd.concat([g1_tr, g1_va], ignore_index=True)
    g1_loc.to_parquet(args.out_dir / "W1_localized_grid.parquet", index=False)
    g2_loc.head(50000).to_parquet(args.out_dir / "W2_localized_grid_sample.parquet", index=False)
    conv_tr.to_csv(args.out_dir / "W1_train_localize_convert_path.csv", index=False)
    share_tr.to_csv(args.out_dir / "W1_train_item_shares.csv", index=False)
    pd.DataFrame({"item_id": loc_items, "rank": np.arange(1, len(loc_items) + 1)}).to_csv(
        args.out_dir / "localize_subset_from_W1_train.csv", index=False
    )

    blob = {
        "protocol": [
            "FE(time_range) on W1 and W2 independently",
            "gap_days >= 30 between W1.end and W2.start",
            "user-grouped split on W1; localization shares from W1-train users only",
            "localization subset = top-k share_linear from W1-train",
            "FSDS select+train on W1-train localized edges",
            "evaluate on W1-holdout users and W2 temporal holdout",
        ],
        "timeline": {
            "t_min": t_min,
            "t_max": t_max,
            "span_days": (t_max - t_min) / DAY,
            "W1": w1.to_dict(),
            "W2": w2.to_dict(),
            "gap_days": gap_d,
        },
        "W1_meta": pack1["meta"],
        "W2_meta": pack2["meta"],
        "localize_k": args.localize_k,
        "n_localize_items": len(loc_items),
        "n_W1_train_loc_edges": int(len(g1_tr)),
        "n_W1_holdout_loc_edges": int(len(g1_va)),
        "n_W2_loc_edges": int(len(g2_loc)),
        "feature_cols": cols,
        "fsds_W1_holdout": res_w1,
        "fsds_W2_temporal": res_w2,
        "leakage_guards": [
            "events outside [t_start,t_end) never enter FE",
            "y / n_clk / n_cnv / ctr / cvr dropped from X",
            "localization shares fit on W1-train users only (frozen onto W2)",
            "SelectKBest fit on W1-train only; W2 never used for selection",
        ],
    }
    (args.out_dir / "summary.json").write_text(json.dumps(blob, indent=2, default=str))

    # plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    labs = ["W1→W1 user-holdout", "W1→W2 temporal"]
    for i, res in enumerate((res_w1, res_w2)):
        if not res.get("ok"):
            continue
        ax.bar(i - 0.15, res["models"]["hgb"]["auc"], width=0.3, color="#4C78A8", label="HGB" if i == 0 else None)
        ax.bar(i + 0.15, res["models"]["logreg"]["auc"], width=0.3, color="#F58518", label="LogReg" if i == 0 else None)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel("AUC")
    ax.set_title("FSDS after localization")
    ax.legend(fontsize=8)
    ax = axes[1]
    if len(share_tr):
        top = share_tr.sort_values("share_linear", ascending=False).head(15)
        ax.barh(range(len(top)), top["share_linear"][::-1], color="#54A24B")
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([str(i) for i in top["item_id"][::-1]], fontsize=7)
        ax.set_xlabel("W1-train share_linear")
        ax.set_title("Localization subset (top shares)")
    fig.suptitle(
        f"TencentGR localize→FSDS  gap={gap_d:.0f}d  k={len(loc_items)}",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(args.out_dir / "localize_fsds_time.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    def _fmt(res: Dict) -> str:
        if not res.get("ok"):
            return f"FAIL ({res.get('reason')})"
        h, l = res["models"]["hgb"], res["models"]["logreg"]
        return (
            f"HGB AUC={h['auc']:.3f} AP={h['ap']:.3f} | "
            f"LR AUC={l['auc']:.3f} AP={l['ap']:.3f} | "
            f"n={res['n_train']}/{res['n_test']} sel={res['n_selected']}"
        )

    md = [
        "# TencentGR prototype: time-window localization → FSDS",
        "",
        "两步：**Localization → FSDS**。FE 带时间范围；两窗间隔 ≥1 个月；无 network。",
        "",
        "> TencentGR-1M labels = exposure(0)/click(1)；terminal success = **click**。",
        "",
        "## Protocol (no leakage)",
        "1. `feature_engineer(root, t_start, t_end)` — 只用窗内事件",
        f"2. W1 early / W2 late，gap = **{gap_d:.1f}** days (≥ {args.gap_days})",
        "3. Localization subset = W1 `share_linear` top-k（三段启发式 aggregate through）",
        "4. FSDS select+train 只在 W1；W2 仅作 temporal holdout",
        "",
        "## Windows",
        f"- timeline span: **{(t_max-t_min)/DAY:.1f}** days",
        f"- W1: `{w1.name}` [{w1.t_start}, {w1.t_end}) ≈ {w1.n_days:.1f}d — "
        f"edges={pack1['meta']['n_edges']} pos={pack1['meta']['pos_rate']:.4f}",
        f"- W2: `{w2.name}` [{w2.t_start}, {w2.t_end}) ≈ {w2.n_days:.1f}d — "
        f"edges={pack2['meta']['n_edges']} pos={pack2['meta']['pos_rate']:.4f}",
        f"- localize-k: **{len(loc_items)}** (W1-train shares) → "
        f"W1-train={len(g1_tr)}, W1-holdout={len(g1_va)}, W2={len(g2_loc)}",
        "",
        "## FSDS",
        f"- W1→W1 user-holdout: {_fmt(res_w1)}",
        f"- W1→W2 temporal: {_fmt(res_w2)}",
        "",
        "## Selected features (W1 fit)",
        ", ".join((res_w2.get("selected") or res_w1.get("selected") or [])[:24]) or "(none)",
        "",
        "## Guards",
        "- no events outside window in FE",
        "- drop `y_convert` / `*_n_cnv` / `*cvr*` from X",
        "- subset + SelectKBest never see W2 labels",
        "",
    ]
    (args.out_dir / "LOCALIZE_FSDS_TIME_REPORT.md").write_text("\n".join(md))
    print("wrote", args.out_dir)
    print(md[md.index("## FSDS") + 1])
    print(md[md.index("## FSDS") + 2])


if __name__ == "__main__":
    main()
