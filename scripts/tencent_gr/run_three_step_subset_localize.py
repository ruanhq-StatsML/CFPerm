#!/usr/bin/env python3
"""Three-step subset localization: merchant → user → order.

Methodology (not causal; modality-agnostic multi-granularity drill-down):
  find the drifting cohort faster via hierarchical subset localization.

  1. Standardization (fit on W1)
  2. Merchant-level MMD² → top merchants
  3. User-level MMD² (within those merchants) → top users
  4. Order-level score (within those users) → top orders
  5. FSDS ranking on the localized order edges

No causality claim — distribution shift ranking only.
Same procedure works at other entity granularities / modalities.

Merchant id = item_feat column (default ``122``, encrypted advertiser/shop proxy).
Order id   = ``{user_id}_{item_id}_{ts}`` for terminal-success edges (click on TencentGR).

  PYTHONPATH=. python3 scripts/tencent_gr/run_three_step_subset_localize.py \\
    --root data/tencent_subset --max-users 20000 --gap-days 30
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

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
from run_standardize_mmd_fsds import (  # noqa: E402
    _matrix,
    fit_standardizer,
    rbf_mmd2,
    run_fsds,
    standardize,
    w1w2_candidate_columns,
)
from direction_report import build_direction_dict  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def load_item_merchant_map(root: Path, *, merchant_col: str = "122") -> pd.DataFrame:
    """item_id → merchant_id from item_feat (encrypted shop/advertiser proxy)."""
    paths = sorted((root / "item_feat").glob("*.parquet"))
    if not paths:
        raise FileNotFoundError(f"no item_feat under {root}")
    frames = []
    for p in paths:
        cols = pq.ParquetFile(p).schema.names
        need = ["item_id", merchant_col]
        if merchant_col not in cols:
            raise KeyError(f"{merchant_col} not in item_feat columns {cols}")
        frames.append(pq.read_table(p, columns=need).to_pandas())
    m = pd.concat(frames, ignore_index=True)
    m = m.dropna(subset=[merchant_col])
    m["item_id"] = m["item_id"].astype(int)
    m["merchant_id"] = pd.to_numeric(m[merchant_col], errors="coerce").astype("Int64")
    m = m.dropna(subset=["merchant_id"])
    m["merchant_id"] = m["merchant_id"].astype(int)
    return m[["item_id", "merchant_id"]].drop_duplicates("item_id")


def attach_merchant(grid: pd.DataFrame, imap: pd.DataFrame) -> pd.DataFrame:
    g = grid.merge(imap, on="item_id", how="left")
    # real map vs singleton proxy (unmapped item acts as its own shop)
    g["merchant_mapped"] = g["merchant_id"].notna().astype(int)
    miss = g["merchant_id"].isna()
    g.loc[miss, "merchant_id"] = g.loc[miss, "item_id"].astype(int)
    g["merchant_id"] = g["merchant_id"].astype(int)
    return g


def _entity_mmd(
    g1: pd.DataFrame,
    g2: pd.DataFrame,
    cols: Sequence[str],
    sc,
    *,
    entity_col: str,
    seed: int,
    max_cand: int,
    min_edges: int,
    mmd_max_n: int,
) -> pd.DataFrame:
    """Rank entities by RBF-MMD² between W1/W2 rows (standardized)."""
    if len(g1) == 0 or len(g2) == 0 or entity_col not in g1.columns:
        return pd.DataFrame(columns=[entity_col, "n_W1", "n_W2", "mmd2", "rank", "selected"])
    c1 = g1[entity_col].value_counts()
    c2 = g2[entity_col].value_counts()
    common = sorted(set(c1[c1 >= min_edges].index) & set(c2[c2 >= min_edges].index))
    vol = {e: int(c1.get(e, 0) + c2.get(e, 0)) for e in common}
    cand = sorted(common, key=lambda e: -vol[e])[:max_cand]
    rng = np.random.default_rng(seed)
    rows = []
    for e in cand:
        a = standardize(sc, _matrix(g1[g1[entity_col] == e], cols))
        b = standardize(sc, _matrix(g2[g2[entity_col] == e], cols))
        rows.append(
            {
                entity_col: int(e) if isinstance(e, (int, np.integer)) else e,
                "n_W1": int(len(a)),
                "n_W2": int(len(b)),
                "mmd2": rbf_mmd2(a, b, max_n=mmd_max_n, rng=rng),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=[entity_col, "n_W1", "n_W2", "mmd2", "rank", "selected"])
    out = out.sort_values("mmd2", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out["selected"] = 0
    return out


def _entity_mean_shift(
    g1: pd.DataFrame,
    g2: pd.DataFrame,
    cols: Sequence[str],
    sc,
    *,
    entity_col: str,
    max_cand: int,
    min_edges: int = 1,
) -> pd.DataFrame:
    """Fallback when MMD lacks paired mass: rank by ‖μ2−μ1‖₂."""
    if len(g1) == 0 or len(g2) == 0:
        return pd.DataFrame(columns=[entity_col, "n_W1", "n_W2", "mmd2", "rank", "selected"])
    c1 = g1[entity_col].value_counts()
    c2 = g2[entity_col].value_counts()
    # union with at least min_edges in one window and ≥1 in both if possible
    both = set(c1[c1 >= min_edges].index) & set(c2[c2 >= 1].index)
    if len(both) < 5:
        both = set(c1[c1 >= min_edges].index) | set(c2[c2 >= min_edges].index)
    vol = {e: int(c1.get(e, 0) + c2.get(e, 0)) for e in both}
    cand = sorted(both, key=lambda e: -vol[e])[:max_cand]
    rows = []
    for e in cand:
        a = standardize(sc, _matrix(g1[g1[entity_col] == e], cols)) if e in c1.index else None
        b = standardize(sc, _matrix(g2[g2[entity_col] == e], cols)) if e in c2.index else None
        if a is None or len(a) == 0:
            score = float(np.linalg.norm(b.mean(axis=0))) if b is not None and len(b) else 0.0
            n1, n2 = 0, int(len(b) if b is not None else 0)
        elif b is None or len(b) == 0:
            score = float(np.linalg.norm(a.mean(axis=0)))
            n1, n2 = int(len(a)), 0
        else:
            score = float(np.linalg.norm(b.mean(axis=0) - a.mean(axis=0)))
            n1, n2 = int(len(a)), int(len(b))
        rows.append(
            {
                entity_col: int(e) if isinstance(e, (int, np.integer)) else e,
                "n_W1": n1,
                "n_W2": n2,
                "mmd2": score,  # reuse column name for ranking/plot
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=[entity_col, "n_W1", "n_W2", "mmd2", "rank", "selected"])
    out = out.sort_values("mmd2", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out["selected"] = 0
    return out


def score_orders(
    g1: pd.DataFrame,
    g2: pd.DataFrame,
    cols: Sequence[str],
    sc,
    *,
    user_set: Set[int],
    seed: int,
    top_k: int,
) -> pd.DataFrame:
    """Order ≈ terminal-success edge (or any edge) for selected users.

    Score = ‖x_W2 − μ_W1_user‖₂ on standardized features when W2 row exists;
    else ‖x_W1 − μ_W1_user‖₂ as within-window extremity (still a localize signal).
    """
    g1u = g1[g1["user_id"].isin(user_set)].copy()
    g2u = g2[g2["user_id"].isin(user_set)].copy()
    # prefer positive terminal rows as "orders"
    def _orders(g: pd.DataFrame) -> pd.DataFrame:
        if "y_convert" in g.columns and int(g["y_convert"].sum()) > 0:
            o = g[g["y_convert"] == 1].copy()
        else:
            o = g.copy()
        ts = o["e_last_ts"] if "e_last_ts" in o.columns else 0
        o = o.copy()
        o["order_id"] = (
            o["user_id"].astype(str)
            + "_"
            + o["item_id"].astype(str)
            + "_"
            + pd.Series(ts, index=o.index).astype(str)
        )
        return o

    o1, o2 = _orders(g1u), _orders(g2u)
    # user means on W1 (standardized)
    X1 = standardize(sc, _matrix(g1u, cols)) if len(g1u) else np.zeros((0, len(cols)))
    uid1 = g1u["user_id"].to_numpy() if len(g1u) else np.array([], dtype=int)
    mu: Dict[int, np.ndarray] = {}
    for u in user_set:
        m = uid1 == u
        if m.any():
            mu[u] = X1[m].mean(axis=0)

    rows = []
    # score W2 orders primarily (drift visible in late window)
    pool = o2 if len(o2) else o1
    X = standardize(sc, _matrix(pool, cols))
    for i, (_, r) in enumerate(pool.iterrows()):
        u = int(r["user_id"])
        x = X[i]
        base = mu.get(u)
        if base is None:
            score = float(np.linalg.norm(x))
        else:
            score = float(np.linalg.norm(x - base))
        rows.append(
            {
                "order_id": r["order_id"],
                "user_id": u,
                "item_id": int(r["item_id"]),
                "merchant_id": int(r.get("merchant_id", -1)),
                "y_convert": int(r.get("y_convert", 0)),
                "shift_l2": score,
                "window": "W2" if len(o2) else "W1",
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values("shift_l2", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    out["selected"] = (out.index < top_k).astype(int)
    return out


def plot_three_step(
    mer: pd.DataFrame,
    usr: pd.DataFrame,
    ord_: pd.DataFrame,
    feat_rank: pd.DataFrame,
    *,
    out_path: Path,
    title: str,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    top = mer.head(12)
    ax.barh(range(len(top)), top["mmd2"][::-1], color="#4C78A8")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([str(x) for x in top["merchant_id"][::-1]], fontsize=7)
    ax.set_xlabel("MMD²")
    ax.set_title("L1 merchant subset")

    ax = axes[0, 1]
    top = usr.head(12)
    ax.barh(range(len(top)), top["mmd2"][::-1], color="#F58518")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([str(x) for x in top["user_id"][::-1]], fontsize=7)
    ax.set_xlabel("MMD²")
    ax.set_title("L2 user subset (in merchants)")

    ax = axes[1, 0]
    top = ord_.head(12)
    ax.barh(range(len(top)), top["shift_l2"][::-1], color="#E45756")
    ax.set_yticks(range(len(top)))
    labels = [f"{int(u)}|{int(i)}" for u, i in zip(top["user_id"], top["item_id"])]
    ax.set_yticklabels(labels[::-1], fontsize=7)
    ax.set_xlabel("‖x − μ_user‖₂")
    ax.set_title("L3 order subset (in users)")

    ax = axes[1, 1]
    top = feat_rank.head(12)
    ax.barh(range(len(top)), top["f_score"][::-1], color="#54A24B")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["feature"][::-1], fontsize=8)
    ax.set_xlabel("FSDS F-score")
    ax.set_title("FSDS ranking on localized orders")

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="merchant → user → order subset localization")
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=20000)
    ap.add_argument("--window-days", type=int, default=45)
    ap.add_argument("--gap-days", type=int, default=30)
    ap.add_argument("--merchant-col", type=str, default="122")
    ap.add_argument("--k-merchant", type=int, default=40)
    ap.add_argument("--k-user", type=int, default=80)
    ap.add_argument("--k-order", type=int, default=200)
    ap.add_argument("--max-cand-merchant", type=int, default=300)
    ap.add_argument("--max-cand-user", type=int, default=400)
    ap.add_argument("--min-edges", type=int, default=3)
    ap.add_argument("--mmd-max-n", type=int, default=128)
    ap.add_argument("--select-k", type=int, default=15)
    ap.add_argument("--co-window", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_three_step_localize",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("load item→merchant map ...", flush=True)
    imap = load_item_merchant_map(args.root, merchant_col=args.merchant_col)
    print(f"  items={len(imap)} merchants={imap['merchant_id'].nunique()} col={args.merchant_col}", flush=True)

    print("scan timeline ...", flush=True)
    t_min, t_max, _ = scan_time_range(args.root, max_users=args.max_users)
    w1, w2 = propose_two_windows(
        t_min, t_max, window_days=args.window_days, gap_days=args.gap_days
    )
    gap_d = (w2.t_start - w1.t_end) / DAY
    print(f"gap={gap_d:.0f}d  W1={w1.n_days:.0f}d  W2={w2.n_days:.0f}d", flush=True)

    print("FE W1 / W2 ...", flush=True)
    p1 = feature_engineer(
        args.root, w1.t_start, w1.t_end,
        max_users=args.max_users, co_window=args.co_window,
        window_name=w1.name, terminal_action=1,
    )
    p2 = feature_engineer(
        args.root, w2.t_start, w2.t_end,
        max_users=args.max_users, co_window=args.co_window,
        window_name=w2.name, terminal_action=1,
    )
    g1 = attach_merchant(p1["grid"], imap)
    g2 = attach_merchant(p2["grid"], imap)
    print(f"W1={len(g1)} W2={len(g2)}", flush=True)

    raw_cols = [c for c in fsds_feature_columns(g1) if c in g2.columns]
    cols = w1w2_candidate_columns(raw_cols)
    cols = [c for c in cols if c not in ("merchant_id", "order_id", "merchant_mapped")]
    print(
        f"features n={len(cols)}  merchant_mapped_rate="
        f"{float(g1['merchant_mapped'].mean()):.3f}",
        flush=True,
    )

    print("1/4 Standardization ...", flush=True)
    sc = fit_standardizer(_matrix(g1, cols))

    print("2/4 L1 merchant MMD ...", flush=True)
    # prefer real mapped merchants; fall back to all if too sparse
    g1_mpref = g1[g1["merchant_mapped"] == 1]
    g2_mpref = g2[g2["merchant_mapped"] == 1]
    use_g1, use_g2 = g1_mpref, g2_mpref
    if len(g1_mpref) < 200 or len(g2_mpref) < 200:
        print("  mapped merchants sparse in windows → use all (item-as-shop proxy ok)", flush=True)
        use_g1, use_g2 = g1, g2
    mer = _entity_mmd(
        use_g1, use_g2, cols, sc,
        entity_col="merchant_id", seed=args.seed,
        max_cand=args.max_cand_merchant, min_edges=args.min_edges,
        mmd_max_n=args.mmd_max_n,
    )
    if len(mer) < 5:
        print("  merchant MMD sparse → mean-shift fallback", flush=True)
        mer = _entity_mean_shift(
            use_g1, use_g2, cols, sc, entity_col="merchant_id",
            max_cand=args.max_cand_merchant, min_edges=1,
        )
    k_m = min(args.k_merchant, len(mer))
    mer["selected"] = 0
    if k_m:
        mer.loc[: k_m - 1, "selected"] = 1
    mer_ids = set(mer.loc[mer["selected"] == 1, "merchant_id"].astype(int)) if len(mer) else set()
    g1m = g1[g1["merchant_id"].isin(mer_ids)] if mer_ids else g1
    g2m = g2[g2["merchant_id"].isin(mer_ids)] if mer_ids else g2
    print(f"  merchants k={len(mer_ids)} edges W1={len(g1m)} W2={len(g2m)}", flush=True)

    print("3/4 L2 user MMD (in merchants) ...", flush=True)
    usr = _entity_mmd(
        g1m, g2m, cols, sc,
        entity_col="user_id", seed=args.seed + 1,
        max_cand=args.max_cand_user, min_edges=2,
        mmd_max_n=args.mmd_max_n,
    )
    if len(usr) < 5:
        print("  user MMD sparse → mean-shift fallback", flush=True)
        usr = _entity_mean_shift(
            g1m, g2m, cols, sc, entity_col="user_id",
            max_cand=args.max_cand_user, min_edges=1,
        )
    k_u = min(args.k_user, len(usr))
    usr["selected"] = 0
    if k_u:
        usr.loc[: k_u - 1, "selected"] = 1
    user_ids = set(usr.loc[usr["selected"] == 1, "user_id"].astype(int)) if len(usr) else set()
    if not user_ids:
        # last resort: top users by volume in merchant subset
        vol = (g1m["user_id"].value_counts().add(g2m["user_id"].value_counts(), fill_value=0))
        user_ids = set(vol.head(args.k_user).index.astype(int).tolist())
        print(f"  user fallback by volume k={len(user_ids)}", flush=True)
    print(f"  users k={len(user_ids)}", flush=True)

    print("4/4 L3 order localize + FSDS ...", flush=True)
    orders = score_orders(
        g1m, g2m, cols, sc, user_set=user_ids, seed=args.seed + 2, top_k=args.k_order
    )
    sel_orders = orders[orders["selected"] == 1] if len(orders) else orders
    # FSDS on user×merchant localized mass (orders drive viz; need enough labels)
    g1_fs = g1m[g1m["user_id"].isin(user_ids)].copy()
    g2_fs = g2m[g2m["user_id"].isin(user_ids)].copy()
    if len(g1_fs) < 80:
        g1_fs, g2_fs = g1m.copy(), g2m.copy()
    if len(g1_fs) >= 40:
        g1_tr = g1_fs.sample(frac=0.75, random_state=args.seed)
        g1_va = g1_fs.drop(g1_tr.index)
    else:
        g1_tr, g1_va = g1_fs, g1_fs

    # drop id cols from modeling frames
    for _g in (g1_tr, g1_va, g2_fs):
        for c in ("merchant_id", "merchant_mapped"):
            if c in _g.columns:
                # keep column in df for reporting but run_fsds only uses cols
                pass

    res_w1 = run_fsds(g1_tr, g1_va, cols, select_k=args.select_k, seed=args.seed)
    res_w2 = run_fsds(g1_tr, g2_fs, cols, select_k=args.select_k, seed=args.seed)
    ranking = res_w1.get("ranking")
    if not isinstance(ranking, pd.DataFrame) or ranking.empty:
        ranking = res_w2.get("ranking")
    if not isinstance(ranking, pd.DataFrame) or ranking.empty:
        ranking = pd.DataFrame(
            {"feature": cols, "f_score": 0.0, "rank": range(1, len(cols) + 1), "selected": 0}
        )

    mer.to_csv(args.out_dir / "L1_merchant_mmd.csv", index=False)
    usr.to_csv(args.out_dir / "L2_user_mmd.csv", index=False)
    orders.to_csv(args.out_dir / "L3_order_scores.csv", index=False)
    sel_orders.to_csv(args.out_dir / "localized_orders.csv", index=False)
    ranking.to_csv(args.out_dir / "fsds_feature_ranking.csv", index=False)

    def _strip(res: Dict) -> Dict:
        out = {k: v for k, v in res.items() if k != "ranking"}
        if isinstance(res.get("ranking"), pd.DataFrame):
            out["top_features"] = res["ranking"].head(args.select_k)["feature"].tolist()
        return out

    tip_feats = ranking.head(args.select_k)["feature"].astype(str).tolist() if len(ranking) else []
    # feature means on localized support for tip signs
    X1 = standardize(sc, _matrix(g1_fs, cols)) if len(g1_fs) else np.zeros((0, len(cols)))
    X2 = standardize(sc, _matrix(g2_fs, cols)) if len(g2_fs) else np.zeros((0, len(cols)))
    feat_rows = []
    for j, c in enumerate(cols):
        m1 = float(X1[:, j].mean()) if len(X1) else 0.0
        m2 = float(X2[:, j].mean()) if len(X2) else 0.0
        feat_rows.append({"feature": c, "mean_W1": m1, "mean_W2": m2, "cmean_abs": abs(m2 - m1)})
    feat_diag = pd.DataFrame(feat_rows)
    direction = build_direction_dict(
        g1_fs,
        g2_fs,
        tip_feats,
        feat_diag=feat_diag,
        y_col="y_convert",
        extra={"support": "merchant_user", "k": {"merchant": len(mer_ids), "user": len(user_ids)}},
    )

    blob = {
        "procedure": [
            "Standardization",
            "L1 merchant MMD",
            "L2 user MMD (in merchants)",
            "L3 order shift (in users)",
            "FSDS on localized merchant×user edges",
            "direction JSON: sign(Δȳ) + tip sign(δ_j) on K*",
        ],
        "merchant_col": args.merchant_col,
        "merchant_mapped_rate_W1": float(g1["merchant_mapped"].mean()),
        "gap_days": gap_d,
        "k": {"merchant": len(mer_ids), "user": len(user_ids), "order": int(len(sel_orders))},
        "n_edges": {"W1_loc": int(len(g1_fs)), "W2_loc": int(len(g2_fs))},
        "fsds_W1": _strip(res_w1),
        "fsds_W2": _strip(res_w2),
        "direction": direction,
        "W1_meta": p1["meta"],
        "W2_meta": p2["meta"],
    }
    (args.out_dir / "summary.json").write_text(json.dumps(blob, indent=2, default=str))
    print("direction:", direction.get("report"), flush=True)

    plot_three_step(
        mer[mer["selected"] == 1] if len(mer) else mer,
        usr[usr["selected"] == 1] if len(usr) else usr,
        sel_orders if len(sel_orders) else orders,
        ranking,
        out_path=args.out_dir / "three_step_subset_localize.png",
        title=f"3-step subset localize  merchant→user→order  gap={gap_d:.0f}d",
    )

    md = [
        "# Three-step subset localization: merchant → user → order",
        "",
        "```",
        "Standardize → L1 merchant MMD → L2 user MMD → L3 order shift → FSDS",
        "```",
        "",
        f"- merchant_col = `{args.merchant_col}` (item_feat; unmapped item → singleton shop proxy)",
        f"- mapped_rate W1 = **{float(g1['merchant_mapped'].mean()):.3f}**",
        f"- gap = **{gap_d:.0f}**d | k_m={len(mer_ids)} k_u={len(user_ids)} k_o={len(sel_orders)}",
        f"- localized edges: W1={len(g1_fs)} W2={len(g2_fs)}",
        "",
        "## L1 merchants (head)",
        "| rank | merchant_id | MMD² | n_W1 | n_W2 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for _, r in mer.head(10).iterrows():
        md.append(
            f"| {int(r['rank'])} | {int(r['merchant_id'])} | {r['mmd2']:.4f} | "
            f"{int(r['n_W1'])} | {int(r['n_W2'])} |"
        )
    md += [
        "",
        "## L2 users (head)",
        "| rank | user_id | MMD² | n_W1 | n_W2 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for _, r in usr.head(10).iterrows():
        md.append(
            f"| {int(r['rank'])} | {int(r['user_id'])} | {r['mmd2']:.4f} | "
            f"{int(r['n_W1'])} | {int(r['n_W2'])} |"
        )
    md += [
        "",
        "## L3 orders (head)",
        "| rank | user_id | item_id | merchant_id | shift_l2 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for _, r in (sel_orders if len(sel_orders) else orders).head(10).iterrows():
        md.append(
            f"| {int(r['rank'])} | {int(r['user_id'])} | {int(r['item_id'])} | "
            f"{int(r['merchant_id'])} | {r['shift_l2']:.3f} |"
        )
    md += [
        "",
        "## FSDS ranking",
        "| rank | feature | F-score |",
        "|---:|---|---:|",
    ]
    for _, r in ranking.head(args.select_k).iterrows():
        md.append(f"| {int(r['rank'])} | `{r['feature']}` | {float(r['f_score']):.3f} |")
    md.append("")
    (args.out_dir / "THREE_STEP_LOCALIZE_REPORT.md").write_text("\n".join(md))
    print("wrote", args.out_dir)
    print("top features:", ", ".join(ranking.head(8)["feature"].astype(str).tolist()))


if __name__ == "__main__":
    main()
