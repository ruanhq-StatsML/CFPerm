#!/usr/bin/env python3
"""TencentGR 图谱特征 → 网格面板（不用图算法）.

主维只有 user / item。所谓「图谱特征」= 关系型聚合列
（度/共现强度/触点 credit/漏斗），全部写成表格字段，再 join 成
``(user_id, item_id)`` 网格行，直接给 sklearn / 下游用。

不做 NetworkX / PageRank / GNN。

Outputs under results/tencent_gr_tabular_ui/:
  - user_features / item_features / edge_ui_features
  - feature_grid.parquet          ← 建模主表（网格）
  - feature_grid_sample.csv
  - localize_convert_path.csv
  - localize_item_covisit.parquet

  PYTHONPATH=. python3 scripts/tencent_gr/run_tabular_user_item_feats.py \\
    --root data/tencent_subset --max-users 3000
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
A_EXP, A_CLK, A_CNV = 0, 1, 2
Event = Tuple[int, int, int]


def _pq(d: Path) -> List[Path]:
    return sorted(p for p in d.glob("*.parquet") if p.is_file())


def parse_events(seq: list) -> List[Event]:
    evs = []
    for e in seq:
        if not isinstance(e, dict):
            e = dict(e)
        evs.append((int(e["item_id"]), int(e["action_type"]), int(e["timestamp"])))
    evs.sort(key=lambda x: x[2])
    return evs


def iter_users(seq_dir: Path, max_users: int):
    n = 0
    for path in _pq(seq_dir):
        pf = pq.ParquetFile(path)
        for rg in range(pf.num_row_groups):
            df = pf.read_row_group(rg).to_pandas()
            for _, row in df.iterrows():
                yield int(row["user_id"]), parse_events(list(row["seq"]))
                n += 1
                if n >= max_users:
                    return


def path_items(evs: Sequence[Event]) -> List[int]:
    out: List[int] = []
    for iid, _, _ in evs:
        if not out or out[-1] != iid:
            out.append(iid)
    return out


def rate(n: float, d: float) -> float:
    return float(n) / float(d) if d > 0 else 0.0


def accumulate(users, *, co_window: int = 8, top_covisit: int = 10):
    """Pure tabular aggregations onto user / item / (u,i) / drill-downs."""
    user_rows: List[dict] = []
    # item accumulators
    item_exp = Counter()
    item_clk = Counter()
    item_cnv = Counter()
    item_users = defaultdict(set)
    item_first_touch = Counter()
    item_last_touch = Counter()
    item_linear = Counter()
    item_as_cnv_terminal = Counter()
    # (u,i)
    ui_exp = Counter()
    ui_clk = Counter()
    ui_cnv = Counter()
    ui_last_ts: Dict[Tuple[int, int], int] = {}
    # localization drill-downs
    convert_path_rows: List[dict] = []
    covisit = defaultdict(Counter)  # item -> Counter(other)

    n_users = 0
    n_conv_users = 0

    for uid, evs in users:
        n_users += 1
        if not evs:
            continue
        n_exp = n_clk = n_cnv = 0
        items_seen = []
        for iid, act, ts in evs:
            items_seen.append(iid)
            item_users[iid].add(uid)
            key = (uid, iid)
            ui_last_ts[key] = ts
            if act == A_EXP:
                n_exp += 1
                item_exp[iid] += 1
                ui_exp[key] += 1
            elif act == A_CLK:
                n_clk += 1
                item_clk[iid] += 1
                ui_clk[key] += 1
            elif act == A_CNV:
                n_cnv += 1
                item_cnv[iid] += 1
                ui_cnv[key] += 1

        uniq_items = list(dict.fromkeys(items_seen))
        user_rows.append(
            {
                "user_id": uid,
                "n_events": len(evs),
                "n_exp": n_exp,
                "n_clk": n_clk,
                "n_cnv": n_cnv,
                "n_uniq_items": len(uniq_items),
                "ctr": rate(n_clk, n_exp),
                "cvr": rate(n_cnv, n_clk),
                "ctcvr": rate(n_cnv, n_exp),
                "has_convert": int(n_cnv > 0),
                "span_sec": int(evs[-1][2] - evs[0][2]) if len(evs) > 1 else 0,
            }
        )

        # co-visit counts within sliding window on path (tabular pair counts only)
        path = path_items(evs)
        for i, a in enumerate(path):
            for b in path[i + 1 : i + 1 + co_window]:
                if a == b:
                    continue
                covisit[a][b] += 1
                covisit[b][a] += 1

        # converting path → touch credit (aggregate to item) + drill-down rows
        if n_cnv > 0:
            n_conv_users += 1
            last_cnv_i = max(i for i, (_, a, _) in enumerate(evs) if a == A_CNV)
            cnv_item = evs[last_cnv_i][0]
            prefix = path_items(evs[: last_cnv_i + 1])
            if not prefix:
                continue
            item_as_cnv_terminal[cnv_item] += 1
            item_first_touch[prefix[0]] += 1.0
            item_last_touch[prefix[-1]] += 1.0
            w = 1.0 / float(len(prefix))
            for pos, iid in enumerate(prefix):
                item_linear[iid] += w
                # localization drill-down: one row per (convert, path item)
                if pos == 0:
                    touch = "first"
                elif pos == len(prefix) - 1:
                    touch = "last"
                else:
                    touch = "mid"
                convert_path_rows.append(
                    {
                        "user_id": uid,
                        "convert_item": cnv_item,
                        "path_item": iid,
                        "path_pos": pos,
                        "path_len": len(prefix),
                        "linear_credit": w,
                        "touch": touch,
                        "is_convert_item": int(iid == cnv_item),
                    }
                )

    # ---- assemble item table ----
    all_items = sorted(
        set(item_exp) | set(item_clk) | set(item_cnv) | set(item_linear) | set(item_users)
    )
    # normalize credit shares
    s_first = sum(item_first_touch.values()) or 1.0
    s_last = sum(item_last_touch.values()) or 1.0
    s_lin = sum(item_linear.values()) or 1.0

    item_rows = []
    for iid in all_items:
        ne, nc, nv = item_exp[iid], item_clk[iid], item_cnv[iid]
        item_rows.append(
            {
                "item_id": iid,
                "n_exp": ne,
                "n_clk": nc,
                "n_cnv": nv,
                "n_users": len(item_users[iid]),
                "ctr": rate(nc, ne),
                "cvr": rate(nv, nc),
                "ctcvr": rate(nv, ne),
                "n_as_convert_terminal": item_as_cnv_terminal[iid],
                "credit_first": float(item_first_touch[iid]),
                "credit_last": float(item_last_touch[iid]),
                "credit_linear": float(item_linear[iid]),
                "share_first": float(item_first_touch[iid]) / s_first,
                "share_last": float(item_last_touch[iid]) / s_last,
                "share_linear": float(item_linear[iid]) / s_lin,
                "n_covisit_neighbors": len(covisit[iid]),
            }
        )
    item_df = pd.DataFrame(item_rows)
    user_df = pd.DataFrame(user_rows)

    # ---- (u,i) edge table (tabular, not a graph API) ----
    edge_rows = []
    for (uid, iid), ne in ui_exp.items():
        edge_rows.append(
            {
                "user_id": uid,
                "item_id": iid,
                "n_exp": ne,
                "n_clk": ui_clk[(uid, iid)],
                "n_cnv": ui_cnv[(uid, iid)],
                "last_ts": ui_last_ts.get((uid, iid), 0),
            }
        )
    # also edges that only clk/cnv without exp key — merge
    for (uid, iid), nc in ui_clk.items():
        if (uid, iid) not in ui_exp:
            edge_rows.append(
                {
                    "user_id": uid,
                    "item_id": iid,
                    "n_exp": 0,
                    "n_clk": nc,
                    "n_cnv": ui_cnv[(uid, iid)],
                    "last_ts": ui_last_ts.get((uid, iid), 0),
                }
            )
    for (uid, iid), nv in ui_cnv.items():
        if (uid, iid) not in ui_exp and (uid, iid) not in ui_clk:
            edge_rows.append(
                {
                    "user_id": uid,
                    "item_id": iid,
                    "n_exp": 0,
                    "n_clk": 0,
                    "n_cnv": nv,
                    "last_ts": ui_last_ts.get((uid, iid), 0),
                }
            )
    edge_df = pd.DataFrame(edge_rows)
    if len(edge_df):
        edge_df["ctr"] = edge_df.apply(lambda r: rate(r["n_clk"], r["n_exp"]), axis=1)
        edge_df["has_convert"] = (edge_df["n_cnv"] > 0).astype(int)

    # ---- covisit drill-down (top-k by count; flat table) ----
    covisit_rows = []
    for a, ctr in covisit.items():
        for b, c in ctr.most_common(top_covisit):
            covisit_rows.append(
                {"item_id": a, "neighbor_item": b, "covisit_count": int(c), "rank": None}
            )
    covisit_df = pd.DataFrame(covisit_rows)
    if len(covisit_df):
        covisit_df["rank"] = (
            covisit_df.groupby("item_id")["covisit_count"]
            .rank(method="first", ascending=False)
            .astype(int)
        )

    convert_df = pd.DataFrame(convert_path_rows)
    meta = {
        "n_users": n_users,
        "n_conv_users": n_conv_users,
        "n_items": len(item_df),
        "n_edges_ui": len(edge_df),
        "n_convert_path_rows": len(convert_df),
        "co_window": co_window,
        "top_covisit": top_covisit,
        "main_dims": ["user_id", "item_id"],
        "note": "No graph/network algorithms — only groupby aggregations into tables.",
    }
    return user_df, item_df, edge_df, convert_df, covisit_df, meta


def enrich_spectrum_cols(user_df: pd.DataFrame, item_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Add lightweight 图谱-style aggregate columns (ranks / log counts). No graph API."""
    u = user_df.copy()
    i = item_df.copy()
    for c in ("n_events", "n_uniq_items", "n_exp", "n_clk", "n_cnv"):
        if c in u.columns:
            u[f"log1p_{c}"] = np.log1p(u[c].astype(float))
    u["user_activity_rank"] = u["n_events"].rank(method="average", ascending=False)
    u["user_ctcvr_rank"] = u["ctcvr"].rank(method="average", ascending=False)

    for c in ("n_exp", "n_clk", "n_cnv", "n_users", "n_covisit_neighbors", "credit_linear"):
        if c in i.columns:
            i[f"log1p_{c}"] = np.log1p(i[c].astype(float))
    i["item_pop_rank"] = i["n_users"].rank(method="average", ascending=False)
    i["item_credit_rank"] = i["share_linear"].rank(method="average", ascending=False)
    i["item_cnv_rank"] = i["n_cnv"].rank(method="average", ascending=False)
    return u, i


def build_feature_grid(
    user_df: pd.DataFrame, item_df: pd.DataFrame, edge_df: pd.DataFrame
) -> pd.DataFrame:
    """Join 图谱特征 onto (user, item) rows → modeling grid / panel."""
    if edge_df is None or len(edge_df) == 0:
        return pd.DataFrame()
    u = user_df.add_prefix("u_").rename(columns={"u_user_id": "user_id"})
    i = item_df.add_prefix("i_").rename(columns={"i_item_id": "item_id"})
    e = edge_df.rename(
        columns={
            "n_exp": "e_n_exp",
            "n_clk": "e_n_clk",
            "n_cnv": "e_n_cnv",
            "last_ts": "e_last_ts",
            "ctr": "e_ctr",
            "has_convert": "y_convert",
        }
    )
    grid = e.merge(u, on="user_id", how="left").merge(i, on="item_id", how="left")
    # simple interaction columns (still tabular)
    grid["e_log1p_exp"] = np.log1p(grid["e_n_exp"].astype(float))
    grid["e_log1p_clk"] = np.log1p(grid["e_n_clk"].astype(float))
    if "u_n_uniq_items" in grid.columns and "i_n_users" in grid.columns:
        grid["ui_pop_mismatch"] = (
            np.log1p(grid["i_n_users"].astype(float))
            - np.log1p(grid["u_n_uniq_items"].astype(float))
        )
    return grid


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=3000)
    ap.add_argument("--co-window", type=int, default=8)
    ap.add_argument("--top-covisit", type=int, default=10)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_tabular_ui",
    )
    args = ap.parse_args()

    print(f"scan ≤{args.max_users} users (图谱特征→网格, no graph algos) ...", flush=True)
    users = iter_users(args.root / "seq", args.max_users)
    user_df, item_df, edge_df, convert_df, covisit_df, meta = accumulate(
        users, co_window=args.co_window, top_covisit=args.top_covisit
    )
    user_df, item_df = enrich_spectrum_cols(user_df, item_df)
    grid = build_feature_grid(user_df, item_df, edge_df)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    user_df.to_parquet(args.out_dir / "user_features.parquet", index=False)
    item_df.to_parquet(args.out_dir / "item_features.parquet", index=False)
    edge_df.to_parquet(args.out_dir / "edge_ui_features.parquet", index=False)
    user_df.to_csv(args.out_dir / "user_features.csv", index=False)
    convert_df.to_csv(args.out_dir / "localize_convert_path.csv", index=False)
    covisit_df.to_parquet(args.out_dir / "localize_item_covisit.parquet", index=False)
    if len(covisit_df):
        top_ids = covisit_df.groupby("item_id").size().nlargest(200).index
        covisit_df[covisit_df["item_id"].isin(top_ids)].to_csv(
            args.out_dir / "localize_item_covisit_sample.csv", index=False
        )

    feat_cols = [
        c
        for c in grid.columns
        if c not in ("user_id", "item_id", "y_convert", "e_last_ts")
    ]
    meta.update(
        {
            "n_grid_rows": int(len(grid)),
            "n_feature_cols": len(feat_cols),
            "feature_cols": feat_cols,
            "label": "y_convert",
            "grid": "(user_id, item_id) panel = edge ⋈ user_图谱特征 ⋈ item_图谱特征",
            "note": "图谱特征=关系聚合列；整合成网格面板；无图算法。",
        }
    )
    if len(grid):
        grid.to_parquet(args.out_dir / "feature_grid.parquet", index=False)
        # stratified-ish sample for eyeballing
        pos = grid[grid["y_convert"] == 1]
        neg = grid[grid["y_convert"] == 0]
        n_pos = min(200, len(pos))
        n_neg = min(800, len(neg))
        sample = pd.concat(
            [
                pos.sample(n_pos, random_state=0) if n_pos else pos,
                neg.sample(n_neg, random_state=0) if n_neg else neg,
            ],
            ignore_index=True,
        )
        sample.to_csv(args.out_dir / "feature_grid_sample.csv", index=False)
        (args.out_dir / "feature_cols.json").write_text(
            json.dumps({"label": "y_convert", "features": feat_cols}, indent=2)
        )

    (args.out_dir / "meta.json").write_text(json.dumps(meta, indent=2))

    top_items = (
        item_df.sort_values("share_linear", ascending=False)
        .head(10)[
            [
                "item_id",
                "n_users",
                "n_cnv",
                "share_linear",
                "share_first",
                "share_last",
                "n_covisit_neighbors",
                "item_credit_rank",
            ]
        ]
        .to_string(index=False)
    )
    (args.out_dir / "README.md").write_text(
        "\n".join(
            [
                "# 图谱特征 → 网格面板（无图算法）",
                "",
                "主维：**user_id / item_id**。",
                "图谱特征 = 漏斗 / 触点 credit / 共现强度 / 活跃度 rank 等聚合列。",
                "网格 = `(user, item)` 行，左连 user 特征、右连 item 特征。",
                "",
                "## 主产出",
                "- `feature_grid.parquet` — 建模网格（label=`y_convert`）",
                "- `feature_grid_sample.csv` / `feature_cols.json`",
                "- `user_features.*` / `item_features.parquet` / `edge_ui_features.parquet`",
                "- `localize_convert_path.csv` — convert 路径下钻",
                "- `localize_item_covisit.parquet` — 共现 top-k 下钻（count only）",
                "",
                f"- grid rows: **{meta.get('n_grid_rows')}**, feature cols: **{meta.get('n_feature_cols')}**",
                "",
                "## Top items by linear credit share",
                "```",
                top_items,
                "```",
                "",
            ]
        )
    )
    print(json.dumps({k: meta[k] for k in meta if k != "feature_cols"}, indent=2))
    print("wrote", args.out_dir)
    print("grid shape:", grid.shape, "label pos rate:", float(grid["y_convert"].mean()) if len(grid) else None)
    print("feature cols:", len(feat_cols))


if __name__ == "__main__":
    main()