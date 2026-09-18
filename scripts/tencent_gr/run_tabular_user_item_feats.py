#!/usr/bin/env python3
"""TencentGR tabular user/item features + localization drill-down (NO graphs).

Main dimensions only: **user** and **item**.
Everything else aggregates into columns on these two tables.
Localization = flat drill-down tables (groupby), not network algorithms.

Outputs (under results/tencent_gr_tabular_ui/):
  - user_features.parquet / .csv
  - item_features.parquet / .csv
  - edge_ui_features.parquet      (optional (u,i) aggregates)
  - localize_convert_path.csv    (convert → path items with touch credit)
  - localize_item_covisit.csv    (item → top co-visited items by count)

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

    print(f"scan ≤{args.max_users} users (tabular only) ...", flush=True)
    users = iter_users(args.root / "seq", args.max_users)
    user_df, item_df, edge_df, convert_df, covisit_df, meta = accumulate(
        users, co_window=args.co_window, top_covisit=args.top_covisit
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)

    user_df.to_parquet(args.out_dir / "user_features.parquet", index=False)
    item_df.to_parquet(args.out_dir / "item_features.parquet", index=False)
    edge_df.to_parquet(args.out_dir / "edge_ui_features.parquet", index=False)
    user_df.to_csv(args.out_dir / "user_features.csv", index=False)
    convert_df.to_csv(args.out_dir / "localize_convert_path.csv", index=False)
    covisit_df.to_parquet(args.out_dir / "localize_item_covisit.parquet", index=False)
    if len(covisit_df):
        top_ids = (
            covisit_df.groupby("item_id").size().nlargest(200).index
            if len(covisit_df)
            else []
        )
        covisit_df[covisit_df["item_id"].isin(top_ids)].to_csv(
            args.out_dir / "localize_item_covisit_sample.csv", index=False
        )
    (args.out_dir / "meta.json").write_text(json.dumps(meta, indent=2))

    # small preview markdown
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
            ]
        ]
        .to_string(index=False)
    )
    (args.out_dir / "README.md").write_text(
        "\n".join(
            [
                "# Tabular user / item features + localization drill-down",
                "",
                "Main dims: **user_id**, **item_id**. No NetworkX / PageRank / Markov graph.",
                "",
                "## Tables",
                "- `user_features.*` — one row per user",
                "- `item_features.parquet` — one row per item (funnel + first/last/linear credit shares)",
                "- `edge_ui_features.parquet` — one row per (user, item)",
                "- `localize_convert_path.csv` — drill-down: convert → path items + linear credit",
                "- `localize_item_covisit.parquet` (+ `_sample.csv`) — item → top co-visited items (count)",
                "",
                f"meta: `{json.dumps(meta)}`",
                "",
                "## Top items by linear credit share",
                "```",
                top_items,
                "```",
                "",
            ]
        )
    )
    print(json.dumps(meta, indent=2))
    print("wrote", args.out_dir)
    print("user cols:", list(user_df.columns))
    print("item cols:", list(item_df.columns))


if __name__ == "__main__":
    main()
