#!/usr/bin/env python3
"""FSDS + networkx 图谱特征融合。PO-risk 没有因果性。

φ=(Y-μ)(W-e) 是早/晚对 Y 的距离，不是 treatment effect。
图谱只是另一包 X：左窗行为 → networkx 度数 / 投影 PageRank / 同场共点。
不是 GNN，不是 DFS。漏斗列仍在；graph 当 hop 给 LOGO/LOCO。

  python3 scripts/tencent_gr/graph_fsds_proto.py
  # pip install networkx
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from feat_proto import split_lr  # noqa: E402
from hop_fsds_proto import ORDER_X, USER_HOP, _md_board, _xy, board  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
EV = ROOT / "results/tencent_gr_fs150/tables/ev.parquet"
POST = ROOT / "results/tencent_gr_fs150/tables/post.parquet"
USER = ROOT / "results/tencent_gr_fs150/tables/user.parquet"
IMAP = ROOT / "results/tencent_gr_fs150/hop/item_merchant.parquet"
OUT = ROOT / "results/tencent_gr_fs150/hop"
SEED = 0
CLK, CNV = 1, 2
SESS_GAP = 30 * 60


def hop_of(name: str) -> str:
    if name.startswith("g_"):
        return "graph"
    if name in USER_HOP:
        return "funnel_user"
    return "funnel_order"


def _uid(x) -> str:
    return f"u{int(x)}"


def _iid(x) -> str:
    return f"i{int(x)}"


def _mid(x) -> str:
    return f"m{int(x)}"


def session_pairs(clk: pd.DataFrame) -> list[tuple[int, int]]:
    """同场共点：一场里点过的货两两连边。"""
    if clk.empty:
        return []
    c = clk.sort_values(["user_id", "ts"]).copy()
    gap = c.groupby("user_id")["ts"].diff().fillna(0)
    c["sess"] = gap.gt(SESS_GAP).groupby(c.user_id).cumsum()
    pairs = []
    for _, g in c.groupby(["user_id", "sess"], sort=False):
        items = g.item_id.unique()
        if len(items) < 2:
            continue
        for a, b in combinations(sorted(int(x) for x in items), 2):
            pairs.append((a, b))
    return pairs


def graph_tables(left: pd.DataFrame, imap: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """左窗点击/转化 → 三张度数表。曝光不连边。"""
    hit = left.loc[left.act.isin([CLK, CNV]), ["user_id", "item_id", "ts"]].copy()
    hit = hit.merge(imap, on="item_id", how="left")
    clk_only = left.loc[left.act.eq(CLK), ["user_id", "item_id", "ts"]]

    ui = hit.drop_duplicates(["user_id", "item_id"])
    G_ui = nx.from_pandas_edgelist(
        pd.DataFrame({"s": ui.user_id.map(_uid), "t": ui.item_id.map(_iid)}),
        "s",
        "t",
    )
    um = hit.dropna(subset=["merchant_id"]).drop_duplicates(["user_id", "merchant_id"])
    G_um = nx.from_pandas_edgelist(
        pd.DataFrame({"s": um.user_id.map(_uid), "t": um.merchant_id.map(_mid)}),
        "s",
        "t",
    )
    m_nodes = [n for n in G_um.nodes if str(n).startswith("m")]
    G_m = nx.bipartite.projected_graph(G_um, m_nodes) if m_nodes else nx.Graph()
    pr = nx.pagerank(G_m, max_iter=50) if G_m.number_of_nodes() else {}
    clust = nx.clustering(G_m) if G_m.number_of_nodes() else {}

    G_co = nx.Graph()
    G_co.add_edges_from((_iid(a), _iid(b)) for a, b in session_pairs(clk_only))

    def deg(G, prefix, key):
        rows = []
        for n, d in G.degree():
            if str(n).startswith(prefix):
                rows.append({key: int(str(n)[1:]), "deg": float(d)})
        return pd.DataFrame(rows)

    users = deg(G_ui, "u", "user_id").rename(columns={"deg": "g_u_item_deg"})
    u2 = deg(G_um, "u", "user_id").rename(columns={"deg": "g_u_merch_deg"})
    users = users.merge(u2, on="user_id", how="outer") if len(u2) else users
    items = deg(G_ui, "i", "item_id").rename(columns={"deg": "g_i_user_deg"})
    i2 = deg(G_co, "i", "item_id").rename(columns={"deg": "g_i_coclick_deg"})
    items = items.merge(i2, on="item_id", how="outer") if len(i2) else items
    merch = deg(G_um, "m", "merchant_id").rename(columns={"deg": "g_m_user_deg"})
    if len(merch):
        merch["g_m_pr"] = merch.merchant_id.map(lambda x: pr.get(_mid(x), 0.0))
        merch["g_m_clust"] = merch.merchant_id.map(lambda x: clust.get(_mid(x), 0.0))
        merch["g_m_proj_deg"] = merch.merchant_id.map(
            lambda x: float(G_m.degree(_mid(x))) if G_m.has_node(_mid(x)) else 0.0
        )
    else:
        merch = pd.DataFrame(
            columns=["merchant_id", "g_m_user_deg", "g_m_pr", "g_m_clust", "g_m_proj_deg"]
        )
    return users.fillna(0.0), items.fillna(0.0), merch.fillna(0.0)


def demo() -> None:
    t0 = 1_700_000_000
    rows = [
        (1, 10, CLK, t0),
        (1, 11, CLK, t0 + 10),
        (2, 10, CLK, t0 + 100),
        (2, 10, CNV, t0 + 120),
        (3, 12, CLK, t0 + 200),
    ]
    left = pd.DataFrame(rows, columns=["user_id", "item_id", "act", "ts"])
    imap = pd.DataFrame({"item_id": [10, 11, 12], "merchant_id": [0, 0, 1]})
    u, i, m = graph_tables(left, imap)
    assert float(u.loc[u.user_id.eq(1), "g_u_item_deg"].iloc[0]) == 2.0
    assert float(i.loc[i.item_id.eq(10), "g_i_user_deg"].iloc[0]) == 2.0
    print("nx demo", u.to_string(index=False))
    print(m.to_string(index=False))


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--no-loco", action="store_true")
    args = ap.parse_args()
    demo()
    if not EV.exists():
        print("no ev, demo only")
        return

    import po_fs_logo

    po_fs_logo.TREES = 20
    po_fs_logo.DEPTH = 5

    ev = pd.read_parquet(EV)
    post = pd.read_parquet(POST)
    users = pd.read_parquet(USER)
    if IMAP.exists():
        imap = pd.read_parquet(IMAP)[["item_id", "merchant_id"]].drop_duplicates()
    else:
        imap = pd.DataFrame(
            {
                "item_id": ev.item_id.unique(),
                "merchant_id": (ev.item_id.unique().astype(np.int64) * (10**9 + 7)) % 400,
            }
        )
    left, _, _ = split_lr(ev)
    print("left", len(left), "clk", int((left.act == CLK).sum()), flush=True)
    gu, gi, gm = graph_tables(left, imap)
    print("graph users", len(gu), "items", len(gi), "shops", len(gm), flush=True)

    o = post.merge(imap, on="item_id", how="left")
    ukeep = [c for c in USER_HOP if c in users.columns]
    o = o.merge(users[["user_id"] + ukeep], on="user_id", how="left")
    o = o.merge(gu, on="user_id", how="left")
    o = o.merge(gi, on="item_id", how="left")
    o = o.merge(gm, on="merchant_id", how="left")
    gcols = [c for c in o.columns if c.startswith("g_")]
    for c in gcols:
        o[c] = o[c].fillna(0.0)
    o["_y"] = pd.to_numeric(o["y_post_clk_1d"], errors="coerce")
    xcols = ukeep + [c for c in ORDER_X if c in o.columns] + gcols
    _, names, X, y, w = _xy(o, xcols, "_y", "t_end", "_y")
    print("=== funnel + graph FSDS (no causality) ===", flush=True)
    rec = board(X, y, w, names, hop_fn=hop_of, seed=SEED, do_loco=not args.no_loco)

    OUT.mkdir(parents=True, exist_ok=True)
    gu.to_parquet(OUT / "graph_user.parquet", index=False)
    gi.to_parquet(OUT / "graph_item.parquet", index=False)
    gm.to_parquet(OUT / "graph_merchant.parquet", index=False)
    (OUT / "GRAPH_FSDS.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    md = [
        "# Funnel FSDS + networkx graph features",
        "",
        "PO-risk **没有因果性**。φ=(Y−μ)(W−e) 是早/晚对买后 Y 的距离，不是谁导致转化。",
        "图谱：左窗 clk/cnv 的 user—item、user—merchant 二部图，店投影 PageRank/聚类，同场共点。",
        "LOGO 把 graph 整包拿掉。不是 GNN，不是 DFS。",
        "",
    ]
    md += _md_board("订单粒 漏斗 ∪ graph", rec)
    (OUT / "GRAPH_FSDS.md").write_text("\n".join(md), encoding="utf-8")
    print("wrote", OUT / "GRAPH_FSDS.md")


if __name__ == "__main__":
    main()
