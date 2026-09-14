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


# 关联定义。先定 G：节点、时间、边。边 = 业务层，不是随便连。
# 节点 u 人 / i 货 / m 店。单号不是节点（单是粒，挂 Y）。店名不是边。
# 时间：只左窗。曝光不是边。
ASSOC = [
    {
        "name": "catalog",
        "edge": "item—merchant",
        "when": "上架，行为之前",
        "rule": "join item→merchant。货在店上，不是点出来的。",
    },
    {
        "name": "hit",
        "edge": "user—item",
        "when": "左窗 clk 或 cnv，去重",
        "rule": "有过点击或转化才连。曝光不是边。没点就买仍因 cnv 连。",
    },
    {
        "name": "hop",
        "edge": "user—merchant",
        "when": "hit ⋈ catalog",
        "rule": "人到店是 hop，不是新观察。",
    },
    {
        "name": "proj",
        "edge": "merchant—merchant",
        "when": "user—merchant 二部投影",
        "rule": "两个店至少有一个共同买家。不是店名相似。",
    },
    {
        "name": "coclick",
        "edge": "item—item",
        "when": "同场共点，clk only",
        "rule": "相邻 Δt≤30min 一场，一场里点过的货两两连。cnv/exp 不进。",
    },
]


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


def assoc_edges(ev: pd.DataFrame, imap: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """五层关联。传入左窗则无未来泄漏。"""
    cat = imap.drop_duplicates(["item_id", "merchant_id"])[["item_id", "merchant_id"]].copy()
    hit = ev.loc[ev.act.isin([CLK, CNV]), ["user_id", "item_id"]].drop_duplicates()
    hit = hit.merge(cat, on="item_id", how="left")
    ui = hit[["user_id", "item_id"]].drop_duplicates()
    um = hit.dropna(subset=["merchant_id"]).drop_duplicates(["user_id", "merchant_id"])[
        ["user_id", "merchant_id"]
    ]
    pairs = []
    for _, g in um.groupby("user_id"):
        ms = sorted({int(x) for x in g.merchant_id})
        for a, b in combinations(ms, 2):
            pairs.append((a, b))
    proj = (
        pd.DataFrame(pairs, columns=["merchant_a", "merchant_b"]).drop_duplicates()
        if pairs
        else pd.DataFrame(columns=["merchant_a", "merchant_b"])
    )
    clk_only = ev.loc[ev.act.eq(CLK), ["user_id", "item_id", "ts"]]
    co = pd.DataFrame(session_pairs(clk_only), columns=["item_a", "item_b"])
    if len(co):
        co = co.drop_duplicates()
    else:
        co = pd.DataFrame(columns=["item_a", "item_b"])
    return {"catalog": cat, "hit": ui, "hop": um, "proj": proj, "coclick": co}


def _nx_edges(df: pd.DataFrame, s: str, t: str, smap, tmap) -> nx.Graph:
    if df is None or len(df) == 0:
        return nx.Graph()
    return nx.from_pandas_edgelist(
        pd.DataFrame({"s": df[s].map(smap), "t": df[t].map(tmap)}),
        "s",
        "t",
    )


def graph_tables(left: pd.DataFrame, imap: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """左窗点击/转化 → 三张度数表。曝光不连边。边定义见 ASSOC / assoc_edges。"""
    E = assoc_edges(left, imap)
    G_ui = _nx_edges(E["hit"], "user_id", "item_id", _uid, _iid)
    G_um = _nx_edges(E["hop"], "user_id", "merchant_id", _uid, _mid)
    m_nodes = [n for n in G_um.nodes if str(n).startswith("m")]
    G_m = nx.bipartite.projected_graph(G_um, m_nodes) if m_nodes else nx.Graph()
    pr = nx.pagerank(G_m, max_iter=50) if G_m.number_of_nodes() else {}
    clust = nx.clustering(G_m) if G_m.number_of_nodes() else {}

    G_co = nx.Graph()
    G_co.add_edges_from((_iid(a), _iid(b)) for a, b in E["coclick"].itertuples(index=False))

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
    E = assoc_edges(left, imap)
    assert float(u.loc[u.user_id.eq(1), "g_u_item_deg"].iloc[0]) == 2.0
    assert float(i.loc[i.item_id.eq(10), "g_i_user_deg"].iloc[0]) == 2.0
    assert len(E["coclick"]) == 1
    assert set(E["coclick"].iloc[0]) == {10, 11}
    assert len(E["hit"]) == 4  # (1,10)(1,11)(2,10)(3,12)
    print("nx demo", u.to_string(index=False))
    print("assoc", {k: len(v) for k, v in E.items()})
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
        "不是 GNN，不是 DFS。LOGO 把 graph 整包拿掉。",
        "",
        "## 关联定义",
        "",
        "节点 `u` 人 / `i` 货 / `m` 店。单号不是节点。店名不是边。只左窗。曝光不是边。",
        "",
    ]
    for a in ASSOC:
        md.append(f"- **{a['name']}** `{a['edge']}` — {a['when']}。{a['rule']}")
    md += ["", ""]
    md += _md_board("订单粒 漏斗 ∪ graph", rec)
    (OUT / "GRAPH_FSDS.md").write_text("\n".join(md), encoding="utf-8")
    print("wrote", OUT / "GRAPH_FSDS.md")


if __name__ == "__main__":
    main()
