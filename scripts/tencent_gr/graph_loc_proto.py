#!/usr/bin/env python3
"""图 localization 原型：关联口径对齐 + 层内 Y 残差分数 + 多种 loc。

不做 PO 上图，不做因果。先回答：边层口径是否和漏斗/左窗一致；
残差分钉在哪层、哪块支撑集。

方法（对照）：
  1) layer_mass     五层分数质量份额（必做基线）
  2) lap_smooth     (L+λI)z=s 拉普拉斯光滑
  3) ppr            种子个性化 PageRank（按层）
  4) local_spectral 两步：PPR 球 → 诱导子图 Fiedler 割
  5) score_scan     种子贪心扩张（分数/√体积）

  python3 scripts/tencent_gr/graph_loc_proto.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from feat_proto import split_lr  # noqa: E402
from graph_fsds_proto import ASSOC, CLK, CNV, SESS_GAP, assoc_edges  # noqa: E402
from onepass_post import onepass_post  # noqa: E402
from synth_recall import EXP, make_synth  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/tencent_gr_fs150/synth/loc"

ALIGN = [
    ("time", "只左窗构图；pad/右窗不得进边"),
    ("expose", "曝光不是边；hit/coclick 不含 exp"),
    ("coclick", f"同场共点 clk-only；Δt≤{SESS_GAP}s 切场"),
    ("hop", "user—merchant = hit ⋈ catalog，不是新观察"),
    ("proj", "店—店 = 共同买家投影，不是店名相似"),
    ("y_grain", "层内 Y 在订单粒；节点分 = 残差聚合，不是把 Y 当边"),
    ("residual", "先减层内可观测结构（漏斗量），再上图 → 看残余共变落点"),
    ("loss", "对齐到店/货粒会抹掉订单内路径 A/B/C/D 与当场/跨场"),
]


def print_align(E_left: dict, E_story: dict, left: pd.DataFrame) -> None:
    print("—— 关联口径对齐 ——")
    for k, msg in ALIGN:
        print(f"  [{k}] {msg}")
    for a in ASSOC:
        print(
            f"  {a['name']:8} {a['edge']:18} "
            f"left={len(E_left[a['name']]):3} story={len(E_story[a['name']]):3}  | {a['when']}"
        )
    left_ui = set(
        zip(
            left.loc[left.act.isin([CLK, CNV]), "user_id"].astype(int),
            left.loc[left.act.isin([CLK, CNV]), "item_id"].astype(int),
        )
    )
    for u, i in zip(E_left["hit"].user_id.astype(int), E_left["hit"].item_id.astype(int)):
        assert (u, i) in left_ui
    print("  ok: left hit ⊆ left clk|cnv")


def order_layer_y(ev: pd.DataFrame) -> pd.DataFrame:
    """转化粒：y=买后当场续逛；跟不满不当 0。"""
    rows = []
    for uid, g in ev.groupby("user_id", sort=False):
        g = g.sort_values("ts")
        seq = [(int(r.item_id), int(r.act), int(r.ts), None) for r in g.itertuples()]
        for rec in onepass_post(seq):
            y = rec.get("next_clk_same_sess")
            if y is None or y != y:
                continue
            cnv_ts = int(rec["cnv_ts"])
            before = g.loc[g.ts < cnv_ts]
            rows.append(
                {
                    "user_id": int(uid),
                    "item_id": int(rec["item_id"]),
                    "cnv_ts": cnv_ts,
                    "y_same_sess": float(y),
                    "n_clk_before": int((before.act == CLK).sum()),
                    "n_exp_before": int((before.act == EXP).sum()),
                    "n_cnv_prior": int((before.act == CNV).sum()),
                }
            )
    return pd.DataFrame(rows)


def residualize(df: pd.DataFrame, ycol: str, xcols: list[str]) -> pd.Series:
    y = df[ycol].to_numpy(dtype=float)
    X = np.column_stack([np.ones(len(df)), df[xcols].to_numpy(dtype=float)])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return pd.Series(y - X @ beta, index=df.index, name="resid")


def node_scores(orders: pd.DataFrame, cat: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    o = orders.merge(cat[["item_id", "merchant_id"]], on="item_id", how="left")
    o["resid"] = residualize(
        o, "y_same_sess", ["n_clk_before", "n_exp_before", "n_cnv_prior"]
    )
    return (
        o.groupby("item_id")["resid"].mean().rename("score"),
        o.groupby("merchant_id")["resid"].mean().rename("score"),
    )


def graph_from_edges(df: pd.DataFrame, a: str, b: str) -> nx.Graph:
    G = nx.Graph()
    if df is None or len(df) == 0:
        return G
    G.add_edges_from((int(u), int(v)) for u, v in zip(df[a], df[b]))
    return G


def laplacian_smooth(G: nx.Graph, score: pd.Series, lam: float = 0.5) -> pd.Series:
    if G.number_of_nodes() == 0:
        return score.astype(float).copy()
    nodes = list(G.nodes())
    s = np.array([float(score.get(n, 0.0)) for n in nodes], dtype=float)
    L = nx.laplacian_matrix(G, nodelist=nodes).astype(float).toarray()
    z = np.linalg.solve(L + lam * np.eye(len(nodes)), s)
    return pd.Series(z, index=nodes, name="lap")


def ppr_from_seeds(
    G: nx.Graph, seeds: list[int], alpha: float = 0.15, topk: int = 8
) -> pd.Series:
    if G.number_of_nodes() == 0 or not seeds:
        return pd.Series(dtype=float, name="ppr")
    personal = {n: 0.0 for n in G.nodes()}
    w = 1.0 / len(seeds)
    for s in seeds:
        if s in personal:
            personal[s] = w
    if sum(personal.values()) <= 0:
        return pd.Series(dtype=float, name="ppr")
    pr = nx.pagerank(G, alpha=1.0 - alpha, personalization=personal, max_iter=100)
    return pd.Series(pr, name="ppr").sort_values(ascending=False).head(topk)


def local_spectral(G: nx.Graph, seeds: list[int], ball: int = 12) -> list[int]:
    pr = ppr_from_seeds(G, seeds, topk=ball)
    if len(pr) < 3:
        return list(pr.index)
    H = G.subgraph(list(pr.index)).copy()
    if H.number_of_edges() == 0:
        return list(pr.index)
    try:
        fiedler = nx.fiedler_vector(H, weight=None)
    except Exception:
        return list(pr.index)
    nodes = list(H.nodes())
    side = [n for n, v in zip(nodes, fiedler) if v >= 0]
    other = [n for n, v in zip(nodes, fiedler) if v < 0]
    seed_set = set(seeds)
    if sum(n in seed_set for n in side) >= sum(n in seed_set for n in other):
        return side
    return other


def score_scan(G: nx.Graph, score: pd.Series, seed: int, budget: int = 6) -> list[int]:
    if seed not in G:
        return [seed] if seed in score.index else []
    S = {seed}

    def util(T: set[int]) -> float:
        sc = sum(float(score.get(x, 0.0)) for x in T)
        vol = sum(dict(G.degree(T)).values()) + 1
        return sc / np.sqrt(vol)

    base = util(S)
    while len(S) < budget:
        border: set[int] = set()
        for u in S:
            border |= set(G.neighbors(u))
        border -= S
        if not border:
            break
        v = max(border, key=lambda x: util(S | {x}))
        if util(S | {v}) < base * 0.05 and len(S) > 1:
            break
        S.add(v)
    return list(S)


def top_support(score: pd.Series, k: int = 5) -> list[int]:
    if score.empty:
        return []
    return list(score.abs().sort_values(ascending=False).head(k).index)


def concentration(score: pd.Series, support: list[int]) -> dict:
    tot = float(score.abs().sum()) + 1e-12
    mass = float(score.reindex(support).abs().fillna(0).sum()) if support else 0.0
    return {"n": len(support), "mass_share": round(mass / tot, 4)}


def layer_mass(E: dict, item_s: pd.Series, merch_s: pd.Series) -> pd.DataFrame:
    rows = []
    specs = [
        ("catalog", E["catalog"], "item_id", "merchant_id", item_s, merch_s),
        ("hit", E["hit"], "item_id", None, item_s, None),
        ("hop", E["hop"], "merchant_id", None, merch_s, None),
        ("proj", E["proj"], "merchant_a", "merchant_b", merch_s, merch_s),
        ("coclick", E["coclick"], "item_a", "item_b", item_s, item_s),
    ]
    for name, df, c1, c2, s1, s2 in specs:
        if df is None or len(df) == 0:
            rows.append((name, 0.0))
            continue
        v = df[c1].map(s1).abs().fillna(0.0)
        if c2 is not None and s2 is not None:
            v = v + df[c2].map(s2).abs().fillna(0.0)
        rows.append((name, float(v.sum())))
    out = pd.DataFrame(rows, columns=["layer", "mass"])
    out["share"] = out["mass"] / (out["mass"].sum() + 1e-12)
    return out.sort_values("share", ascending=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ev, story, cat, _orders, _stories = make_synth()
    left, _right, _mm = split_lr(story)

    oy = order_layer_y(ev)
    if oy.empty:
        raise SystemExit("no observable same-sess y")
    item_s, merch_s = node_scores(oy, cat)
    xcols = ["n_clk_before", "n_exp_before", "n_cnv_prior"]
    resid = residualize(oy, "y_same_sess", xcols)

    imap = cat[["item_id", "merchant_id"]]
    E_left = assoc_edges(left, imap)
    E_story = assoc_edges(story, imap)
    print_align(E_left, E_story, left)

    print("\n—— 层内 Y 残差（买后当场续逛 | OLS ← clk/exp/prior_cnv）——")
    print(
        f"n_orders_obs={len(oy)}  y_mean={oy.y_same_sess.mean():.3f}  "
        f"resid_std={float(resid.std()):.4f}  "
        f"n_item={int(item_s.notna().sum())} n_merch={int(merch_s.notna().sum())}"
    )

    mass = layer_mass(E_left, item_s, merch_s)
    print("\n—— layer_mass（残余 |score| 沿边摊到层）——")
    print(mass.to_string(index=False))

    G_proj = graph_from_edges(E_left["proj"], "merchant_a", "merchant_b")
    G_co = graph_from_edges(E_left["coclick"], "item_a", "item_b")
    seeds_m = top_support(merch_s, 3)
    seeds_i = top_support(item_s, 3)

    print("\n—— localization 对照（店=proj，货=coclick）——")
    z_m = laplacian_smooth(G_proj, merch_s.fillna(0.0), lam=0.5)
    z_i = laplacian_smooth(G_co, item_s.fillna(0.0), lam=0.5)
    supp_m_lap = top_support(z_m, 5)
    supp_i_lap = top_support(z_i, 5)
    print("lap_smooth merch", supp_m_lap, concentration(merch_s, supp_m_lap))
    print("lap_smooth item ", supp_i_lap, concentration(item_s, supp_i_lap))

    ppr_m = ppr_from_seeds(G_proj, seeds_m, topk=5)
    ppr_i = ppr_from_seeds(G_co, seeds_i, topk=5)
    print(
        "ppr merch seeds", seeds_m, "→", list(ppr_m.index),
        concentration(merch_s, list(ppr_m.index)),
    )
    print(
        "ppr item  seeds", seeds_i, "→", list(ppr_i.index),
        concentration(item_s, list(ppr_i.index)),
    )

    ls_m = local_spectral(G_proj, seeds_m, ball=8) if seeds_m else []
    ls_i = local_spectral(G_co, seeds_i, ball=8) if seeds_i else []
    print("local_spectral merch", ls_m, concentration(merch_s, ls_m))
    print("local_spectral item ", ls_i, concentration(item_s, ls_i))

    scan_m = (
        score_scan(G_proj, merch_s, seeds_m[0], budget=6)
        if seeds_m and G_proj.number_of_nodes()
        else []
    )
    scan_i = (
        score_scan(G_co, item_s, seeds_i[0], budget=6)
        if seeds_i and G_co.number_of_nodes()
        else []
    )
    print("score_scan merch", scan_m, concentration(merch_s, scan_m))
    print("score_scan item ", scan_i, concentration(item_s, scan_i))

    print("\n—— 对齐丢信息 + 其他方面（原型显式承认）——")
    print("  grain: 订单→店/货均值抹平 A/B/C/D、当场/跨场")
    print("  multi-order: 多单用户残差进同一店点")
    print("  empty-layer: 左窗无边时 loc 退回种子，不报假支撑")
    print("  collinear layers: hop⊂hit 几何，layer_mass 不能当独立贡献")
    print("  sign: 残差可正可负；mass 用绝对值，方向要另报")
    print("  coverage: 无转化用户不上节点分，bounce 只影响构图不进 Y")

    rec = {
        "layer_mass": mass.to_dict(orient="records"),
        "seeds_merch": seeds_m,
        "seeds_item": seeds_i,
        "lap_merch": supp_m_lap,
        "lap_item": supp_i_lap,
        "ppr_merch": list(ppr_m.index),
        "ppr_item": list(ppr_i.index),
        "local_spectral_merch": ls_m,
        "local_spectral_item": ls_i,
        "scan_merch": scan_m,
        "scan_item": scan_i,
        "proj_nodes": G_proj.number_of_nodes(),
        "proj_edges": G_proj.number_of_edges(),
        "coclick_nodes": G_co.number_of_nodes(),
        "coclick_edges": G_co.number_of_edges(),
        "n_oy": len(oy),
    }
    oy.assign(resid=resid).to_parquet(OUT / "order_y_resid.parquet", index=False)
    item_s.to_frame().to_parquet(OUT / "score_item.parquet")
    merch_s.to_frame().to_parquet(OUT / "score_merch.parquet")
    mass.to_parquet(OUT / "layer_mass.parquet", index=False)
    (OUT / "loc_summary.json").write_text(
        json.dumps(rec, indent=2, default=str), encoding="utf-8"
    )
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
