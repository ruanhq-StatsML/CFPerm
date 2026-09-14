#!/usr/bin/env python3
"""小合成集：回忆漏斗 + 场 + 店跳 + 左窗构图。防泄漏。

每人一条能点名的故事（uid%6）。店名是标签不是 X。单号只在 cnv_ts 出现。
构图只用左窗 clk/cnv；曝光不是边；右窗不得进 G。

  python3 scripts/tencent_gr/synth_recall.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from feat_proto import split_lr  # noqa: E402
from graph_fsds_proto import graph_tables  # noqa: E402
from merchant_name import generate_catalog_names  # noqa: E402
from onepass_post import onepass_post, user_post  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/tencent_gr_fs150/synth"
EXP, CLK, CNV = 0, 1, 2
SESS = 30 * 60
T0 = 1_700_000_000
SEED = 0
FOLLOW_PAD = 2 * 86400  # 垫 t_end，买后窗跟满；不参与 t_cut


def catalog(n_item: int, n_merch: int, seed: int) -> pd.DataFrame:
    names = generate_catalog_names(n_merch, seed, backend="faker")
    rows = []
    for i in range(n_item):
        mid = i % n_merch
        rows.append({"item_id": i, "merchant_id": mid, "merchant_name": names[mid]})
    return pd.DataFrame(rows)


def mint_orders(cnv: pd.DataFrame, seed: int) -> pd.DataFrame:
    n = len(cnv)
    rng = np.random.default_rng(seed + 9)
    nums = rng.choice(np.arange(100_000_000, 200_000_000), size=n, replace=False)
    out = cnv.copy()
    out["order_number"] = nums.astype(str)
    return out


def _user_story(uid: int, items: np.ndarray) -> tuple[list, str]:
    """每人一条可点名的漏斗故事。pad 另加，不写进故事。"""
    t = T0 + uid * 86400
    a, b = int(items[0]), int(items[1])
    kind = uid % 6
    evs: list[tuple[int, int, int, int]] = []

    def add(iid, act, dt):
        nonlocal t
        t += dt
        evs.append((uid, iid, act, t))

    if kind == 0:
        tag = "bounce_exp"
        add(a, EXP, 0)
        add(b, EXP, SESS + 60)  # 下一场，看一眼走
    elif kind == 1:
        tag = "clk_no_cnv"
        add(a, EXP, 0)
        add(a, CLK, 20)
        add(b, EXP, 400)
        add(b, CLK, 15)
    elif kind == 2:
        tag = "same_item_cnv"  # 同品点完买 + 当场续逛
        add(a, EXP, 0)
        add(a, CLK, 20)
        add(a, CNV, 600)
        add(a, CLK, 120)
    elif kind == 3:
        tag = "other_clk_cnv"  # 先点别的再买这件 + 跨场回访
        add(b, EXP, 0)
        add(b, CLK, 20)
        add(a, EXP, 300)
        add(a, CNV, 40)
        add(b, CLK, SESS + 90)
    elif kind == 4:
        tag = "empty_path_cnv"  # 没点就买
        add(a, EXP, 0)
        add(a, CNV, 80)
    else:
        tag = "two_cnv"
        add(a, EXP, 0)
        add(a, CLK, 10)
        add(a, CNV, 50)
        add(b, EXP, SESS + 30)
        add(b, CLK, 20)
        add(b, CNV, 30)
        add(b, CLK, 200)
    return evs, tag


def make_synth(n_users: int = 60, n_item: int = 24, n_merch: int = 6, seed: int = SEED):
    rng = np.random.default_rng(seed)
    cat = catalog(n_item, n_merch, seed)
    stories = []
    rows = []
    for uid in range(n_users):
        pair = rng.choice(cat.item_id.to_numpy(), size=2, replace=False)
        evs, tag = _user_story(uid, pair)
        stories.append(
            {"user_id": uid, "story": tag, "item_a": int(pair[0]), "item_b": int(pair[1])}
        )
        rows.extend(evs)
    story = pd.DataFrame(rows, columns=["user_id", "item_id", "act", "ts"])
    last = story.groupby("user_id")["ts"].max().rename("t_story")
    pad = last.reset_index()
    pad["item_id"] = 0
    pad["act"] = EXP
    pad["ts"] = pad["t_story"] + FOLLOW_PAD
    pad = pad[["user_id", "item_id", "act", "ts"]]
    ev = pd.concat([story.assign(is_pad=False), pad.assign(is_pad=True)], ignore_index=True)
    ev = ev.merge(cat, on="item_id", how="left")
    story = story.merge(cat, on="item_id", how="left")
    cnv = story.loc[story.act.eq(CNV), ["user_id", "item_id", "merchant_id", "ts"]].rename(
        columns={"ts": "cnv_ts"}
    )
    orders = mint_orders(cnv, seed)
    orders = orders.merge(cat[["item_id", "merchant_name"]], on="item_id", how="left")
    return ev, story, cat, orders, pd.DataFrame(stories)


def path_mix(ev: pd.DataFrame) -> pd.DataFrame:
    """转化路径：同品点过 / 只点过别的 / 空路径。"""
    rows = []
    for _, g in ev.groupby("user_id", sort=False):
        g = g.sort_values("ts")
        for _, r in g.loc[g.act.eq(CNV)].iterrows():
            before = g.loc[g.ts < r.ts]
            clk = before.loc[before.act.eq(CLK)]
            same = bool((clk.item_id == r.item_id).any())
            any_clk = len(clk) > 0
            if same:
                mode = "A_same_clk"
            elif any_clk:
                mode = "B_other_clk"
            else:
                mode = "D_empty"
            rows.append(mode)
    return pd.Series(rows).value_counts().rename_axis("path").reset_index(name="n")


def hit_pairs(df: pd.DataFrame) -> set[tuple[int, int]]:
    h = df.loc[df.act.isin([CLK, CNV]), ["user_id", "item_id"]]
    return {(int(u), int(i)) for u, i in zip(h.user_id, h.item_id)}


def _post_rows(ev: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for uid, g in ev.groupby("user_id", sort=False):
        seq = [(int(r.item_id), int(r.act), int(r.ts), None) for r in g.itertuples()]
        for rec in onepass_post(seq):
            rec = {k: v for k, v in rec.items() if not isinstance(v, dict)}
            rec["user_id"] = int(uid)
            rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ev, story, cat, orders, stories = make_synth()
    # t_cut 只看故事，pad 一律进右窗（跟满，不当行为）
    left, right, mm = split_lr(story)
    right = pd.concat([right, ev.loc[ev.is_pad, right.columns]], ignore_index=True)
    imap = cat[["item_id", "merchant_id"]]
    gu, gi, gm = graph_tables(left, imap)
    gu_all, _, _ = graph_tables(story, imap)  # 反例：故事全量构图会吃到右窗

    left_p = hit_pairs(left)
    right_only = hit_pairs(right) - left_p
    # 左窗图的 user—item 度 = 左窗 clk/cnv 去重货数
    deg = left.loc[left.act.isin([CLK, CNV])].drop_duplicates(["user_id", "item_id"])
    expect = deg.groupby("user_id").size().rename("expect")
    got = gu.set_index("user_id")["g_u_item_deg"]
    chk = expect.to_frame().join(got, how="outer").fillna(0.0)
    assert (chk["expect"] == chk["g_u_item_deg"]).all()
    leak_n = len(right_only)

    post = _post_rows(ev)
    up = []
    for uid, g in ev.groupby("user_id", sort=False):
        seq = [(int(r.item_id), int(r.act), int(r.ts), None) for r in g.itertuples()]
        row = user_post(onepass_post(seq))
        row["user_id"] = int(uid)
        up.append(row)
    userp = pd.DataFrame(up).merge(stories, on="user_id")

    ev.to_parquet(OUT / "ev.parquet", index=False)
    cat.to_parquet(OUT / "item_merchant.parquet", index=False)
    orders.to_parquet(OUT / "orders.parquet", index=False)
    stories.to_parquet(OUT / "stories.parquet", index=False)
    gu.to_parquet(OUT / "graph_user.parquet", index=False)
    gi.to_parquet(OUT / "graph_item.parquet", index=False)
    gm.to_parquet(OUT / "graph_merchant.parquet", index=False)
    mm.to_parquet(OUT / "split.parquet", index=False)
    post.to_parquet(OUT / "post.parquet", index=False)
    userp.to_parquet(OUT / "user_post.parquet", index=False)

    print("wrote", OUT)
    print("n_ev", len(ev), "n_story", len(story), "n_cnv", int((story.act == CNV).sum()),
          "n_orders", len(orders), "n_shop", cat.merchant_id.nunique())
    print("left", len(left), "right", len(right),
          "graph_users", len(gu), "graph_items", len(gi), "graph_shops", len(gm))
    print("future_clkcnv_pairs", leak_n, "(right-only; must not enter G)")
    print("naive_full_graph extra user-deg",
          float((gu_all.set_index("user_id").g_u_item_deg
                 - gu.set_index("user_id").g_u_item_deg.reindex(gu_all.user_id).fillna(0)).clip(lower=0).sum()))
    print("stories", stories.story.value_counts().sort_index().to_dict())
    print("cnv paths\n", path_mix(story).to_string(index=False))
    print("same_sess vs cross (user_post by story)")
    print(userp.groupby("story")[["post_clk_same_sess_rate", "post_clk_cross_sess_rate"]]
          .mean().round(3).to_string())
    print("orders e.g.", orders.order_number.head(3).tolist())
    print("shops\n",
          cat.drop_duplicates("merchant_id")[["merchant_id", "merchant_name"]].head(6).to_string(index=False))

    u2 = story.loc[story.user_id.eq(2), ["ts", "item_id", "act"]].copy()
    u2["act"] = u2.act.map({0: "exp", 1: "clk", 2: "cnv"})
    print("user 2 (same_item_cnv)\n", u2.to_string(index=False))
    p2 = post.loc[post.user_id.eq(2)]
    if len(p2):
        print("post dt_next", int(p2.iloc[0].dt_next_clk_sec),
              "same_sess", p2.iloc[0].next_clk_same_sess)
    u3 = story.loc[story.user_id.eq(3), ["ts", "item_id", "act"]].copy()
    u3["act"] = u3.act.map({0: "exp", 1: "clk", 2: "cnv"})
    print("user 3 (other_clk_cnv)\n", u3.to_string(index=False))
    p3 = post.loc[post.user_id.eq(3)]
    if len(p3):
        print("post dt_next", int(p3.iloc[0].dt_next_clk_sec),
              "same_sess", p3.iloc[0].next_clk_same_sess,
              "cross_sess", p3.iloc[0].next_clk_cross_sess)


if __name__ == "__main__":
    main()
