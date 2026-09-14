#!/usr/bin/env python3
"""好几年不写就忘。跑一遍，对着六条故事把口径找回来。

  python3 scripts/tencent_gr/synth_recall.py

时间不能倒。一层漏斗一问，不要共用一个 Y。
曝光 → 点击 → 转化 → 买后点击
 CTR     CVR    路径     续逛 ≠ 转化

相邻 Δt>30min 切断一场。同一 W：(t−W,t] 是 X，(t,t+W] 是 Y。
跟不满不当 0。空 asof 不填 0。n_clk 是量，cnv_share 是结构。

店名是标签不是 X。先有货在店上，行为序列，cnv_ts 才有单号。
图：左窗 clk/cnv 才是边，曝光不是边，右窗不得进 G。不是 GNN，不是 DFS。
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

# uid%6。忘了就从这六个人往下对。
STORY = {
    0: "bounce_exp",       # 看一眼走。有场切，没点。曝光层。图上没边。
    1: "clk_no_cnv",       # 点了不买。有 CTR，没有 CVR，没有单。
    2: "same_item_cnv",    # 同品点完买。路径 A。当场续逛 ≠ 第二笔转化。
    3: "other_clk_cnv",    # 先点别的再买这件。路径 B。跨场回访。
    4: "empty_path_cnv",   # 没点就买。路径 D。dt_item 空，别填 0。
    5: "two_cnv",          # 两单。第一单 next 跨场，第二单当场。
}


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
        add(a, EXP, 0)
        add(b, EXP, SESS + 60)
    elif kind == 1:
        add(a, EXP, 0)
        add(a, CLK, 20)
        add(b, EXP, 400)
        add(b, CLK, 15)
    elif kind == 2:
        add(a, EXP, 0)
        add(a, CLK, 20)
        add(a, CNV, 600)
        add(a, CLK, 120)
    elif kind == 3:
        add(b, EXP, 0)
        add(b, CLK, 20)
        add(a, EXP, 300)
        add(a, CNV, 40)
        add(b, CLK, SESS + 90)
    elif kind == 4:
        add(a, EXP, 0)
        add(a, CNV, 80)
    else:
        add(a, EXP, 0)
        add(a, CLK, 10)
        add(a, CNV, 50)
        add(b, EXP, SESS + 30)
        add(b, CLK, 20)
        add(b, CNV, 30)
        add(b, CLK, 200)
    return evs, STORY[kind]


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


def _fmt_seq(g: pd.DataFrame) -> str:
    g = g.sort_values("ts")
    prev = None
    bits = []
    actn = {0: "exp", 1: "clk", 2: "cnv"}
    for r in g.itertuples():
        dt = "" if prev is None else f"+{int(r.ts - prev)}s"
        if prev is not None and r.ts - prev > SESS:
            dt += "|sess"
        bits.append(f"{dt} i{int(r.item_id)}:{actn[int(r.act)]}".strip())
        prev = r.ts
    return " ".join(bits)


def walk(story: pd.DataFrame, post: pd.DataFrame) -> None:
    print("—— 六条故事（user 0..5）——")
    for uid in range(6):
        g = story.loc[story.user_id.eq(uid)]
        print(f"u{uid} {STORY[uid]}: {_fmt_seq(g)}")
        rows = post.loc[post.user_id.eq(uid)]
        if rows.empty:
            print("    无单")
            continue
        for i, r in enumerate(rows.itertuples()):
            dt_item = r.dt_item_min
            dt_item_s = "empty" if dt_item != dt_item else f"{dt_item:.1f}min"
            nxt = r.dt_next_clk_sec
            nxt_s = "none" if nxt is None or nxt != nxt else f"{int(nxt)}s"
            print(
                f"    cnv{i} item={int(r.item_id)} wo_prior_clk={bool(r.wo_prior_clk)} "
                f"dt_item={dt_item_s} next={nxt_s} "
                f"same={r.next_clk_same_sess} cross={r.next_clk_cross_sess}"
            )


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
    print(__doc__)
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
    print("shops", cat.drop_duplicates("merchant_id").merchant_name.tolist())
    walk(story, post)

    p2 = post.loc[post.user_id.eq(2)]
    p3 = post.loc[post.user_id.eq(3)]
    p4 = post.loc[post.user_id.eq(4)]
    p5 = post.loc[post.user_id.eq(5)]
    assert float(p2.iloc[0].next_clk_same_sess) == 1.0
    assert float(p3.iloc[0].next_clk_cross_sess) == 1.0
    assert bool(p4.iloc[0].wo_prior_clk)
    assert p4.iloc[0].dt_item_min != p4.iloc[0].dt_item_min  # empty ≠ 0
    assert float(p5.iloc[0].next_clk_cross_sess) == 1.0
    assert float(p5.iloc[1].next_clk_same_sess) == 1.0
    bounce = story.loc[story.user_id.eq(0)].sort_values("ts")
    assert int(bounce.act.max()) == EXP
    assert int(bounce.ts.diff().iloc[-1]) > SESS
    assert 0 not in set(gu.user_id.astype(int))  # 曝光不是边
    empty = path_mix(story)
    assert int(empty.loc[empty.path.eq("D_empty"), "n"].iloc[0]) > 0


if __name__ == "__main__":
    main()
