#!/usr/bin/env python3
"""Block tables → polish → merge. Attribution is merge_asof, not a dict walk.

Grain
-----
  ev     user × event
  attr   user × conversion   (asof 中间表，归因在这里细化)
  sess   user × session
  user_* 每块聚合成 user 一行，再 outer join

  python3 scripts/tencent_gr/block_tables.py
"""
from __future__ import annotations

import math
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

EXP, CLK, CNV = 0, 1, 2
WINS = {"3d": 3 * 86400, "7d": 7 * 86400, "14d": 14 * 86400, "30d": 30 * 86400}
HLS = {"hl1d": 86400.0, "hl3d": 3 * 86400.0, "hl7d": 7 * 86400.0}
SESS_GAP = 30 * 60
ATTR_BUCKETS = [(5, "5m"), (30, "30m"), (60, "1h"), (1440, "1d")]
Event = Tuple[int, int, int, Optional[float]]


def rate(n, d):
    return float(n) / float(d) if d > 0 else 0.0


def events_frame(uid: int, evs: List[Event]) -> pd.DataFrame:
    if not evs:
        return pd.DataFrame(columns=["user_id", "item_id", "act", "ts", "price"])
    evs = sorted(evs, key=lambda e: e[2])
    return pd.DataFrame(
        {
            "user_id": uid,
            "item_id": [e[0] for e in evs],
            "act": [e[1] for e in evs],
            "ts": [e[2] for e in evs],
            "price": [e[3] for e in evs],
        }
    )


def add_age(ev: pd.DataFrame) -> pd.DataFrame:
    t_end = ev.groupby("user_id")["ts"].transform("max")
    out = ev.copy()
    out["t_end"] = t_end
    out["age"] = (t_end - out["ts"]).clip(lower=0)
    return out


# ----- 块：漏斗 -----
def tab_funnel(ev: pd.DataFrame) -> pd.DataFrame:
    ev = add_age(ev)
    rows = []
    for uid, g in ev.groupby("user_id", sort=False):
        rec = {"user_id": uid, "hist_len": float(len(g)), "active_days": float(g["ts"].div(86400).astype(int).nunique())}
        for w, sec in {"life": 10**18, **WINS}.items():
            sub = g if w == "life" else g[g["age"] < sec]
            n_exp = int((sub.act == EXP).sum())
            n_clk = int((sub.act == CLK).sum())
            n_cnv = int((sub.act == CNV).sum())
            tot = int(len(sub))
            rec[f"{w}_n_exp"] = n_exp
            rec[f"{w}_n_clk"] = n_clk
            rec[f"{w}_n_cnv"] = n_cnv
            rec[f"{w}_ctr"] = rate(n_clk, n_exp)
            rec[f"{w}_cvr"] = rate(n_cnv, n_clk)
            rec[f"{w}_ctcvr"] = rate(n_cnv, n_exp)
            rec[f"{w}_clk_share"] = rate(n_clk, tot)
            rec[f"{w}_cnv_share"] = rate(n_cnv, tot)
        rec["pay_cnt"] = rec["life_n_cnv"]
        rec["trend_n_cnv_3d_minus_14d"] = rec["3d_n_cnv"] - rec["14d_n_cnv"]
        rec["trend_cvr_3d_minus_14d"] = rec["3d_cvr"] - rec["14d_cvr"]
        rows.append(rec)
    return pd.DataFrame(rows)


# ----- 块：衰减 -----
def tab_decay(ev: pd.DataFrame) -> pd.DataFrame:
    ev = add_age(ev)
    rows = []
    for uid, g in ev.groupby("user_id", sort=False):
        rec = {"user_id": uid}
        for h, sec in HLS.items():
            lam = math.log(2.0) / sec
            w = np.exp(-lam * g["age"].to_numpy())
            rec[f"dec_{h}_dec_clk"] = float(w[g.act.eq(CLK).to_numpy()].sum())
            rec[f"dec_{h}_dec_cnv"] = float(w[g.act.eq(CNV).to_numpy()].sum())
        rows.append(rec)
    return pd.DataFrame(rows)


# ----- 块：场次 -----
def tab_session(ev: pd.DataFrame) -> pd.DataFrame:
    ev = ev.sort_values(["user_id", "ts"]).copy()
    gap = ev.groupby("user_id")["ts"].diff().fillna(0)
    ev["sess"] = gap.gt(SESS_GAP).groupby(ev["user_id"]).cumsum()
    s = ev.groupby(["user_id", "sess"], sort=False).agg(
        n=("act", "size"),
        n_clk=("act", lambda a: int((a == CLK).sum())),
        n_cnv=("act", lambda a: int((a == CNV).sum())),
    )
    u = s.groupby("user_id").agg(
        sess_n=("n", "size"),
        sess_bounce_rate=("n", lambda x: float((x <= 1).mean())),
        sess_depth_clk_mean=("n_clk", "mean"),
        sess_depth_cnv_mean=("n_cnv", "mean"),
        sess_clk_sess_rate=("n_clk", lambda x: float((x > 0).mean())),
        sess_cnv_sess_rate=("n_cnv", lambda x: float((x > 0).mean())),
    )
    return u.reset_index()


# ----- 块：归因（中间表 = 每笔转化一行）-----
def tab_attr_events(ev: pd.DataFrame) -> pd.DataFrame:
    """每笔 cnv 一行。同商品 last-click / 任意 last-click / 同商品 first-click。"""
    clk = (
        ev.loc[ev.act == CLK, ["user_id", "item_id", "ts"]]
        .rename(columns={"ts": "clk_ts"})
        .sort_values(["user_id", "item_id", "clk_ts"])
    )
    cnv = (
        ev.loc[ev.act == CNV, ["user_id", "item_id", "ts"]]
        .rename(columns={"ts": "cnv_ts"})
        .sort_values(["user_id", "item_id", "cnv_ts"])
    )
    if cnv.empty:
        return cnv.assign(
            dt_item_min=np.nan,
            dt_any_min=np.nan,
            dt_first_min=np.nan,
            wo_prior_clk=True,
        )

    same = pd.merge_asof(
        cnv,
        clk,
        by=["user_id", "item_id"],
        left_on="cnv_ts",
        right_on="clk_ts",
        direction="backward",
    )
    clk_any = clk.drop(columns="item_id").sort_values(["user_id", "clk_ts"])
    anyj = pd.merge_asof(
        cnv.sort_values(["user_id", "cnv_ts"]),
        clk_any,
        by="user_id",
        left_on="cnv_ts",
        right_on="clk_ts",
        direction="backward",
    )
    first = clk.drop_duplicates(["user_id", "item_id"], keep="first")
    firstj = pd.merge_asof(
        cnv,
        first.rename(columns={"clk_ts": "first_clk_ts"}),
        by=["user_id", "item_id"],
        left_on="cnv_ts",
        right_on="first_clk_ts",
        direction="backward",
    )
    keys = ["user_id", "item_id", "cnv_ts"]
    out = same[keys + ["clk_ts"]].copy()
    out = out.merge(
        anyj[keys + ["clk_ts"]].rename(columns={"clk_ts": "any_clk_ts"}),
        on=keys,
        how="left",
    )
    out = out.merge(
        firstj[keys + ["first_clk_ts"]],
        on=keys,
        how="left",
    )
    out["wo_prior_clk"] = out["clk_ts"].isna()
    out["dt_item_min"] = (out["cnv_ts"] - out["clk_ts"]) / 60.0
    out["dt_any_min"] = (out["cnv_ts"] - out["any_clk_ts"]) / 60.0
    out["dt_first_min"] = (out["cnv_ts"] - out["first_clk_ts"]) / 60.0
    for b, name in ATTR_BUCKETS:
        out[f"item_within_{name}"] = out["dt_item_min"].le(b) & ~out["wo_prior_clk"]
    return out


def tab_attr_user(attr_ev: pd.DataFrame, users: pd.Index) -> pd.DataFrame:
    if attr_ev.empty:
        return pd.DataFrame({"user_id": users}).assign(
            attr_cnv_n=0,
            attr_wo_clk_cnt=0,
            attr_wo_clk_rate=0.0,
            attr_clk2cnv_min_p50=-1.0,
            attr_anyclk2cnv_min_p50=-1.0,
            attr_anyclk2cnv_min_std=-1.0,
            attr_firstclk2cnv_min_p50=-1.0,
            **{f"attr_clk2cnv_within_{n}_rate": 0.0 for _, n in ATTR_BUCKETS},
        )
    g = attr_ev.groupby("user_id")
    rec = pd.DataFrame(
        {
            "attr_cnv_n": g.size(),
            "attr_wo_clk_cnt": g["wo_prior_clk"].sum(),
            "attr_wo_clk_rate": g["wo_prior_clk"].mean(),
            "attr_clk2cnv_min_p50": g["dt_item_min"].median(),
            "attr_anyclk2cnv_min_p50": g["dt_any_min"].median(),
            "attr_anyclk2cnv_min_std": g["dt_any_min"].std(ddof=0),
            "attr_firstclk2cnv_min_p50": g["dt_first_min"].median(),
        }
    )
    for _, name in ATTR_BUCKETS:
        rec[f"attr_clk2cnv_within_{name}_rate"] = g[f"item_within_{name}"].mean()
    rec = rec.reset_index()
    rec[["attr_clk2cnv_min_p50", "attr_anyclk2cnv_min_p50", "attr_anyclk2cnv_min_std", "attr_firstclk2cnv_min_p50"]] = (
        rec[["attr_clk2cnv_min_p50", "attr_anyclk2cnv_min_p50", "attr_anyclk2cnv_min_std", "attr_firstclk2cnv_min_p50"]]
        .fillna(-1.0)
    )
    return rec


# ----- 块：转移 -----
def tab_trans(ev: pd.DataFrame) -> pd.DataFrame:
    ev = ev.sort_values(["user_id", "ts"]).copy()
    ev["prev"] = ev.groupby("user_id")["act"].shift(1)
    ev["prev2"] = ev.groupby("user_id")["act"].shift(2)
    pairs = ev.dropna(subset=["prev"])
    n1 = pairs.groupby("user_id").size().rename("n_pair")
    def share(a, b):
        hit = pairs[pairs.prev.eq(a) & pairs.act.eq(b)].groupby("user_id").size()
        return (hit / (n1 + 1e-6)).fillna(0.0)

    triples = ev.dropna(subset=["prev2"])
    n2 = triples.groupby("user_id").size().rename("n_triple")
    t2 = triples[triples.prev2.eq(EXP) & triples.prev.eq(CNV) & triples.act.eq(EXP)].groupby("user_id").size()
    out = pd.DataFrame(
        {
            "trans_exp_to_clk": share(EXP, CLK),
            "trans_exp_to_exp": share(EXP, EXP),
            "trans_exp_to_cnv": share(EXP, CNV),
            "trans_cnv_to_exp": share(CNV, EXP),
            "trans_clk_to_cnv": share(CLK, CNV),
            "trans2_exp_cnv_exp": (t2 / (n2 + 1e-6)).fillna(0.0),
        }
    )
    return out.reset_index().rename(columns={"index": "user_id"})


# ----- 块：点击多样性 -----
def tab_diversity(ev: pd.DataFrame) -> pd.DataFrame:
    clk = ev[ev.act == CLK]
    rows = []
    for uid, g in clk.groupby("user_id"):
        vc = g["item_id"].value_counts()
        tot = float(vc.sum())
        h = 0.0
        if tot:
            for c in vc:
                p = c / tot
                h -= p * math.log(p + 1e-6)
        rows.append({"user_id": uid, "item_entropy_clk": h, "n_uniq_clk_item": float(len(vc))})
    div = pd.DataFrame(rows)
    all_u = ev[["user_id"]].drop_duplicates()
    return all_u.merge(div, on="user_id", how="left").fillna({"item_entropy_clk": 0.0, "n_uniq_clk_item": 0.0})


def polish_merge(blocks: List[pd.DataFrame]) -> pd.DataFrame:
    out = blocks[0]
    for b in blocks[1:]:
        out = out.merge(b, on="user_id", how="outer")
    return out


def tab_cross(u: pd.DataFrame) -> pd.DataFrame:
    """原子齐了再乘。"""
    x = u[["user_id"]].copy()

    def col(name, default=0.0):
        return u[name] if name in u.columns else default

    pairs = [
        ("life_ctcvr", "sess_bounce_rate"),
        ("life_cvr", "sess_bounce_rate"),
        ("life_ctr", "active_days"),
        ("life_ctcvr", "active_days"),
        ("life_ctr", "sess_n"),
        ("life_ctcvr", "sess_n"),
        ("pay_cnt", "active_days"),
        ("life_cvr", "7d_ctr"),
        ("life_ctr", "attr_anyclk2cnv_min_p50"),
        ("sess_n", "item_entropy_clk"),
        ("hist_len", "item_entropy_clk"),
        ("life_ctcvr", "pay_cnt"),
    ]
    for a, b in pairs:
        x[f"x_{a}__{b}"] = col(a) * col(b)
    x["log1p_abs_life_ctcvr"] = np.log1p(np.abs(col("life_ctcvr")))
    return x


def n_selected(k: int, n_feat: int, n_row: int) -> int:
    """就要 min(想要的 k, 现有列)。别拿 0.5n 卡选择。"""
    del n_row
    return max(1, min(int(k), int(n_feat)))


def build_user_table(ev: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    funnel = tab_funnel(ev)
    decay = tab_decay(ev)
    sess = tab_session(ev)
    attr_ev = tab_attr_events(ev)
    attr_u = tab_attr_user(attr_ev, funnel["user_id"])
    trans = tab_trans(ev)
    div = tab_diversity(ev)
    atoms = polish_merge([funnel, decay, sess, attr_u, trans, div])
    atoms = atoms.fillna(0)
    user = polish_merge([atoms, tab_cross(atoms)])
    return user, attr_ev


def _demo() -> None:
    t0 = 1_700_000_000
    evs = [
        (11, EXP, t0, 99.0),
        (11, CLK, t0 + 20, 99.0),
        (11, CNV, t0 + 3600, 99.0),
        (22, EXP, t0 + 86400 + 10, 50.0),
        (22, CLK, t0 + 86400 + 30, 50.0),
        (22, CNV, t0 + 86400 + 40, 50.0),  # 10s 就买：同商品 within 5m
    ]
    ev = events_frame(1, evs)
    user, attr = build_user_table(ev)
    print("--- attr 中间表（每笔转化一行）---")
    cols = ["item_id", "wo_prior_clk", "dt_item_min", "dt_any_min", "item_within_5m", "item_within_1h"]
    print(attr[cols].to_string(index=False))
    print("--- user 一块 ---")
    show = [
        "life_cnv_share",
        "sess_n",
        "sess_bounce_rate",
        "attr_wo_clk_rate",
        "attr_clk2cnv_min_p50",
        "attr_clk2cnv_within_5m_rate",
        "trans_exp_to_clk",
        "x_life_ctcvr__sess_bounce_rate",
    ]
    print(user[show].T.to_string(header=False))
    print("n_selected(32, p, n)=", n_selected(32, user.shape[1] - 1, len(user)))


if __name__ == "__main__":
    _demo()
