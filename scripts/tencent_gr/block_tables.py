#!/usr/bin/env python3
"""Block tables → polish → merge. Attribution is merge_asof, not a dict walk.

Grain
-----
  ev     user × event
  attr   user × conversion   backward asof：买之前怎么点过来
  post   user × conversion   forward asof：买完会不会接着点
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
# 买后点击窗：5m 还在当场、1h 余热、1d 回访、7d 中期
POST_WINS = {"5m": 5 * 60, "1h": 3600, "1d": 86400, "7d": 7 * 86400}
# 预测「这一单之后会不会点」：只许 cnv_ts 已知的列（Y / lift / n_after 不准进）
POST_PREDICT_X = [
    "wo_prior_clk",
    "dt_item_min",
    "dt_any_min",
    "item_within_5m",
    "item_within_1h",
    "n_clk_before_1h",
    "n_clk_before_1d",
    "n_clk_before_7d",
    "n_clk_same_before",
    "sess_pos",
    "sess_clk_before",
    "price",
    "log1p_price",
    "n_prior_cnv",
    "lag_post_clk_1d_rate",
]
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


def _asof_cum(left: pd.DataFrame, clk: pd.DataFrame, by: List[str], left_on: str, cum_col: str) -> np.ndarray:
    """每行 left：by 组里 ts<=left_on 的点击累计。空 → 0。"""
    n = len(left)
    if n == 0 or clk.empty or cum_col not in clk.columns:
        return np.zeros(n, dtype=float)
    right = (
        clk[list(by) + ["ts", cum_col]]
        .rename(columns={"ts": "_rts", cum_col: "_cum"})
        .sort_values(list(by) + ["_rts"])
    )
    tmp = left[list(by) + [left_on]].copy()
    tmp["_i"] = np.arange(n)
    tmp = tmp.sort_values(list(by) + [left_on])
    j = pd.merge_asof(
        tmp,
        right,
        by=list(by),
        left_on=left_on,
        right_on="_rts",
        direction="backward",
        allow_exact_matches=True,
    )
    j = j.sort_values("_i")
    return j["_cum"].fillna(0.0).to_numpy(dtype=float)


def _asof_next(left: pd.DataFrame, clk: pd.DataFrame, by: List[str], left_on: str) -> pd.DataFrame:
    """每行 left：by 组里严格晚于 left_on 的下一次点击。"""
    want_item = "item_id" in clk.columns and "item_id" not in by
    empty = pd.DataFrame({"_next_ts": np.full(len(left), np.nan)})
    if want_item:
        empty["_next_item"] = np.nan
    if left.empty or clk.empty:
        return empty
    keep = list(by) + ["ts"]
    rename = {"ts": "_next_ts"}
    if want_item:
        keep = list(by) + ["item_id", "ts"]
        rename["item_id"] = "_next_item"
    right = clk[keep].rename(columns=rename).sort_values(list(by) + ["_next_ts"])
    tmp = left[list(by) + [left_on]].copy()
    tmp["_i"] = np.arange(len(left))
    tmp = tmp.sort_values(list(by) + [left_on])
    j = pd.merge_asof(
        tmp,
        right,
        by=list(by),
        left_on=left_on,
        right_on="_next_ts",
        direction="forward",
        allow_exact_matches=False,
    )
    j = j.sort_values("_i")
    out = pd.DataFrame({"_next_ts": j["_next_ts"].to_numpy()})
    if want_item:
        out["_next_item"] = j["_next_item"].to_numpy()
    return out


def _y_hit(dt_sec: pd.Series, t_end: pd.Series, cnv_ts: pd.Series, win: int) -> pd.Series:
    """删失：窗内已点到 → 1；跟满窗没点 → 0；跟不满且没点 → NaN。"""
    hit = dt_sec.notna() & dt_sec.le(win)
    follow = t_end - cnv_ts
    y = np.where(hit, 1.0, np.where(follow.ge(win), 0.0, np.nan))
    return pd.Series(y, index=dt_sec.index)


# ----- 块：买后点击（中间表 = 每笔转化一行，forward asof）-----
def tab_post_cnv_events(ev: pd.DataFrame, attr_ev: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """买完会不会接着点。Markov 邻接转移不够：这里是窗内次数 + 下一次点击时延。

    Y（预测目标，不能进 X）：y_post_clk_* / y_post_same_* / n_clk_after_* / lift_*
    X（cnv_ts 已知）：买前路径、买前窗点击量、当场深度、价格、此前各单的滞后买后点击率。
    """
    cnv = (
        ev.loc[ev.act == CNV, ["user_id", "item_id", "ts", "price"]]
        .rename(columns={"ts": "cnv_ts"})
        .reset_index(drop=True)
    )
    t_end = ev.groupby("user_id")["ts"].max().rename("t_end")
    if cnv.empty:
        return cnv.assign(t_end=np.nan, wo_prior_clk=True)

    clk = ev.loc[ev.act == CLK, ["user_id", "item_id", "ts"]].copy()
    if clk.empty:
        clk = clk.assign(cum_any=pd.Series(dtype=float), cum_item=pd.Series(dtype=float))
    else:
        clk = clk.sort_values(["user_id", "ts"])
        clk["cum_any"] = clk.groupby("user_id").cumcount() + 1
        clk = clk.sort_values(["user_id", "item_id", "ts"])
        clk["cum_item"] = clk.groupby(["user_id", "item_id"]).cumcount() + 1

    nxt_any = _asof_next(cnv, clk, ["user_id"], "cnv_ts")
    nxt_item = _asof_next(cnv, clk, ["user_id", "item_id"], "cnv_ts")
    out = cnv.copy()
    out["t_end"] = out["user_id"].map(t_end)
    out["next_clk_ts"] = nxt_any["_next_ts"].to_numpy()
    if "_next_item" in nxt_any.columns:
        out["next_clk_item"] = nxt_any["_next_item"].to_numpy()
    out["next_same_ts"] = nxt_item["_next_ts"].to_numpy()
    out["dt_next_clk_sec"] = out["next_clk_ts"] - out["cnv_ts"]
    out["dt_next_same_sec"] = out["next_same_ts"] - out["cnv_ts"]
    out["dt_next_clk_min"] = out["dt_next_clk_sec"] / 60.0
    out["dt_next_same_min"] = out["dt_next_same_sec"] / 60.0
    out["next_clk_same_sess"] = np.where(
        out["next_clk_ts"].notna(),
        out["dt_next_clk_sec"].le(SESS_GAP).astype(float),
        np.nan,
    )

    # 买前/买后同长窗点击量：cum(t+W)-cum(t) vs cum(t)-cum(t-W)
    out["n_clk_le"] = _asof_cum(out, clk, ["user_id"], "cnv_ts", "cum_any")
    out["n_clk_same_le"] = _asof_cum(out, clk, ["user_id", "item_id"], "cnv_ts", "cum_item")
    out["n_clk_same_before"] = out["n_clk_same_le"]
    for name, sec in POST_WINS.items():
        out[f"ts_p_{name}"] = out["cnv_ts"] + sec
        out[f"ts_m_{name}"] = out["cnv_ts"] - sec
        c_after = _asof_cum(out, clk, ["user_id"], f"ts_p_{name}", "cum_any")
        c_before = _asof_cum(out, clk, ["user_id"], f"ts_m_{name}", "cum_any")
        c_after_item = _asof_cum(out, clk, ["user_id", "item_id"], f"ts_p_{name}", "cum_item")
        c_before_item = _asof_cum(out, clk, ["user_id", "item_id"], f"ts_m_{name}", "cum_item")
        out[f"n_clk_after_{name}"] = c_after - out["n_clk_le"]
        out[f"n_clk_before_{name}"] = out["n_clk_le"] - c_before
        out[f"n_same_after_{name}"] = c_after_item - out["n_clk_same_le"]
        out[f"n_same_before_{name}"] = out["n_clk_same_le"] - c_before_item
        out[f"lift_{name}"] = out[f"n_clk_after_{name}"] / (out[f"n_clk_before_{name}"] + 1.0)
        out[f"delta_{name}"] = out[f"n_clk_after_{name}"] - out[f"n_clk_before_{name}"]
        out[f"lift_same_{name}"] = out[f"n_same_after_{name}"] / (out[f"n_same_before_{name}"] + 1.0)
        out[f"y_post_clk_{name}"] = _y_hit(out["dt_next_clk_sec"], out["t_end"], out["cnv_ts"], sec)
        out[f"y_post_same_{name}"] = _y_hit(out["dt_next_same_sec"], out["t_end"], out["cnv_ts"], sec)
        obs = (out["t_end"] - out["cnv_ts"]) >= sec
        for col in (
            f"n_clk_after_{name}",
            f"n_same_after_{name}",
            f"lift_{name}",
            f"delta_{name}",
            f"lift_same_{name}",
        ):
            out.loc[~obs, col] = np.nan

    out["log1p_price"] = np.log1p(out["price"].fillna(0.0).clip(lower=0.0))

    # 当场（买这一刻之前的场深）——热场续点 vs 隔天回访
    ev_s = ev.sort_values(["user_id", "ts"]).copy()
    gap = ev_s.groupby("user_id")["ts"].diff().fillna(0)
    ev_s["sess"] = gap.gt(SESS_GAP).groupby(ev_s["user_id"]).cumsum()
    ev_s["sess_pos"] = ev_s.groupby(["user_id", "sess"]).cumcount() + 1
    ev_s["_clk"] = ev_s.act.eq(CLK).astype(int)
    ev_s["sess_clk_before"] = ev_s.groupby(["user_id", "sess"])["_clk"].cumsum()
    # 当前行若是点击，cumsum 已含自己；转化行不是点击，就是场内已有点击数
    hit = ev_s.loc[ev_s.act == CNV, ["user_id", "item_id", "ts", "sess_pos", "sess_clk_before"]]
    hit = hit.rename(columns={"ts": "cnv_ts"}).drop_duplicates(["user_id", "item_id", "cnv_ts"])
    out = out.merge(hit, on=["user_id", "item_id", "cnv_ts"], how="left")

    if attr_ev is not None and len(attr_ev):
        akeys = ["user_id", "item_id", "cnv_ts"]
        acols = [c for c in attr_ev.columns if c in akeys or c in (
            "wo_prior_clk", "dt_item_min", "dt_any_min", "dt_first_min",
            *(f"item_within_{n}" for _, n in ATTR_BUCKETS),
        )]
        out = out.merge(attr_ev[acols], on=akeys, how="left")
    if "wo_prior_clk" not in out.columns:
        out["wo_prior_clk"] = False

    out = out.sort_values(["user_id", "cnv_ts"]).reset_index(drop=True)
    out["n_prior_cnv"] = out.groupby("user_id").cumcount()
    shifted = out.groupby("user_id")["y_post_clk_1d"].shift(1)
    ok = shifted.notna()
    csum = shifted.fillna(0.0).groupby(out["user_id"]).cumsum()
    ccnt = ok.groupby(out["user_id"]).cumsum()
    out["lag_post_clk_1d_rate"] = np.where(ccnt > 0, csum / ccnt, np.nan)

    drop_tmp = [c for c in out.columns if c.startswith("ts_p_") or c.startswith("ts_m_")]
    return out.drop(columns=drop_tmp)


def tab_post_cnv_user(post_ev: pd.DataFrame, users) -> pd.DataFrame:
    """用户一行：历史上「买完还点」的倾向。只平均未删失的转化。"""
    users = pd.Index(users)
    base = pd.DataFrame({"user_id": users})
    zeros = {
        "post_n_obs_1d": 0.0,
        "post_clk_5m_rate": 0.0,
        "post_clk_1h_rate": 0.0,
        "post_clk_1d_rate": 0.0,
        "post_clk_7d_rate": 0.0,
        "post_same_1h_rate": 0.0,
        "post_same_1d_rate": 0.0,
        "post_clk_same_sess_rate": 0.0,
        "post_clk_dt_p50": -1.0,
        "post_lift_1d_p50": 0.0,
        "post_delta_1d_p50": 0.0,
        "post_lift_same_1d_p50": 0.0,
        "post_n_clk_after_1d_mean": 0.0,
        "post_n_clk_before_1d_mean": 0.0,
    }
    if post_ev is None or post_ev.empty:
        return base.assign(**zeros)

    def mean_obs(col: str) -> pd.Series:
        s = post_ev.groupby("user_id")[col].mean()
        return s

    rec = pd.DataFrame(
        {
            "post_n_obs_1d": post_ev.groupby("user_id")["y_post_clk_1d"].count(),
            "post_clk_5m_rate": mean_obs("y_post_clk_5m"),
            "post_clk_1h_rate": mean_obs("y_post_clk_1h"),
            "post_clk_1d_rate": mean_obs("y_post_clk_1d"),
            "post_clk_7d_rate": mean_obs("y_post_clk_7d"),
            "post_same_1h_rate": mean_obs("y_post_same_1h"),
            "post_same_1d_rate": mean_obs("y_post_same_1d"),
            "post_clk_same_sess_rate": mean_obs("next_clk_same_sess"),
            "post_clk_dt_p50": post_ev.groupby("user_id")["dt_next_clk_min"].median(),
            "post_lift_1d_p50": post_ev.groupby("user_id")["lift_1d"].median(),
            "post_delta_1d_p50": post_ev.groupby("user_id")["delta_1d"].median(),
            "post_lift_same_1d_p50": post_ev.groupby("user_id")["lift_same_1d"].median(),
            "post_n_clk_after_1d_mean": mean_obs("n_clk_after_1d"),
            "post_n_clk_before_1d_mean": mean_obs("n_clk_before_1d"),
        }
    ).reset_index()
    rec["post_clk_dt_p50"] = rec["post_clk_dt_p50"].fillna(-1.0)
    rec["post_clk_same_sess_rate"] = rec["post_clk_same_sess_rate"].fillna(0.0)
    out = base.merge(rec, on="user_id", how="left")
    for k, v in zeros.items():
        out[k] = out[k].fillna(v)
    return out


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
        ("post_clk_1d_rate", "life_ctcvr"),
        ("post_lift_1d_p50", "life_cvr"),
        ("post_clk_1d_rate", "attr_wo_clk_rate"),
        ("post_same_1d_rate", "life_cvr"),
    ]
    for a, b in pairs:
        x[f"x_{a}__{b}"] = col(a) * col(b)
    x["log1p_abs_life_ctcvr"] = np.log1p(np.abs(col("life_ctcvr")))
    return x


def n_selected(k: int, n_feat: int, n_row: int) -> int:
    """就要 min(想要的 k, 现有列)。别拿 0.5n 卡选择。"""
    del n_row
    return max(1, min(int(k), int(n_feat)))


def build_user_table(ev: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    funnel = tab_funnel(ev)
    decay = tab_decay(ev)
    sess = tab_session(ev)
    attr_ev = tab_attr_events(ev)
    attr_u = tab_attr_user(attr_ev, funnel["user_id"])
    post_ev = tab_post_cnv_events(ev, attr_ev)
    post_u = tab_post_cnv_user(post_ev, funnel["user_id"])
    trans = tab_trans(ev)
    div = tab_diversity(ev)
    atoms = polish_merge([funnel, decay, sess, attr_u, post_u, trans, div])
    atoms = atoms.fillna(0)
    user = polish_merge([atoms, tab_cross(atoms)])
    return user, attr_ev, post_ev


def _demo() -> None:
    t0 = 1_700_000_000
    evs = [
        (11, EXP, t0, 99.0),
        (11, CLK, t0 + 20, 99.0),
        (11, CNV, t0 + 3600, 99.0),
        (11, CLK, t0 + 3600 + 600, 99.0),  # 买后 10min 还点同品 → 当场续点
        (22, EXP, t0 + 86400 + 10, 50.0),
        (22, CLK, t0 + 86400 + 30, 50.0),
        (22, CNV, t0 + 86400 + 40, 50.0),  # 10s 就买；后面不再点
        (99, EXP, t0 + 20 * 86400, None),  # 垫 t_end，两单的 7d 窗都跟满
    ]
    ev = events_frame(1, evs)
    user, attr, post = build_user_table(ev)
    print("--- attr 中间表（买之前）---")
    cols = ["item_id", "wo_prior_clk", "dt_item_min", "dt_any_min", "item_within_5m", "item_within_1h"]
    print(attr[cols].to_string(index=False))
    print("--- post 中间表（买之后）---")
    pcols = [
        "item_id",
        "y_post_clk_1h",
        "y_post_clk_1d",
        "y_post_same_1d",
        "n_clk_before_1d",
        "n_clk_after_1d",
        "lift_1d",
        "next_clk_same_sess",
        "lag_post_clk_1d_rate",
    ]
    print(post[pcols].to_string(index=False))
    print("--- user 一块 ---")
    show = [
        "life_cnv_share",
        "sess_n",
        "sess_bounce_rate",
        "attr_wo_clk_rate",
        "attr_clk2cnv_min_p50",
        "attr_clk2cnv_within_5m_rate",
        "post_clk_1d_rate",
        "post_same_1d_rate",
        "post_lift_1d_p50",
        "trans_exp_to_clk",
        "x_life_ctcvr__sess_bounce_rate",
        "x_post_clk_1d_rate__life_ctcvr",
    ]
    print(user[show].T.to_string(header=False))
    print("n_selected(32, p, n)=", n_selected(32, user.shape[1] - 1, len(user)))
    assert float(post.iloc[0]["y_post_clk_1h"]) == 1.0
    assert float(post.iloc[1]["y_post_clk_1d"]) == 0.0
    assert float(post.iloc[0]["next_clk_same_sess"]) == 1.0


if __name__ == "__main__":
    _demo()
