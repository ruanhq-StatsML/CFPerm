#!/usr/bin/env python3
"""特征计算 prototype。和 VIMP 管道同一口径，可单独跑。

ev 列：user_id, item_id, act∈{0=exp,1=clk,2=cnv}, ts=unix 秒。
每人 [t0, t_end] 时间中位切开：左窗 → X，右窗 → Y。
age 相对**左窗最后一条**，不是全局 t_end。

  python3 scripts/tencent_gr/feat_proto.py
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

EXP, CLK, CNV = 0, 1, 2
HL7 = 7 * 86400.0
SESS_GAP = 30 * 60
MIN_EXP = 3
SPLIT = 0.5
LN2 = math.log(2.0)


def split_lr(ev: pd.DataFrame, frac: float = SPLIT) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mm = ev.groupby("user_id")["ts"].agg(t0="min", t_end="max").reset_index()
    mm["t_cut"] = mm["t0"] + (mm["t_end"] - mm["t0"]) * frac
    e = ev.merge(mm, on="user_id")
    left = e.loc[e.ts <= e.t_cut, ["user_id", "item_id", "act", "ts"]].copy()
    right = e.loc[e.ts > e.t_cut, ["user_id", "item_id", "act", "ts"]].copy()
    return left, right, mm


def x_left(left: pd.DataFrame) -> pd.DataFrame:
    """左窗 X。点击热度 / 点率惯性 / share / 量 / 场 / 只曝不点。"""
    left = left.sort_values(["user_id", "ts"])
    t_left = left.groupby("user_id")["ts"].transform("max")
    age = (t_left - left["ts"]).clip(lower=0)

    is_exp = left.act.eq(EXP)
    is_clk = left.act.eq(CLK)
    is_cnv = left.act.eq(CNV)
    w7 = np.exp(-(LN2 / HL7) * age.to_numpy())  # 2^{-age/7d}

    g = left.groupby("user_id", sort=False)
    n_exp = is_exp.groupby(left.user_id).sum()
    n_clk = is_clk.groupby(left.user_id).sum()
    n_cnv = is_cnv.groupby(left.user_id).sum()
    n_tot = g.size()

    # 点击热度：左窗每个点击 2^{-age/7d} 求和。昨天的点 ≈ 0.91，7 天前 = 0.5。
    dec_clk = pd.Series(w7, index=left.index)[is_clk].groupby(left.loc[is_clk, "user_id"]).sum()
    dec_cnv = pd.Series(w7, index=left.index)[is_cnv].groupby(left.loc[is_cnv, "user_id"]).sum()

    # 场：相邻 >30min 切断。bounce = 一场只有 1 条。
    gap = g["ts"].diff().fillna(0)
    sess = gap.gt(SESS_GAP).groupby(left.user_id).cumsum()
    sn = left.assign(sess=sess).groupby(["user_id", "sess"]).size()
    bounce = sn.le(1).groupby(level=0).mean()
    depth_cnv = (
        left.assign(sess=sess, cnv=is_cnv.astype(int))
        .groupby(["user_id", "sess"])["cnv"]
        .sum()
        .groupby(level=0)
        .mean()
    )

    # 只曝不点：该 item 在左窗 max(act)==0 的占比。
    only_exp = left.groupby(["user_id", "item_id"])["act"].max().eq(0).groupby(level=0).mean()

    x = pd.DataFrame(
        {
            "hist_len": n_tot.astype(float),
            "active_days": left.assign(d=left.ts.floordiv(86400)).groupby("user_id")["d"].nunique(),
            "life_n_exp": n_exp,
            "life_n_clk": n_clk,
            "pay_cnt": n_cnv,  # 量，不是人
            "life_clk_share": n_clk / n_tot.replace(0, np.nan),
            "life_cnv_share": n_cnv / n_tot.replace(0, np.nan),
            "life_ctr": np.where(n_exp >= 1, n_clk / n_exp, 0.0),  # 点率惯性（左窗）
            "dec_hl7d_dec_clk": dec_clk,  # 点击热度
            "dec_hl7d_dec_cnv": dec_cnv,
            "sess_bounce_rate": bounce,
            "sess_depth_cnv_mean": depth_cnv,
            "ui_only_exp_share": only_exp,
        }
    )
    return x.fillna(0.0).reset_index().rename(columns={"index": "user_id"})


def y_right(right: pd.DataFrame) -> pd.DataFrame:
    """右窗 Y。CTR 分母 < MIN_EXP → NaN，不当 0。"""
    y = right.groupby("user_id").agg(
        y_n_exp=("act", lambda a: int((a == EXP).sum())),
        y_n_clk=("act", lambda a: int((a == CLK).sum())),
        y_n_cnv=("act", lambda a: int((a == CNV).sum())),
    )
    y["y_ctr"] = np.where(y.y_n_exp >= MIN_EXP, y.y_n_clk / y.y_n_exp, np.nan)
    y["y_cvr"] = np.where(y.y_n_clk >= 1, y.y_n_cnv / y.y_n_clk, np.nan)
    y["y_any_cnv"] = (y.y_n_cnv > 0).astype(float)
    return y.reset_index()


def demo() -> pd.DataFrame:
    t = 1_700_000_000
    # 一人：左窗热点击 + 右窗还点；一人：左窗只曝、右窗几乎不点
    rows = []
    for i in range(8):
        rows.append((1, 11, EXP, t + i * 3600))
        rows.append((1, 11, CLK, t + i * 3600 + 20))
    rows.append((1, 11, CNV, t + 8 * 3600))
    for i in range(8):
        rows.append((1, 22, EXP, t + 20 * 86400 + i * 1800))
        if i < 3:
            rows.append((1, 22, CLK, t + 20 * 86400 + i * 1800 + 10))
    for i in range(20):
        rows.append((2, 99, EXP, t + i * 86400))
    rows.append((2, 99, CLK, t + 21 * 86400))
    ev = pd.DataFrame(rows, columns=["user_id", "item_id", "act", "ts"])
    left, right, mm = split_lr(ev)
    x = x_left(left)
    y = y_right(right)
    return mm.merge(x, on="user_id").merge(y, on="user_id", how="left")


if __name__ == "__main__":
    out = demo()
    cols = [
        "user_id",
        "dec_hl7d_dec_clk",
        "life_ctr",
        "life_clk_share",
        "pay_cnt",
        "y_ctr",
        "y_n_clk",
    ]
    print(out[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    p = Path("results/tencent_gr_fs150/tables/ev.parquet")
    if p.exists():
        ev = pd.read_parquet(p)
        left, right, mm = split_lr(ev)
        x = x_left(left)
        print("real left X head")
        print(x.head(3).to_string(index=False, float_format=lambda v: f"{v:.4f}"))
