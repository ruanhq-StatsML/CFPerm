"""User × Ctx 交叉。SKU 交叉跟 GATE。

三塔 fusion 的显式版：人怎么来 × 这单怎么来。
不要 heat × empty：empty_any=1 ⇒ prefix 无点击 ⇒ n_clk_7d=0，乘积恒 0。

  from cross_feats import add_cross
  x = add_cross(prep(post, cols), sku_on=False)
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# (user, ctx) 名字只是注释。乘的是已经填好的列。
PAIRS_UC = [
    ("lag_post_clk_1d_rate", "empty_any"),       # 习惯会点 × 这单空路径
    ("lag_empty_any", "empty_any"),              # 人经常空 × 这单空
    ("lag_post_clk_1d_rate", "dt_any_log"),      # 习惯会点 × 路径远近
    ("n_clk_7d_log", "dt_any_log"),              # 7d 热 × 路径远近
    ("lag_post_clk_1d_rate", "sess_clk_before"), # 习惯会点 × 当场在逛
    ("n_clk_7d_log", "sess_clk_before"),         # 7d 热 × 当场
]

# GATE=0 不建。empty × wo_prior_clk 几乎共线，不开。
PAIRS_SKU = [
    ("lag_post_clk_1d_rate", "wo_prior_clk"),
    ("n_clk_7d_log", "n_clk_same_log"),
]


def _num(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(default)
    return pd.Series(default, index=df.index, dtype=float)


def bases(df: pd.DataFrame) -> pd.DataFrame:
    """交叉用的派生列。dt miss → dt_any_log=0（没有路径可乘）。"""
    b = pd.DataFrame(index=df.index)
    b["empty_any"] = _num(df, "empty_any")
    b["lag_post_clk_1d_rate"] = _num(df, "lag_post_clk_1d_rate")
    b["lag_empty_any"] = _num(df, "lag_empty_any")
    b["sess_clk_before"] = _num(df, "sess_clk_before")
    b["n_clk_7d_log"] = np.log1p(_num(df, "n_clk_before_7d").clip(lower=0))
    b["n_clk_same_log"] = np.log1p(_num(df, "n_clk_same_before").clip(lower=0))
    b["wo_prior_clk"] = _num(df, "wo_prior_clk")
    if "dt_any_min" in df.columns:
        dt_raw = pd.to_numeric(df["dt_any_min"], errors="coerce")
    else:
        dt_raw = pd.Series(np.nan, index=df.index)
    if "dt_any_min_miss" in df.columns:
        miss = _num(df, "dt_any_min_miss") > 0
    else:
        miss = dt_raw.isna()
    miss = miss.astype(bool) | (dt_raw < 0).fillna(False)
    dt = dt_raw.clip(lower=0).fillna(0.0)
    b["dt_any_log"] = np.where(miss.to_numpy(), 0.0, np.log1p(dt.to_numpy()))
    return b


def mul_pairs(b: pd.DataFrame, pairs: list[tuple[str, str]], prefix: str) -> pd.DataFrame:
    out = pd.DataFrame(index=b.index)
    for a, c in pairs:
        out[f"{prefix}{a}__x__{c}"] = b[a].astype(float) * b[c].astype(float)
    return out


def add_cross(x: pd.DataFrame, sku_on: bool = False) -> pd.DataFrame:
    b = bases(x)
    parts = [x, mul_pairs(b, PAIRS_UC, "uc_")]
    if sku_on:
        parts.append(mul_pairs(b, PAIRS_SKU, "us_"))
    return pd.concat(parts, axis=1)
