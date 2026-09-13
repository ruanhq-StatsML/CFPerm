"""Flatten nested seq → event table, then user-level behavior counts.

This is Deep Feature Synthesis by hand (entity: events → users). Counts and
rates are the right objects to StandardScaler for a tabular / DR path. They
are *not* a unique CS/CD decomposition and *not* a CATE.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

ACTION_EXPOSURE = 0
ACTION_CLICK = 1
ACTION_CONVERSION = 2

# Numeric columns that are counts / rates / times — scale these.
SCALE_COLS: tuple[str, ...] = (
    "n_events",
    "n_exposure",
    "n_click",
    "n_conversion",
    "n_unique_items",
    "n_unique_click_items",
    "n_unique_conv_items",
    "repeat_rate",
    "click_rate",
    "conversion_rate",
    "engage_rate",
    "span_days",
    "mean_gap_hours",
    "median_gap_hours",
    "recency_to_end_days",
    "pos_last_click",
    "pos_last_conv",
    "n_item_format",
    "n_clicks_last20",
    "n_conv_last20",
    "engage_last20",
)


def explode_seq(seq_df: pd.DataFrame) -> pd.DataFrame:
    """Nested `seq` list[dict] → long event table."""
    rows: list[dict] = []
    for rec in seq_df.itertuples(index=False):
        uid = int(rec.user_id)
        seq = rec.seq if rec.seq is not None else []
        n = len(seq)
        for pos, ev in enumerate(seq):
            if not isinstance(ev, dict):
                continue
            ts = ev.get("timestamp")
            ts = int(ts) if ts is not None else 0
            rows.append(
                {
                    "user_id": uid,
                    "pos": int(pos),
                    "pos_from_end": int(n - 1 - pos),
                    "item_id": int(ev.get("item_id", 0) or 0),
                    "action_type": int(ev.get("action_type", -1)),
                    "timestamp": ts,
                }
            )
    events = pd.DataFrame(rows)
    if events.empty:
        return events
    events["is_click"] = (events["action_type"] == ACTION_CLICK).astype(np.int8)
    events["is_conversion"] = (events["action_type"] == ACTION_CONVERSION).astype(np.int8)
    events["is_engage"] = ((events["is_click"] == 1) | (events["is_conversion"] == 1)).astype(np.int8)
    events["is_exposure"] = (events["action_type"] == ACTION_EXPOSURE).astype(np.int8)
    return events


def _gap_hours(ts: np.ndarray) -> tuple[float, float]:
    if ts.size < 2:
        return 0.0, 0.0
    order = np.sort(ts.astype(np.float64))
    gaps = np.diff(order) / 3600.0
    return float(np.mean(gaps)), float(np.median(gaps))


def synthesize_user_behavior(
    events: pd.DataFrame,
    *,
    user_feat: Optional[pd.DataFrame] = None,
    item_feat: Optional[pd.DataFrame] = None,
    global_tmax: Optional[int] = None,
) -> pd.DataFrame:
    """DFS-style primitives: COUNT / NUM_UNIQUE / MEAN / LAST / TIME_SINCE on events.

    Second-level: same primitives on click-only and conversion-only slices,
    plus last-20 window (recency-weighted counts).
    """
    if events.empty:
        return pd.DataFrame()

    tmax = int(global_tmax) if global_tmax is not None else int(events["timestamp"].max())
    item_format = None
    if item_feat is not None and "item_id" in item_feat.columns and "100" in item_feat.columns:
        item_format = item_feat[["item_id", "100"]].drop_duplicates("item_id")
        item_format = item_format.rename(columns={"100": "item_format"})

    recs: list[dict] = []
    for uid, g in events.groupby("user_id", sort=False):
        g = g.sort_values("pos")
        ts = g["timestamp"].to_numpy()
        acts = g["action_type"].to_numpy()
        items = g["item_id"].to_numpy()
        n = int(len(g))
        n_exp = int((acts == ACTION_EXPOSURE).sum())
        n_clk = int((acts == ACTION_CLICK).sum())
        n_cnv = int((acts == ACTION_CONVERSION).sum())
        n_uni = int(np.unique(items).size)
        click_items = items[acts == ACTION_CLICK]
        conv_items = items[acts == ACTION_CONVERSION]
        mean_gap, med_gap = _gap_hours(ts)
        tmin, tmx = int(ts.min()), int(ts.max())
        last20 = g[g["pos_from_end"] < 20]
        last_clk = np.where(acts == ACTION_CLICK)[0]
        last_cnv = np.where(acts == ACTION_CONVERSION)[0]
        n_fmt = 0
        if item_format is not None:
            merged = g.merge(item_format, on="item_id", how="left")
            n_fmt = int(merged["item_format"].nunique(dropna=True))
        recs.append(
            {
                "user_id": int(uid),
                "n_events": n,
                "n_exposure": n_exp,
                "n_click": n_clk,
                "n_conversion": n_cnv,
                "n_unique_items": n_uni,
                "n_unique_click_items": int(np.unique(click_items).size) if click_items.size else 0,
                "n_unique_conv_items": int(np.unique(conv_items).size) if conv_items.size else 0,
                "repeat_rate": float(1.0 - n_uni / n) if n else 0.0,
                "click_rate": float(n_clk / n) if n else 0.0,
                "conversion_rate": float(n_cnv / n) if n else 0.0,
                "engage_rate": float((n_clk + n_cnv) / n) if n else 0.0,
                "span_days": float((tmx - tmin) / 86400.0),
                "mean_gap_hours": mean_gap,
                "median_gap_hours": med_gap,
                "recency_to_end_days": float((tmax - tmx) / 86400.0),
                "tmin": tmin,
                "tmax": tmx,
                "last_action": int(acts[-1]),
                "pos_last_click": int(last_clk[-1]) if last_clk.size else -1,
                "pos_last_conv": int(last_cnv[-1]) if last_cnv.size else -1,
                "any_click": int(n_clk > 0),
                "any_conversion": int(n_cnv > 0),
                "n_item_format": n_fmt,
                "n_clicks_last20": int((last20["action_type"] == ACTION_CLICK).sum()),
                "n_conv_last20": int((last20["action_type"] == ACTION_CONVERSION).sum()),
                "engage_last20": int((last20["is_engage"] == 1).sum()) if "is_engage" in last20 else int(
                    ((last20["action_type"] == ACTION_CLICK) | (last20["action_type"] == ACTION_CONVERSION)).sum()
                ),
            }
        )
    users = pd.DataFrame(recs)
    if user_feat is not None and "user_id" in user_feat.columns:
        uids = set(user_feat["user_id"].astype(int).tolist())
        users["has_user_feat"] = users["user_id"].isin(uids).astype(np.int8)
        for col in ("103", "104", "105", "109"):
            if col in user_feat.columns:
                tmp = user_feat[["user_id", col]].copy()
                tmp["user_id"] = tmp["user_id"].astype(int)
                tmp = tmp.rename(columns={col: f"user_{col}"})
                users = users.merge(tmp, on="user_id", how="left")
    else:
        users["has_user_feat"] = np.int8(0)
    mid = float(users["tmax"].median())
    users["time_late"] = (users["tmax"] >= mid).astype(np.int8)
    users["time_split_median"] = mid
    return users


def standard_scale_behavior(
    users: pd.DataFrame,
    cols: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, StandardScaler, list[str]]:
    """z-score the count/rate columns. Fit on the table you pass (document leak if you later split)."""
    use = [c for c in (cols or SCALE_COLS) if c in users.columns]
    scaler = StandardScaler()
    mat = users[use].to_numpy(dtype=np.float64)
    mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
    z = scaler.fit_transform(mat)
    scaled = users.copy()
    for i, c in enumerate(use):
        scaled[f"z_{c}"] = z[:, i]
    return scaled, scaler, use


def z_matrix(scaled: pd.DataFrame, cols: Iterable[str]) -> np.ndarray:
    names = [f"z_{c}" for c in cols if f"z_{c}" in scaled.columns]
    return scaled[names].to_numpy(dtype=np.float64)
