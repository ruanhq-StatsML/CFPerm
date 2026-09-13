"""Flatten nested seq → event table, then user-level sequence summaries.

Pandas Deep Feature Synthesis (events → users). Primitives map to Featuretools:

    COUNT, NUM_UNIQUE, SUM, MEAN, MAX, PERCENT_TRUE, TIME_SINCE,
    TIME_SINCE_PREVIOUS, LAST, plus last-k / exponential-decay windows.

Counts and rates are the objects to StandardScaler for a tabular / DR path.
They rank *summaries* of past sequence behavior. They are not a unique
CS/CD decomposition and not a CATE.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
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
    "click_rate_last20",
    "decay_engage",
    "n_exp_to_click",
    "mean_item_log_freq",
)

# Features that are *not* a direct rewrite of any_click / n_click.
STRUCTURAL_COLS: tuple[str, ...] = (
    "n_events",
    "n_unique_items",
    "repeat_rate",
    "span_days",
    "mean_gap_hours",
    "median_gap_hours",
    "recency_to_end_days",
    "n_item_format",
    "mean_item_log_freq",
    "n_exp_to_click",
)


def explode_seq(seq_df: pd.DataFrame) -> pd.DataFrame:
    """Nested `seq` list[dict] → long event table (one row per action)."""
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
    events["event_id"] = np.arange(len(events), dtype=np.int64)
    events["is_click"] = (events["action_type"] == ACTION_CLICK).astype(np.int8)
    events["is_conversion"] = (events["action_type"] == ACTION_CONVERSION).astype(np.int8)
    events["is_engage"] = ((events["is_click"] == 1) | (events["is_conversion"] == 1)).astype(np.int8)
    events["is_exposure"] = (events["action_type"] == ACTION_EXPOSURE).astype(np.int8)
    return events


def _group_nunique(df: pd.DataFrame, mask: pd.Series, key: str, col: str) -> pd.Series:
    sub = df.loc[mask, ["user_id", col]]
    if sub.empty:
        return pd.Series(dtype=np.int64, name=key)
    return sub.groupby("user_id")[col].nunique().rename(key)


def synthesize_user_behavior(
    events: pd.DataFrame,
    *,
    user_feat: Optional[pd.DataFrame] = None,
    item_feat: Optional[pd.DataFrame] = None,
    global_tmax: Optional[int] = None,
    decay_halflife: float = 10.0,
) -> pd.DataFrame:
    """DFS-style primitives: COUNT / NUM_UNIQUE / MEAN / LAST / TIME_SINCE on events.

    Second-level: click-only and conversion-only slices, last-20 window,
    exponential recency weights, exposure→click transitions, item popularity.
    """
    if events.empty:
        return pd.DataFrame()

    ev = events.copy()
    tmax = int(global_tmax) if global_tmax is not None else int(ev["timestamp"].max())

    item_freq = ev.groupby("item_id").size().rename("item_freq")
    ev = ev.join(item_freq, on="item_id")
    ev["item_log_freq"] = np.log1p(ev["item_freq"].to_numpy(dtype=np.float64))
    ev["decay_w"] = np.exp(-ev["pos_from_end"].to_numpy(dtype=np.float64) / float(decay_halflife))
    ev["decay_engage"] = ev["is_engage"].to_numpy(dtype=np.float64) * ev["decay_w"]

    ev = ev.sort_values(["user_id", "pos"], kind="mergesort")
    prev = ev.groupby("user_id", sort=False)["action_type"].shift(1)
    ev["exp_to_click"] = ((prev == ACTION_EXPOSURE) & (ev["action_type"] == ACTION_CLICK)).astype(np.int8)
    gap_h = ev.groupby("user_id", sort=False)["timestamp"].diff() / 3600.0
    ev["gap_hours"] = gap_h

    if item_feat is not None and "item_id" in item_feat.columns and "100" in item_feat.columns:
        fmt = item_feat[["item_id", "100"]].drop_duplicates("item_id").rename(columns={"100": "item_format"})
        ev = ev.merge(fmt, on="item_id", how="left")
    else:
        ev["item_format"] = np.nan

    g = ev.groupby("user_id", sort=False)
    users = g.agg(
        n_events=("event_id", "size") if "event_id" in ev.columns else ("item_id", "size"),
        n_exposure=("is_exposure", "sum"),
        n_click=("is_click", "sum"),
        n_conversion=("is_conversion", "sum"),
        n_unique_items=("item_id", "nunique"),
        tmin=("timestamp", "min"),
        tmax=("timestamp", "max"),
        last_action=("action_type", "last"),
        n_item_format=("item_format", "nunique"),
        n_exp_to_click=("exp_to_click", "sum"),
        mean_item_log_freq=("item_log_freq", "mean"),
        decay_engage=("decay_engage", "sum"),
        mean_gap_hours=("gap_hours", "mean"),
        median_gap_hours=("gap_hours", "median"),
    ).reset_index()

    clk_u = _group_nunique(ev, ev["is_click"] == 1, "n_unique_click_items", "item_id")
    cnv_u = _group_nunique(ev, ev["is_conversion"] == 1, "n_unique_conv_items", "item_id")
    users = users.merge(clk_u, left_on="user_id", right_index=True, how="left")
    users = users.merge(cnv_u, left_on="user_id", right_index=True, how="left")
    users["n_unique_click_items"] = users["n_unique_click_items"].fillna(0).astype(np.int64)
    users["n_unique_conv_items"] = users["n_unique_conv_items"].fillna(0).astype(np.int64)

    pos_clk = ev.loc[ev["is_click"] == 1].groupby("user_id")["pos"].max().rename("pos_last_click")
    pos_cnv = ev.loc[ev["is_conversion"] == 1].groupby("user_id")["pos"].max().rename("pos_last_conv")
    users = users.merge(pos_clk, left_on="user_id", right_index=True, how="left")
    users = users.merge(pos_cnv, left_on="user_id", right_index=True, how="left")
    users["pos_last_click"] = users["pos_last_click"].fillna(-1).astype(np.int64)
    users["pos_last_conv"] = users["pos_last_conv"].fillna(-1).astype(np.int64)

    last20 = ev.loc[ev["pos_from_end"] < 20]
    l20 = last20.groupby("user_id").agg(
        n_clicks_last20=("is_click", "sum"),
        n_conv_last20=("is_conversion", "sum"),
        engage_last20=("is_engage", "sum"),
        n_last20=("item_id", "size"),
    )
    users = users.merge(l20, left_on="user_id", right_index=True, how="left")
    for c in ("n_clicks_last20", "n_conv_last20", "engage_last20", "n_last20"):
        users[c] = users[c].fillna(0)
    n = users["n_events"].to_numpy(dtype=np.float64)
    n = np.maximum(n, 1.0)
    users["repeat_rate"] = 1.0 - users["n_unique_items"].to_numpy(dtype=np.float64) / n
    users["click_rate"] = users["n_click"].to_numpy(dtype=np.float64) / n
    users["conversion_rate"] = users["n_conversion"].to_numpy(dtype=np.float64) / n
    users["engage_rate"] = (users["n_click"] + users["n_conversion"]).to_numpy(dtype=np.float64) / n
    users["span_days"] = (users["tmax"] - users["tmin"]).to_numpy(dtype=np.float64) / 86400.0
    users["recency_to_end_days"] = (tmax - users["tmax"]).to_numpy(dtype=np.float64) / 86400.0
    users["click_rate_last20"] = users["n_clicks_last20"].to_numpy(dtype=np.float64) / np.maximum(
        users["n_last20"].to_numpy(dtype=np.float64), 1.0
    )
    users["any_click"] = (users["n_click"] > 0).astype(np.int8)
    users["any_conversion"] = (users["n_conversion"] > 0).astype(np.int8)
    users["mean_gap_hours"] = users["mean_gap_hours"].fillna(0.0)
    users["median_gap_hours"] = users["median_gap_hours"].fillna(0.0)
    users["mean_item_log_freq"] = users["mean_item_log_freq"].fillna(0.0)
    users["n_item_format"] = users["n_item_format"].fillna(0).astype(np.int64)

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


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return 0.0
    va, vb = float(a.var(ddof=1)), float(b.var(ddof=1))
    pooled = np.sqrt(((na - 1) * va + (nb - 1) * vb) / max(na + nb - 2, 1))
    if pooled < 1e-12:
        return 0.0
    return float((a.mean() - b.mean()) / pooled)


def rank_against_binary(scaled: pd.DataFrame, cols: Sequence[str], y: np.ndarray) -> pd.DataFrame:
    """Univariate localization: AUC and Cohen's d of each z-feature vs a binary label.

    Ranking of summaries, not unique CS/CD. AUC uses the raw z-score (not flipped).
    """
    y = np.asarray(y).astype(int)
    rows = []
    for c in cols:
        name = f"z_{c}" if f"z_{c}" in scaled.columns else c
        if name not in scaled.columns:
            continue
        x = scaled[name].to_numpy(dtype=np.float64)
        x = np.nan_to_num(x, nan=0.0)
        if len(np.unique(y)) < 2:
            auc = 0.5
        else:
            auc = float(roc_auc_score(y, x))
        d = cohens_d(x[y == 1], x[y == 0])
        rows.append(
            {
                "feature": c,
                "auc": auc,
                "auc_abs": float(max(auc, 1.0 - auc)),
                "cohens_d": d,
                "mean_y0": float(x[y == 0].mean()) if (y == 0).any() else 0.0,
                "mean_y1": float(x[y == 1].mean()) if (y == 1).any() else 0.0,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("auc_abs", ascending=False).reset_index(drop=True)
        out["rank"] = np.arange(1, len(out) + 1)
    return out


def try_featuretools_dfs(events: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Optional Featuretools EntitySet path. Returns None if the package is absent."""
    try:
        import featuretools as ft  # type: ignore
    except ImportError:
        return None
    try:
        ev = events.copy()
        if "event_id" not in ev.columns:
            ev["event_id"] = np.arange(len(ev), dtype=np.int64)
        users = pd.DataFrame({"user_id": ev["user_id"].drop_duplicates().to_numpy()})
        es = ft.EntitySet(id="tencentgr_seq")
        es = es.add_dataframe(dataframe_name="events", dataframe=ev, index="event_id", time_index="timestamp")
        es = es.add_dataframe(dataframe_name="users", dataframe=users, index="user_id")
        es = es.add_relationship("users", "user_id", "events", "user_id")
        fm, _defs = ft.dfs(
            entityset=es,
            target_dataframe_name="users",
            agg_primitives=["count", "mean", "num_unique", "percent_true", "max", "min"],
            trans_primitives=[],
            max_depth=2,
            verbose=False,
        )
        return fm.reset_index()
    except Exception:
        return None
