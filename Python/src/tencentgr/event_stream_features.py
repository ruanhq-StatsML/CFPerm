"""Request-level features on a long event table.

Input columns (minimum): user_id, item_id, action_type, timestamp
Optional: item side info joined beforehand (item_format from field 100, etc.)

All window counts use events *strictly before* the current row (no leakage).
y_ctr / y_cvr look *strictly after* the current row (labels, not features).

Not unique CS/CD. Not CATE. time_late is still a batch index if you join it later.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

ACTION_EXPOSURE = 0
ACTION_CLICK = 1
ACTION_CONVERSION = 2

HOUR = 3600
DAY = 86400
WEEK = 7 * DAY
DEFAULT_SESSION_GAP = 30 * 60  # 30 min silent → new session
CTR_WINDOW = DAY
CVR_WINDOW = 7 * DAY


def _as_int_action(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(-1).astype(np.int16)


def prepare_events(
    df: pd.DataFrame,
    *,
    item_feat: Optional[pd.DataFrame] = None,
    format_col: str = "100",
) -> pd.DataFrame:
    """Normalize the exploded table. Does not require pos / event_id."""
    ev = df.copy()
    need = {"user_id", "item_id", "action_type", "timestamp"}
    missing = need - set(ev.columns)
    if missing:
        raise ValueError(f"need columns {sorted(need)}, missing {sorted(missing)}")
    ev["user_id"] = ev["user_id"].astype(np.int64)
    ev["item_id"] = ev["item_id"].astype(np.int64)
    ev["action_type"] = _as_int_action(ev["action_type"])
    ev["timestamp"] = ev["timestamp"].astype(np.int64)
    if "time" not in ev.columns:
        ev["time"] = pd.to_datetime(ev["timestamp"], unit="s", utc=True)
    else:
        ev["time"] = pd.to_datetime(ev["time"], utc=True, errors="coerce")
        if ev["time"].isna().any():
            ev.loc[ev["time"].isna(), "time"] = pd.to_datetime(
                ev.loc[ev["time"].isna(), "timestamp"], unit="s", utc=True
            )
    ev["hour"] = ev["time"].dt.hour.astype(np.int8)
    ev["dow"] = ev["time"].dt.dayofweek.astype(np.int8)
    ev["is_click"] = (ev["action_type"] == ACTION_CLICK).astype(np.int8)
    ev["is_conversion"] = (ev["action_type"] == ACTION_CONVERSION).astype(np.int8)
    ev["is_exposure"] = (ev["action_type"] == ACTION_EXPOSURE).astype(np.int8)
    ev["is_engage"] = ((ev["is_click"] == 1) | (ev["is_conversion"] == 1)).astype(np.int8)
    if item_feat is not None and "item_id" in item_feat.columns and format_col in item_feat.columns:
        fmt = item_feat[["item_id", format_col]].drop_duplicates("item_id")
        fmt = fmt.rename(columns={format_col: "item_format"})
        ev = ev.merge(fmt, on="item_id", how="left")
    elif "item_format" not in ev.columns:
        ev["item_format"] = np.nan
    if "event_id" not in ev.columns:
        ev["event_id"] = np.arange(len(ev), dtype=np.int64)
    return ev


def _window_prev_sum(sorted_df: pd.DataFrame, group_cols: Sequence[str], flag: str, window: int) -> np.ndarray:
    """Two-pointer: sum of `flag` in (t - window, t) within each group. Current row excluded."""
    n = len(sorted_df)
    out = np.zeros(n, dtype=np.float64)
    if n == 0:
        return out
    ts = sorted_df["timestamp"].to_numpy(dtype=np.int64)
    fl = sorted_df[flag].to_numpy(dtype=np.float64)
    keys = [sorted_df[c].to_numpy() for c in group_cols]
    start = 0
    while start < n:
        end = start + 1
        while end < n and all(keys[j][end] == keys[j][start] for j in range(len(keys))):
            end += 1
        left = start
        acc = 0.0
        for i in range(start, end):
            while ts[i] - ts[left] > window:
                acc -= fl[left]
                left += 1
            out[i] = acc
            acc += fl[i]
        start = end
    return out


def add_gaps_and_sessions(ev: pd.DataFrame, session_gap_sec: int = DEFAULT_SESSION_GAP) -> pd.DataFrame:
    """Inter-event gap, session cut, position in session. Gap on the first event of a user is NaN."""
    d = ev.sort_values(["user_id", "timestamp", "event_id"], kind="mergesort").copy()
    d["gap_sec"] = d.groupby("user_id", sort=False)["timestamp"].diff()
    prev_a = d.groupby("user_id", sort=False)["action_type"].shift(1)
    d["prev_action"] = prev_a
    d["gap_exp_to_click_sec"] = np.where(
        (prev_a == ACTION_EXPOSURE) & (d["action_type"] == ACTION_CLICK),
        d["gap_sec"],
        np.nan,
    )
    new_sess = d["gap_sec"].isna() | (d["gap_sec"] > float(session_gap_sec))
    d["session_id"] = new_sess.groupby(d["user_id"], sort=False).cumsum().astype(np.int32)
    d["pos_in_session"] = d.groupby(["user_id", "session_id"], sort=False).cumcount().astype(np.int32)
    d["session_len"] = d.groupby(["user_id", "session_id"], sort=False)["event_id"].transform("size")
    return d


def add_recency(ev: pd.DataFrame) -> pd.DataFrame:
    """Seconds since last click / conversion / same item. First time → NaN (not zero)."""
    d = ev.sort_values(["user_id", "timestamp", "event_id"], kind="mergesort").copy()
    last_click = np.where(d["is_click"] == 1, d["timestamp"].to_numpy(), np.nan)
    last_conv = np.where(d["is_conversion"] == 1, d["timestamp"].to_numpy(), np.nan)
    # shift so current click does not zero out its own recency
    d["_lc"] = pd.Series(last_click, index=d.index).groupby(d["user_id"]).ffill()
    d["_lc"] = d.groupby("user_id")["_lc"].shift(1)
    d["_lv"] = pd.Series(last_conv, index=d.index).groupby(d["user_id"]).ffill()
    d["_lv"] = d.groupby("user_id")["_lv"].shift(1)
    d["sec_since_click"] = d["timestamp"] - d["_lc"]
    d["sec_since_conv"] = d["timestamp"] - d["_lv"]
    last_item = d.groupby(["user_id", "item_id"], sort=False)["timestamp"].shift(1)
    d["sec_since_same_item"] = d["timestamp"] - last_item
    d = d.drop(columns=["_lc", "_lv"])
    return d


def add_fatigue(
    ev: pd.DataFrame,
    windows: Iterable[tuple[str, int]] = (
        ("1h", HOUR),
        ("1d", DAY),
        ("7d", WEEK),
    ),
) -> pd.DataFrame:
    """Frequency in lookback windows, current row excluded.

    item_exp_*  same user × item exposures (creative-level fatigue)
    fmt_exp_*   same user × item_format (morphology fatigue; needs item_format)
    user_exp_*  all exposures for the user (session pressure)
    user_clk_*  clicks in the window (recent engagement, not the label)
    """
    d = ev.sort_values(["user_id", "timestamp", "event_id"], kind="mergesort").copy()
    idx = d.index.to_numpy()
    work = d.reset_index(drop=True)
    for suffix, w in windows:
        work[f"item_exp_{suffix}"] = _window_prev_sum(work, ["user_id", "item_id"], "is_exposure", w)
        work[f"user_exp_{suffix}"] = _window_prev_sum(work, ["user_id"], "is_exposure", w)
        work[f"user_clk_{suffix}"] = _window_prev_sum(work, ["user_id"], "is_click", w)
        if work["item_format"].notna().any():
            filled = work.copy()
            filled["item_format"] = filled["item_format"].where(filled["item_format"].notna(), other=-1)
            work[f"fmt_exp_{suffix}"] = _window_prev_sum(filled, ["user_id", "item_format"], "is_exposure", w)
        else:
            work[f"fmt_exp_{suffix}"] = 0.0
    work.index = idx
    return work


def add_multitask_labels(
    ev: pd.DataFrame,
    *,
    ctr_window_sec: int = CTR_WINDOW,
    cvr_window_sec: int = CVR_WINDOW,
) -> pd.DataFrame:
    """Same-item future click / conversion. Attribution windows, not CATE.

    y_ctr: exposure (or any row) followed by a click on the same item within ctr_window.
    y_cvr: followed by a conversion on the same item within cvr_window.
    On a click row, y_cvr is the ESMM-style conversion head (post-click).
    """
    d = ev.sort_values(["user_id", "item_id", "timestamp", "event_id"], kind="mergesort").copy()

    def _next_time(mask: pd.Series) -> pd.Series:
        t = np.where(mask.to_numpy(), d["timestamp"].to_numpy(dtype=np.float64), np.inf)
        s = pd.Series(t, index=d.index)
        rev = s.iloc[::-1]
        g_user = d["user_id"].iloc[::-1]
        g_item = d["item_id"].iloc[::-1]
        cum = rev.groupby([g_user, g_item], sort=False).cummin()
        nxt = cum.iloc[::-1]
        nxt = nxt.groupby([d["user_id"], d["item_id"]], sort=False).shift(-1)
        return nxt.replace([np.inf, -np.inf], np.nan)

    d["t_next_click"] = _next_time(d["is_click"] == 1)
    d["t_next_conv"] = _next_time(d["is_conversion"] == 1)
    dc = d["t_next_click"] - d["timestamp"]
    dv = d["t_next_conv"] - d["timestamp"]
    d["y_ctr"] = ((dc > 0) & (dc <= ctr_window_sec)).astype(np.int8)
    d["y_cvr"] = ((dv > 0) & (dv <= cvr_window_sec)).astype(np.int8)
    d["y_cvr_on_click"] = np.where(d["is_click"] == 1, d["y_cvr"], np.nan)
    return d


def split_behavior_channels(ev: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Multi-behavior: three parallel sequences, same columns."""
    return {
        "exposure": ev.loc[ev["action_type"] == ACTION_EXPOSURE].copy(),
        "click": ev.loc[ev["action_type"] == ACTION_CLICK].copy(),
        "conversion": ev.loc[ev["action_type"] == ACTION_CONVERSION].copy(),
    }


def annotate_event_stream(
    df: pd.DataFrame,
    *,
    item_feat: Optional[pd.DataFrame] = None,
    session_gap_sec: int = DEFAULT_SESSION_GAP,
    ctr_window_sec: int = CTR_WINDOW,
    cvr_window_sec: int = CVR_WINDOW,
) -> pd.DataFrame:
    """Full pass: prepare → gap/session → recency → fatigue → multi-task labels."""
    ev = prepare_events(df, item_feat=item_feat)
    ev = add_gaps_and_sessions(ev, session_gap_sec=session_gap_sec)
    ev = add_recency(ev)
    ev = add_fatigue(ev)
    ev = add_multitask_labels(ev, ctr_window_sec=ctr_window_sec, cvr_window_sec=cvr_window_sec)
    return ev.sort_values(["user_id", "timestamp", "event_id"], kind="mergesort").reset_index(drop=True)


def exposure_rank_table(ev: pd.DataFrame) -> pd.DataFrame:
    """Rows that look like pCTR samples: exposures only, with features + y_ctr/y_cvr."""
    cols = [c for c in ev.columns if c.startswith(("item_exp_", "fmt_exp_", "user_exp_", "user_clk_")) or c in {
        "user_id", "item_id", "timestamp", "time", "hour", "dow", "item_format",
        "gap_sec", "session_id", "pos_in_session", "session_len",
        "sec_since_click", "sec_since_conv", "sec_since_same_item",
        "y_ctr", "y_cvr",
    }]
    return ev.loc[ev["is_exposure"] == 1, cols].copy()
