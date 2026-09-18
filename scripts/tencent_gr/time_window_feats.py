#!/usr/bin/env python3
"""Time-windowed 图谱特征 engineering (no network / graph algos).

``feature_engineer(root, t_start, t_end, ...)`` only sees events with
``t_start <= timestamp < t_end``. Call twice on distant windows (≥1 month gap)
to build two independent panels, then localize → FSDS.

Outputs per window:
  - edge / user / item tables
  - convert-path localization rows (first / mid / last + linear credit)
  - (user, item) feature grid with ``y_convert`` = in-window conversion
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

A_EXP, A_CLK, A_CNV = 0, 1, 2
Event = Tuple[int, int, int]  # item_id, action_type, timestamp
DAY = 86400
# TencentGR-1M public logs label exposure(0)/click(1) only — no action=2.
# Terminal success for localization / y defaults to click; override if cnv exists.
DEFAULT_TERMINAL_ACTION = A_CLK


@dataclass
class TimeWindow:
    name: str
    t_start: int
    t_end: int

    @property
    def n_days(self) -> float:
        return max(0.0, (self.t_end - self.t_start) / DAY)

    def to_dict(self) -> dict:
        return asdict(self)


def _pq(d: Path) -> List[Path]:
    return sorted(p for p in d.glob("*.parquet") if p.suffix == ".parquet" or p.name.endswith(".parquet"))


def parse_events(seq: list) -> List[Event]:
    evs: List[Event] = []
    for e in seq:
        if not isinstance(e, dict):
            e = dict(e)
        act = e.get("action_type")
        if act is None:
            continue
        try:
            evs.append((int(e["item_id"]), int(float(act)), int(e["timestamp"])))
        except (TypeError, ValueError):
            continue
    evs.sort(key=lambda x: x[2])
    return evs


def filter_events(evs: Sequence[Event], t_start: int, t_end: int) -> List[Event]:
    return [e for e in evs if t_start <= e[2] < t_end]


def path_items(evs: Sequence[Event]) -> List[int]:
    out: List[int] = []
    for iid, _, _ in evs:
        if not out or out[-1] != iid:
            out.append(iid)
    return out


def rate(n: float, d: float) -> float:
    return float(n) / float(d) if d > 0 else 0.0


def scan_time_range(root: Path, *, max_users: Optional[int] = None) -> Tuple[int, int, int]:
    """Return (t_min, t_max, n_users_seen)."""
    seq_dir = root / "seq"
    t_min, t_max = None, None
    n = 0
    for path in _pq(seq_dir):
        pf = pq.ParquetFile(path)
        for rg in range(pf.num_row_groups):
            df = pf.read_row_group(rg, columns=["user_id", "seq"]).to_pandas()
            for seq in df["seq"]:
                for e in seq:
                    e = dict(e)
                    t = int(e["timestamp"])
                    t_min = t if t_min is None else min(t_min, t)
                    t_max = t if t_max is None else max(t_max, t)
                n += 1
                if max_users is not None and n >= max_users:
                    return int(t_min), int(t_max), n
    if t_min is None:
        raise FileNotFoundError(f"no seq under {seq_dir}")
    return int(t_min), int(t_max), n


def propose_two_windows(
    t_min: int,
    t_max: int,
    *,
    window_days: int = 45,
    gap_days: int = 30,
) -> Tuple[TimeWindow, TimeWindow]:
    """Early / late windows with ≥ ``gap_days`` between them.

    Anchors the late window near ``t_max`` (denser tail on TencentGR-1M).
    """
    span = t_max - t_min
    need = (window_days * 2 + gap_days) * DAY
    if span < need:
        rem = span - gap_days * DAY
        if rem < 2 * DAY:
            raise ValueError(
                f"timeline too short for gap={gap_days}d: span={span/DAY:.1f}d"
            )
        window_days = max(1, int(rem / (2 * DAY)))
    w = window_days * DAY
    g = gap_days * DAY
    w2 = TimeWindow("W2_late", int(t_max - w), int(t_max))
    w1 = TimeWindow("W1_early", int(w2.t_start - g - w), int(w2.t_start - g))
    if w1.t_start < t_min:
        w1 = TimeWindow("W1_early", int(t_min), int(t_min + w))
        w2 = TimeWindow("W2_late", int(w1.t_end + g), int(w1.t_end + g + w))
    if w2.t_start - w1.t_end < gap_days * DAY - 1:
        raise ValueError("failed to place windows with required gap")
    return w1, w2


def iter_users(root: Path, max_users: int) -> Iterable[Tuple[int, List[Event]]]:
    n = 0
    for path in _pq(root / "seq"):
        pf = pq.ParquetFile(path)
        for rg in range(pf.num_row_groups):
            df = pf.read_row_group(rg, columns=["user_id", "seq"]).to_pandas()
            for _, row in df.iterrows():
                yield int(row["user_id"]), parse_events(list(row["seq"]))
                n += 1
                if n >= max_users:
                    return


def feature_engineer(
    root: Path,
    t_start: int,
    t_end: int,
    *,
    max_users: int = 4000,
    co_window: int = 8,
    top_covisit: int = 10,
    window_name: str = "window",
    terminal_action: int = DEFAULT_TERMINAL_ACTION,
) -> Dict:
    """图谱特征 FE restricted to ``[t_start, t_end)``.

    No events outside the range enter aggregates, localization, or labels.

    ``terminal_action``: path success used for 三段启发式 + ``y_convert``
    (click=1 on TencentGR-1M; set 2 if conversion codes exist).
    """
    if t_end <= t_start:
        raise ValueError("t_end must be > t_start")

    user_rows: List[dict] = []
    item_exp: Counter = Counter()
    item_clk: Counter = Counter()
    item_cnv: Counter = Counter()
    item_users: Dict[int, set] = defaultdict(set)
    item_first: Counter = Counter()
    item_last: Counter = Counter()
    item_linear: Counter = Counter()
    item_as_cnv: Counter = Counter()
    ui_exp: Counter = Counter()
    ui_clk: Counter = Counter()
    ui_cnv: Counter = Counter()
    ui_last_ts: Dict[Tuple[int, int], int] = {}
    convert_rows: List[dict] = []
    covisit: Dict[int, Counter] = defaultdict(Counter)

    n_users = 0
    n_users_active = 0
    n_conv_users = 0
    n_events = 0
    term = int(terminal_action)

    for uid, evs_all in iter_users(root, max_users):
        n_users += 1
        evs = filter_events(evs_all, t_start, t_end)
        if not evs:
            continue
        n_users_active += 1
        n_events += len(evs)
        n_exp = n_clk = n_term = 0
        items_seen: List[int] = []
        for iid, act, ts in evs:
            items_seen.append(iid)
            item_users[iid].add(uid)
            key = (uid, iid)
            ui_last_ts[key] = ts
            if act == A_EXP:
                n_exp += 1
                item_exp[iid] += 1
                ui_exp[key] += 1
            elif act == A_CLK:
                n_clk += 1
                item_clk[iid] += 1
                ui_clk[key] += 1
            elif act == A_CNV:
                item_cnv[iid] += 1
            if act == term:
                n_term += 1
                ui_cnv[key] += 1  # reuse ui_cnv as terminal-success counter

        user_rows.append(
            {
                "user_id": uid,
                "n_events": len(evs),
                "n_exp": n_exp,
                "n_clk": n_clk,
                "n_cnv": n_term,  # terminal successes in-window
                "n_uniq_items": len(set(items_seen)),
                "ctr": rate(n_clk, n_exp),
                "cvr": rate(n_term, n_clk) if term != A_CLK else rate(n_clk, n_exp),
                "has_convert": int(n_term > 0),
                "span_sec": int(evs[-1][2] - evs[0][2]) if len(evs) > 1 else 0,
            }
        )

        path = path_items(evs)
        for i, a in enumerate(path):
            for b in path[i + 1 : i + 1 + co_window]:
                if a == b:
                    continue
                covisit[a][b] += 1
                covisit[b][a] += 1

        if n_term > 0:
            n_conv_users += 1
            last_term_i = max(i for i, (_, a, _) in enumerate(evs) if a == term)
            term_item = evs[last_term_i][0]
            prefix = path_items(evs[: last_term_i + 1])
            if not prefix:
                continue
            item_as_cnv[term_item] += 1
            item_first[prefix[0]] += 1.0
            item_last[prefix[-1]] += 1.0
            w = 1.0 / float(len(prefix))
            for pos, iid in enumerate(prefix):
                item_linear[iid] += w
                if pos == 0:
                    touch = "first"
                elif pos == len(prefix) - 1:
                    touch = "last"
                else:
                    touch = "mid"
                convert_rows.append(
                    {
                        "user_id": uid,
                        "convert_item": term_item,
                        "path_item": iid,
                        "path_pos": pos,
                        "path_len": len(prefix),
                        "linear_credit": w,
                        "touch": touch,
                        "is_convert_item": int(iid == term_item),
                    }
                )

    s_first = sum(item_first.values()) or 1.0
    s_last = sum(item_last.values()) or 1.0
    s_lin = sum(item_linear.values()) or 1.0
    all_items = sorted(
        set(item_exp) | set(item_clk) | set(item_cnv) | set(item_linear) | set(item_users)
    )
    item_rows = []
    for iid in all_items:
        ne, nc = item_exp[iid], item_clk[iid]
        n_term_item = item_clk[iid] if term == A_CLK else item_cnv[iid]
        item_rows.append(
            {
                "item_id": iid,
                "n_exp": ne,
                "n_clk": nc,
                "n_cnv": n_term_item,
                "n_users": len(item_users[iid]),
                "ctr": rate(nc, ne),
                "cvr": rate(n_term_item, nc) if term != A_CLK else rate(nc, ne),
                "n_as_convert_terminal": item_as_cnv[iid],
                "credit_first": float(item_first[iid]),
                "credit_last": float(item_last[iid]),
                "credit_linear": float(item_linear[iid]),
                "share_first": float(item_first[iid]) / s_first,
                "share_last": float(item_last[iid]) / s_last,
                "share_linear": float(item_linear[iid]) / s_lin,
                "n_covisit_neighbors": len(covisit[iid]),
            }
        )
    item_df = pd.DataFrame(item_rows)
    user_df = pd.DataFrame(user_rows)
    if len(item_df):
        item_df["log1p_n_exp"] = np.log1p(item_df["n_exp"])
        item_df["log1p_n_clk"] = np.log1p(item_df["n_clk"])
        item_df["log1p_n_users"] = np.log1p(item_df["n_users"])
        item_df["log1p_n_covisit"] = np.log1p(item_df["n_covisit_neighbors"])
        item_df["item_credit_rank"] = item_df["share_linear"].rank(
            method="average", ascending=False
        )
        item_df["item_pop_rank"] = item_df["n_users"].rank(method="average", ascending=False)

    if len(user_df):
        user_df["log1p_n_events"] = np.log1p(user_df["n_events"])
        user_df["log1p_n_uniq"] = np.log1p(user_df["n_uniq_items"])
        user_df["user_activity_rank"] = user_df["n_events"].rank(
            method="average", ascending=False
        )

    # edges: label = in-window convert; FE features for modeling drop raw cnv counts
    keys = set(ui_exp) | set(ui_clk) | set(ui_cnv)
    edge_rows = []
    for uid, iid in keys:
        edge_rows.append(
            {
                "user_id": uid,
                "item_id": iid,
                "n_exp": ui_exp.get((uid, iid), 0),
                "n_clk": ui_clk.get((uid, iid), 0),
                "n_cnv": ui_cnv.get((uid, iid), 0),
                "last_ts": ui_last_ts.get((uid, iid), 0),
            }
        )
    edge_df = pd.DataFrame(edge_rows)
    if len(edge_df):
        edge_df["ctr"] = edge_df.apply(lambda r: rate(r["n_clk"], r["n_exp"]), axis=1)
        # CTR-style label: only edges with ≥1 exposure; y = terminal success.
        # Dropping click-only edges avoids the trivial e_n_exp=0 ⇒ y=1 leak.
        if term == A_CLK:
            edge_df = edge_df[edge_df["n_exp"] > 0].copy()
            edge_df["y_convert"] = (edge_df["n_clk"] > 0).astype(int)
        else:
            edge_df["y_convert"] = (edge_df["n_cnv"] > 0).astype(int)

    covisit_rows = []
    for a, ctr in covisit.items():
        for b, c in ctr.most_common(top_covisit):
            covisit_rows.append(
                {"item_id": a, "neighbor_item": b, "covisit_count": int(c)}
            )
    covisit_df = pd.DataFrame(covisit_rows)
    convert_df = pd.DataFrame(convert_rows)

    grid = build_grid(user_df, item_df, edge_df)

    meta = {
        "window_name": window_name,
        "t_start": int(t_start),
        "t_end": int(t_end),
        "n_days": float((t_end - t_start) / DAY),
        "terminal_action": term,
        "terminal_name": {0: "exp", 1: "click", 2: "cnv"}.get(term, str(term)),
        "max_users_scan": max_users,
        "n_users_scanned": n_users,
        "n_users_active": n_users_active,
        "n_conv_users": n_conv_users,
        "n_events": n_events,
        "n_items": int(len(item_df)),
        "n_edges": int(len(edge_df)),
        "n_grid": int(len(grid)),
        "pos_rate": float(edge_df["y_convert"].mean()) if len(edge_df) else float("nan"),
        "co_window": co_window,
    }
    return {
        "user_df": user_df,
        "item_df": item_df,
        "edge_df": edge_df,
        "convert_df": convert_df,
        "covisit_df": covisit_df,
        "grid": grid,
        "meta": meta,
    }


def build_grid(
    user_df: pd.DataFrame, item_df: pd.DataFrame, edge_df: pd.DataFrame
) -> pd.DataFrame:
    if edge_df is None or len(edge_df) == 0:
        return pd.DataFrame()
    e = edge_df.rename(
        columns={
            "n_exp": "e_n_exp",
            "n_clk": "e_n_clk",
            "n_cnv": "e_n_cnv",
            "ctr": "e_ctr",
            "last_ts": "e_last_ts",
        }
    )
    u = user_df.add_prefix("u_").rename(columns={"u_user_id": "user_id"})
    i = item_df.add_prefix("i_").rename(columns={"i_item_id": "item_id"})
    g = e.merge(u, on="user_id", how="left").merge(i, on="item_id", how="left")
    g["e_log1p_exp"] = np.log1p(g["e_n_exp"])
    g["e_log1p_clk"] = np.log1p(g["e_n_clk"])
    if "i_item_pop_rank" in g.columns and "u_user_activity_rank" in g.columns:
        g["ui_pop_mismatch"] = (
            g["i_item_pop_rank"].fillna(0) - g["u_user_activity_rank"].fillna(0)
        ).astype(float)
    return g


def localize_item_subset(item_df: pd.DataFrame, *, top_k: int = 200) -> List[int]:
    """Top-k items by linear credit share (localization subset)."""
    if item_df is None or len(item_df) == 0:
        return []
    sub = item_df.sort_values("share_linear", ascending=False).head(top_k)
    return [int(x) for x in sub["item_id"].tolist()]


# Hard leaks for FSDS on y_convert / terminal success
LEAK_EXACT = {
    "y_convert",
    "e_n_cnv",
    "e_n_clk",
    "e_ctr",
    "e_log1p_clk",
    "u_n_cnv",
    "u_n_clk",
    "u_has_convert",
    "u_cvr",
    "u_ctr",
    "i_n_cnv",
    "i_n_clk",
    "i_ctr",
    "i_cvr",
    "i_n_as_convert_terminal",
}
LEAK_SUBSTR = (
    "n_cnv",
    "n_clk",
    "cvr",
    "ctr",
    "ctcvr",
    "has_convert",
    "convert_terminal",
    "y_convert",
    "log1p_clk",
)


def fsds_feature_columns(grid: pd.DataFrame) -> List[str]:
    cols = []
    for c in grid.columns:
        if c in ("user_id", "item_id", "e_last_ts"):
            continue
        if c in LEAK_EXACT:
            continue
        if any(s in c for s in LEAK_SUBSTR):
            continue
        if not np.issubdtype(grid[c].dtype, np.number):
            continue
        cols.append(c)
    return cols
