#!/usr/bin/env python3
"""买后点击：正向扫一遍的状态机。和 merge_asof 同构。

last-touch（已在 onepass_feats）：CLK 记下 ts，CNV 回头减。
post-cnv：CNV 推进 pending，CLK 把还没闭上的单闭上，并给窗计数 +1。

  C(t)       = #{ clk | ts <= t }
  n_before_W = C(t) - C(t-W)          # (t-W, t]
  n_after_W  = C(t+W) - C(t)          # (t, t+W]
  next_clk   = min { clk.ts | ts > t }   # 严格晚于，同秒不算
  y_W        = 1{dt<=W} ; 跟不满 W 且没点到 → NaN

  python3 scripts/tencent_gr/onepass_post.py
"""
from __future__ import annotations

import bisect
import math
from typing import Dict, List, Optional, Tuple

EXP, CLK, CNV = 0, 1, 2
SESS_GAP = 30 * 60
POST_WINS = {"5m": 5 * 60, "1h": 3600, "1d": 86400, "7d": 7 * 86400}
ATTR_BUCKETS = [(5, "5m"), (30, "30m"), (60, "1h"), (1440, "1d")]
Event = Tuple[int, int, int, Optional[float]]


def _n_in(times: List[int], lo: int, hi: int) -> int:
    """#{ x | lo < x <= hi }，times 升序。"""
    return bisect.bisect_right(times, hi) - bisect.bisect_right(times, lo)


def _y(dt: Optional[int], follow: int, w: int) -> Optional[float]:
    if dt is not None and dt <= w:
        return 1.0
    if follow >= w:
        return 0.0
    return None


def onepass_post(evs: List[Event]) -> List[Dict]:
    """一条用户序列 → 每笔 cnv 一行。先按 (ts, act) 排：同秒 CLK 在 CNV 前。"""
    if not evs:
        return []
    evs = sorted(evs, key=lambda e: (e[2], e[1]))
    t_end = evs[-1][2]

    clk_any: List[int] = []
    clk_item: Dict[int, List[int]] = {}
    last_any: Optional[int] = None
    last_item: Dict[int, int] = {}
    first_item: Dict[int, int] = {}

    prev_ts: Optional[int] = None
    sess_pos = 0
    sess_clk = 0

    pending: List[Dict] = []
    rows: List[Dict] = []

    for iid, act, ts, price in evs:
        if prev_ts is not None and ts - prev_ts > SESS_GAP:
            sess_pos = 0
            sess_clk = 0
        sess_pos += 1

        if act == CLK:
            for p in pending:
                dt = ts - p["cnv_ts"]
                if p["next_clk_ts"] is None:
                    p["next_clk_ts"] = ts
                    p["next_clk_item"] = iid
                if iid == p["item_id"] and p["next_same_ts"] is None:
                    p["next_same_ts"] = ts
                if dt <= 0:
                    continue
                if dt <= SESS_GAP:
                    p["n_clk_after_same_sess"] += 1
                elif dt <= POST_WINS["1d"]:
                    p["n_clk_after_cross_1d"] += 1
                for name, w in POST_WINS.items():
                    if dt <= w:
                        p["n_clk_after"][name] += 1
                        if iid == p["item_id"]:
                            p["n_same_after"][name] += 1
            clk_any.append(ts)
            clk_item.setdefault(iid, []).append(ts)
            last_any = ts
            last_item[iid] = ts
            first_item.setdefault(iid, ts)
            sess_clk += 1

        elif act == CNV:
            rec: Dict = {
                "item_id": iid,
                "cnv_ts": ts,
                "price": price,
                "sess_pos": sess_pos,
                "sess_clk_before": sess_clk,
                "n_prior_cnv": len(rows),
                "wo_prior_clk": iid not in last_item,
                "dt_item_min": (ts - last_item[iid]) / 60.0 if iid in last_item else float("nan"),
                "dt_any_min": (ts - last_any) / 60.0 if last_any is not None else float("nan"),
                "dt_first_min": (ts - first_item[iid]) / 60.0 if iid in first_item else float("nan"),
                "n_clk_same_before": _n_in(clk_item.get(iid, []), -1, ts),
                "n_clk_before": {},
                "n_same_before": {},
                "n_clk_after": {k: 0 for k in POST_WINS},
                "n_same_after": {k: 0 for k in POST_WINS},
                "n_clk_after_same_sess": 0,
                "n_clk_after_cross_1d": 0,
                "next_clk_ts": None,
                "next_clk_item": None,
                "next_same_ts": None,
            }
            for name, w in POST_WINS.items():
                rec["n_clk_before"][name] = _n_in(clk_any, ts - w, ts)
                rec["n_same_before"][name] = _n_in(clk_item.get(iid, []), ts - w, ts)
            if iid in last_item:
                dtm = rec["dt_item_min"]
                for b, name in ATTR_BUCKETS:
                    rec[f"item_within_{name}"] = bool(dtm <= b)
            else:
                for _, name in ATTR_BUCKETS:
                    rec[f"item_within_{name}"] = False
            pending.append(rec)
            rows.append(rec)

        prev_ts = ts

    # 收口：y / lift / 删失。lag 只用此前已经收口的 y_1d。
    y1d_cnt = 0
    y1d_sum = 0.0
    for rec in rows:
        follow = t_end - rec["cnv_ts"]
        dt_any = None if rec["next_clk_ts"] is None else rec["next_clk_ts"] - rec["cnv_ts"]
        dt_same = None if rec["next_same_ts"] is None else rec["next_same_ts"] - rec["cnv_ts"]
        rec["dt_next_clk_sec"] = dt_any
        rec["dt_next_same_sec"] = dt_same
        rec["dt_next_clk_min"] = None if dt_any is None else dt_any / 60.0
        rec["next_clk_same_sess"] = None if dt_any is None else float(dt_any <= SESS_GAP)
        rec["next_clk_cross_sess"] = None if dt_any is None else float(dt_any > SESS_GAP)
        rec["log1p_price"] = math.log1p(max(float(rec["price"] or 0.0), 0.0))
        if follow >= SESS_GAP:
            rec["n_clk_after_same_sess"] = rec["n_clk_after_same_sess"]
        else:
            rec["n_clk_after_same_sess"] = float("nan")
        if follow >= POST_WINS["1d"]:
            rec["n_clk_after_cross_1d"] = rec["n_clk_after_cross_1d"]
        else:
            rec["n_clk_after_cross_1d"] = float("nan")
        rec["lag_post_clk_1d_rate"] = (y1d_sum / y1d_cnt) if y1d_cnt else float("nan")

        for name, w in POST_WINS.items():
            rec[f"y_post_clk_{name}"] = _y(dt_any, follow, w)
            rec[f"y_post_same_{name}"] = _y(dt_same, follow, w)
            rec[f"n_clk_before_{name}"] = rec["n_clk_before"][name]
            rec[f"n_same_before_{name}"] = rec["n_same_before"][name]
            if follow >= w:
                rec[f"n_clk_after_{name}"] = rec["n_clk_after"][name]
                rec[f"n_same_after_{name}"] = rec["n_same_after"][name]
                rec[f"lift_{name}"] = rec[f"n_clk_after_{name}"] / (rec[f"n_clk_before_{name}"] + 1.0)
                rec[f"delta_{name}"] = rec[f"n_clk_after_{name}"] - rec[f"n_clk_before_{name}"]
                rec[f"lift_same_{name}"] = rec[f"n_same_after_{name}"] / (rec[f"n_same_before_{name}"] + 1.0)
            else:
                rec[f"n_clk_after_{name}"] = float("nan")
                rec[f"n_same_after_{name}"] = float("nan")
                rec[f"lift_{name}"] = float("nan")
                rec[f"delta_{name}"] = float("nan")
                rec[f"lift_same_{name}"] = float("nan")

        y1d = rec["y_post_clk_1d"]
        if y1d is not None and y1d == y1d:
            y1d_sum += float(y1d)
            y1d_cnt += 1
    return rows


def user_post(rows: List[Dict]) -> Dict[str, float]:
    """用户一行：未删失 y 的均值；lift / dt 的中位数。"""
    import numpy as np

    def mean_obs(key: str) -> float:
        xs = [r[key] for r in rows if r.get(key) is not None and r[key] == r[key]]
        return float(np.mean(xs)) if xs else 0.0

    def med(key: str, empty: float) -> float:
        xs = [r[key] for r in rows if r.get(key) is not None and r[key] == r[key]]
        return float(np.median(xs)) if xs else empty

    if not rows:
        return {
            "post_n_obs_1d": 0.0,
            "post_clk_1d_rate": 0.0,
            "post_same_1d_rate": 0.0,
            "post_lift_1d_p50": 0.0,
            "post_clk_dt_p50": -1.0,
            "post_clk_same_sess_rate": 0.0,
            "post_clk_cross_sess_rate": 0.0,
            "post_n_clk_after_same_sess_mean": 0.0,
            "post_n_clk_after_cross_1d_mean": 0.0,
        }
    return {
        "post_n_obs_1d": float(sum(r["y_post_clk_1d"] == r["y_post_clk_1d"] for r in rows if r["y_post_clk_1d"] is not None)),
        "post_clk_5m_rate": mean_obs("y_post_clk_5m"),
        "post_clk_1h_rate": mean_obs("y_post_clk_1h"),
        "post_clk_1d_rate": mean_obs("y_post_clk_1d"),
        "post_clk_7d_rate": mean_obs("y_post_clk_7d"),
        "post_same_1h_rate": mean_obs("y_post_same_1h"),
        "post_same_1d_rate": mean_obs("y_post_same_1d"),
        "post_clk_same_sess_rate": mean_obs("next_clk_same_sess"),
        "post_clk_cross_sess_rate": mean_obs("next_clk_cross_sess"),
        "post_clk_dt_p50": med("dt_next_clk_min", -1.0),
        "post_lift_1d_p50": med("lift_1d", 0.0),
        "post_delta_1d_p50": med("delta_1d", 0.0),
        "post_n_clk_after_1d_mean": mean_obs("n_clk_after_1d"),
        "post_n_clk_before_1d_mean": mean_obs("n_clk_before_1d"),
        "post_n_clk_after_same_sess_mean": mean_obs("n_clk_after_same_sess"),
        "post_n_clk_after_cross_1d_mean": mean_obs("n_clk_after_cross_1d"),
    }


def _close(a, b) -> bool:
    def nanish(x) -> bool:
        if x is None:
            return True
        try:
            return x != x
        except Exception:
            return False

    if nanish(a) and nanish(b):
        return True
    if nanish(a) or nanish(b):
        return False
    try:
        return abs(float(a) - float(b)) < 1e-9
    except (TypeError, ValueError):
        return bool(a) == bool(b)


def assert_match_asof(evs: List[Event], uid: int = 1) -> None:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from block_tables import events_frame, tab_attr_events, tab_post_cnv_events

    rows = onepass_post(evs)
    ev = events_frame(uid, evs)
    post = tab_post_cnv_events(ev, tab_attr_events(ev))
    assert len(rows) == len(post), (len(rows), len(post))
    keys = [
        "y_post_clk_1h",
        "y_post_clk_1d",
        "y_post_same_1d",
        "n_clk_before_1d",
        "n_clk_after_1d",
        "lift_1d",
        "n_clk_same_before",
        "dt_item_min",
        "dt_any_min",
        "sess_pos",
        "sess_clk_before",
        "lag_post_clk_1d_rate",
        "wo_prior_clk",
        "next_clk_cross_sess",
    ]
    for i, r in enumerate(rows):
        for k in keys:
            av, bv = r[k], post.iloc[i][k]
            if not _close(av, bv):
                raise AssertionError(f"row{i} {k}: sm={av!r} asof={bv!r}")


def _demo() -> None:
    t0 = 1_700_000_000
    evs: List[Event] = [
        (11, EXP, t0, 99.0),
        (11, CLK, t0 + 20, 99.0),
        (11, CNV, t0 + 3600, 99.0),
        (11, CLK, t0 + 3600 + 600, 99.0),
        (22, EXP, t0 + 86400 + 10, 50.0),
        (22, CLK, t0 + 86400 + 30, 50.0),
        (22, CNV, t0 + 86400 + 40, 50.0),
        (99, EXP, t0 + 20 * 86400, None),
    ]
    rows = onepass_post(evs)
    keys = [
        "item_id",
        "y_post_clk_1h",
        "y_post_clk_1d",
        "n_clk_before_1d",
        "n_clk_after_1d",
        "lift_1d",
        "next_clk_same_sess",
        "lag_post_clk_1d_rate",
        "dt_item_min",
        "sess_clk_before",
    ]
    for r in rows:
        print(" ".join(f"{k}={r[k]}" for k in keys))
    u = user_post(rows)
    print("user", {k: round(v, 4) for k, v in u.items()})
    assert rows[0]["y_post_clk_1h"] == 1.0
    assert rows[1]["y_post_clk_1d"] == 0.0
    assert rows[0]["n_clk_after_1d"] == 2
    assert rows[0]["n_clk_before_1d"] == 1
    assert rows[1]["lag_post_clk_1d_rate"] == 1.0
    assert rows[0]["next_clk_same_sess"] == 1.0
    assert rows[0]["next_clk_cross_sess"] == 0.0
    assert rows[0]["n_clk_after_same_sess"] == 1
    assert rows[0]["n_clk_after_cross_1d"] == 1
    assert_match_asof(evs)

    rng = __import__("random").Random(0)
    t = t0
    rnd: List[Event] = []
    for _ in range(80):
        t += rng.choice([10, 60, 400, 2000, 9000, 86400])
        rnd.append((rng.choice([11, 22, 33]), rng.choice([EXP, CLK, CNV]), t, 10.0))
    rnd.append((0, EXP, t + 30 * 86400, None))
    assert_match_asof(rnd)
    print("asof match ok")


if __name__ == "__main__":
    _demo()
