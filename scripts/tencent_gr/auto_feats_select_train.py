#!/usr/bin/env python3
"""Auto-generate ~1000 TencentGR seq features → select → train.

Pipeline
--------
1. Expand combinatorial sequence features (~1000 dims):
   multi-window funnel × rates × log1p × decays × session × attribution
   × price/ARPU/deal × Markov × TOD/weekday × side-user.
2. Feature selection: variance filter → MI / F-score → optional L1.
3. Train: HistGradientBoosting (or logreg) on label = pay_user / has_cnv.

Example::

  python3 scripts/tencent_gr/auto_feats_select_train.py \\
    --root data/tencent_subset --max-users 8000 --target-dim 1000 \\
    --select-k 128 --out results/tencent_gr_fs
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import (
    SelectKBest,
    VarianceThreshold,
    f_classif,
    mutual_info_classif,
)
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

A_EXP, A_CLK, A_CNV = 0, 1, 2
A_NAME = {0: "exp", 1: "clk", 2: "cnv"}
# denser windows to explode feature count
WINS = {
    "1h": 3600,
    "6h": 6 * 3600,
    "12h": 12 * 3600,
    "1d": 86400,
    "3d": 3 * 86400,
    "7d": 7 * 86400,
    "14d": 14 * 86400,
    "30d": 30 * 86400,
}
DECAYS = {"hl6h": 6 * 3600, "hl1d": 86400, "hl3d": 3 * 86400, "hl7d": 7 * 86400}
SESS_GAP = 30 * 60
EPS = 1e-6
Event = Tuple[int, int, int, Optional[float]]


def parquets(d: Path) -> List[Path]:
    return sorted(p for p in d.glob("*.parquet") if p.is_file())


def load_price(item_dir: Path) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for p in parquets(item_dir):
        df = pq.read_table(p, columns=["item_id", "115"]).to_pandas()
        for iid, v in zip(df["item_id"].to_numpy(), df["115"].to_numpy()):
            if v is None:
                continue
            try:
                fv = float(v)
            except Exception:
                continue
            if math.isfinite(fv) and fv >= 0:
                out[int(iid)] = fv
    return out


def load_users(user_dir: Path, keep: Optional[set] = None) -> pd.DataFrame:
    parts = []
    for p in parquets(user_dir):
        df = pq.read_table(p).to_pandas()
        if keep is not None:
            df = df[df["user_id"].isin(keep)]
        parts.append(df)
    if not parts:
        return pd.DataFrame(columns=["user_id"])
    return pd.concat(parts, ignore_index=True).drop_duplicates("user_id")


def iter_users(seq_paths: Sequence[Path], max_users: int) -> Iterable[Tuple[int, list]]:
    n = 0
    for path in seq_paths:
        pf = pq.ParquetFile(path)
        for rg in range(pf.num_row_groups):
            df = pf.read_row_group(rg).to_pandas()
            for _, row in df.iterrows():
                yield int(row["user_id"]), list(row["seq"])
                n += 1
                if n >= max_users:
                    return


def parse_seq(seq: list, price: Mapping[int, float]) -> List[Event]:
    evs: List[Event] = []
    for e in seq:
        if not isinstance(e, dict):
            e = dict(e)
        iid = int(e["item_id"])
        act = int(e["action_type"])
        ts = int(e["timestamp"])
        evs.append((iid, act, ts, price.get(iid)))
    evs.sort(key=lambda x: x[2])
    return evs


def rate(n: float, d: float) -> float:
    if d <= 0:
        return 0.0
    return float(n) / float(d)


def entropy(counts: Mapping[Any, int]) -> float:
    tot = float(sum(counts.values()))
    if tot <= 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c <= 0:
            continue
        p = c / tot
        h -= p * math.log(p + EPS)
    return float(h)


def stats(xs: List[float], prefix: str) -> Dict[str, float]:
    if not xs:
        return {
            f"{prefix}_n": 0.0,
            f"{prefix}_min": -1.0,
            f"{prefix}_p50": -1.0,
            f"{prefix}_mean": -1.0,
            f"{prefix}_p90": -1.0,
            f"{prefix}_std": -1.0,
        }
    a = np.asarray(xs, dtype=float)
    return {
        f"{prefix}_n": float(len(xs)),
        f"{prefix}_min": float(a.min()),
        f"{prefix}_p50": float(np.median(a)),
        f"{prefix}_mean": float(a.mean()),
        f"{prefix}_p90": float(np.percentile(a, 90)),
        f"{prefix}_std": float(a.std()),
    }


def sessions(evs: Sequence[Event]) -> List[List[Event]]:
    if not evs:
        return []
    out = [[evs[0]]]
    for e in evs[1:]:
        if e[2] - out[-1][-1][2] > SESS_GAP:
            out.append([e])
        else:
            out[-1].append(e)
    return out


def funnel(evs: Sequence[Event], t_end: int, win: int) -> Dict[str, float]:
    t0 = t_end - win
    sub = [e for e in evs if t0 < e[2] <= t_end]
    n_exp = sum(1 for e in sub if e[1] == A_EXP)
    n_clk = sum(1 for e in sub if e[1] == A_CLK)
    n_cnv = sum(1 for e in sub if e[1] == A_CNV)
    n_tot = len(sub)
    prices = [e[3] for e in sub if e[3] is not None]
    clk_prices = [e[3] for e in sub if e[1] == A_CLK and e[3] is not None]
    cnv_prices = [e[3] for e in sub if e[1] == A_CNV and e[3] is not None]
    uniq = len({e[0] for e in sub})
    return {
        "n_tot": float(n_tot),
        "n_exp": float(n_exp),
        "n_clk": float(n_clk),
        "n_cnv": float(n_cnv),
        "n_uniq_item": float(uniq),
        "ctr": rate(n_clk, n_exp),
        "cvr": rate(n_cnv, n_clk),
        "ctcvr": rate(n_cnv, n_exp),
        "clk_share": rate(n_clk, n_tot),
        "cnv_share": rate(n_cnv, n_tot),
        "uniq_share": rate(uniq, n_tot),
        "price_mean": float(np.mean(prices)) if prices else 0.0,
        "price_p50": float(np.median(prices)) if prices else 0.0,
        "clk_price_mean": float(np.mean(clk_prices)) if clk_prices else 0.0,
        "cnv_price_mean": float(np.mean(cnv_prices)) if cnv_prices else 0.0,
        "cnv_price_sum": float(np.sum(cnv_prices)) if cnv_prices else 0.0,
        "log1p_n_tot": float(math.log1p(n_tot)),
        "log1p_n_clk": float(math.log1p(n_clk)),
        "log1p_n_cnv": float(math.log1p(n_cnv)),
    }


def decay_counts(evs: Sequence[Event], t_end: int, hl: float) -> Dict[str, float]:
    """Exponentially decayed action counts (half-life = hl seconds)."""
    lam = math.log(2.0) / max(hl, 1.0)
    w = {a: 0.0 for a in (A_EXP, A_CLK, A_CNV)}
    w_price_cnv = 0.0
    for iid, act, ts, price in evs:
        age = max(t_end - ts, 0)
        ww = math.exp(-lam * age)
        w[act] += ww
        if act == A_CNV and price is not None:
            w_price_cnv += ww * price
    return {
        "dec_exp": w[A_EXP],
        "dec_clk": w[A_CLK],
        "dec_cnv": w[A_CNV],
        "dec_ctr": rate(w[A_CLK], w[A_EXP]),
        "dec_cvr": rate(w[A_CNV], w[A_CLK]),
        "dec_ctcvr": rate(w[A_CNV], w[A_EXP]),
        "dec_arpu": w_price_cnv,
    }


def attribution(evs: Sequence[Event]) -> Dict[str, float]:
    last_exp: Dict[int, int] = {}
    first_exp: Dict[int, int] = {}
    last_clk: Dict[int, int] = {}
    first_clk: Dict[int, int] = {}
    last_any_clk: Optional[int] = None
    exp_cnt: Dict[int, int] = defaultdict(int)

    exp2clk: List[float] = []
    clk2cnv: List[float] = []
    fclk2cnv: List[float] = []
    exp2cnv: List[float] = []
    any2cnv: List[float] = []
    exp_before: List[float] = []
    cnv_wo_clk = 0
    cnv_n = 0

    # attribution window flags (minutes buckets)
    buckets = [5, 30, 60, 360, 1440, 10080]  # 5m..7d
    clk2cnv_bucket = {b: 0.0 for b in buckets}

    for iid, act, ts, _ in evs:
        if act == A_EXP:
            last_exp[iid] = ts
            first_exp.setdefault(iid, ts)
            exp_cnt[iid] += 1
        elif act == A_CLK:
            if iid in last_exp:
                exp2clk.append((ts - last_exp[iid]) / 60.0)
            last_clk[iid] = ts
            first_clk.setdefault(iid, ts)
            last_any_clk = ts
        elif act == A_CNV:
            cnv_n += 1
            exp_before.append(float(exp_cnt.get(iid, 0)))
            if iid in last_clk:
                dt = (ts - last_clk[iid]) / 60.0
                clk2cnv.append(dt)
                for b in buckets:
                    if dt <= b:
                        clk2cnv_bucket[b] += 1.0
            else:
                cnv_wo_clk += 1
            if iid in first_clk:
                fclk2cnv.append((ts - first_clk[iid]) / 60.0)
            if iid in first_exp:
                exp2cnv.append((ts - first_exp[iid]) / 60.0)
            if last_any_clk is not None:
                any2cnv.append((ts - last_any_clk) / 60.0)

    out: Dict[str, float] = {}
    out.update(stats(exp2clk, "attr_exp2clk_min"))
    out.update(stats(clk2cnv, "attr_clk2cnv_min"))
    out.update(stats(fclk2cnv, "attr_firstclk2cnv_min"))
    out.update(stats(exp2cnv, "attr_exp2cnv_min"))
    out.update(stats(any2cnv, "attr_anyclk2cnv_min"))
    out["attr_cnv_wo_prior_clk_cnt"] = float(cnv_wo_clk)
    out["attr_cnv_wo_prior_clk_rate"] = rate(cnv_wo_clk, cnv_n)
    out["attr_exp_before_cnv_mean"] = float(np.mean(exp_before)) if exp_before else 0.0
    for b, v in clk2cnv_bucket.items():
        out[f"attr_clk2cnv_within_{b}m_cnt"] = v
        out[f"attr_clk2cnv_within_{b}m_rate"] = rate(v, cnv_n)
    return out


def monetization(evs: Sequence[Event], price_p50: float) -> Dict[str, float]:
    convs = [e for e in evs if e[1] == A_CNV]
    clicks = [e for e in evs if e[1] == A_CLK]
    expos = [e for e in evs if e[1] == A_EXP]
    pay_cnt = float(len(convs))
    prices = [e[3] for e in convs if e[3] is not None]
    viewed = [e[3] for e in evs if e[3] is not None]
    arpu_sum = float(np.sum(prices)) if prices else 0.0
    arpu_mean = float(np.mean(prices)) if prices else 0.0
    arpu_p50 = float(np.median(prices)) if prices else 0.0
    cheap = sum(1 for p in prices if p <= price_p50)
    if len(prices) >= 2:
        med = float(np.median(prices))
        vs_self = rate(sum(1 for p in prices if p < med), len(prices))
    else:
        vs_self = 0.0
    view_mean = float(np.mean(viewed)) if viewed else 0.0
    # price histogram bins (relative to global p50)
    bins = [0.25, 0.5, 1.0, 1.5, 2.0, 4.0]
    hist = {b: 0.0 for b in bins}
    hist["gt_last"] = 0.0
    for p in prices:
        ratio = p / max(price_p50, 1.0)
        placed = False
        for b in bins:
            if ratio <= b:
                hist[b] += 1.0
                placed = True
                break
        if not placed:
            hist["gt_last"] += 1.0
    out = {
        "pay_cnt": pay_cnt,
        "pay_user": 1.0 if pay_cnt > 0 else 0.0,
        "arpu_sum_proxy": arpu_sum,
        "arpu_mean_proxy": arpu_mean,
        "arpu_p50_proxy": arpu_p50,
        "deal_cheap_cnv_share": rate(cheap, len(prices)) if prices else 0.0,
        "deal_vs_self_price_share": vs_self,
        "deal_price_gap_view_minus_pay": float(view_mean - arpu_mean),
        "cvr_pay_per_click": rate(pay_cnt, len(clicks)),
        "ctvr_pay_per_expose": rate(pay_cnt, len(expos)),
        "pay_price_known_rate": rate(len(prices), pay_cnt),
        "log1p_pay_cnt": float(math.log1p(pay_cnt)),
        "log1p_arpu_sum": float(math.log1p(arpu_sum)),
    }
    for b, v in hist.items():
        out[f"pay_price_bin_le_{b}_cnt"] = v
        out[f"pay_price_bin_le_{b}_share"] = rate(v, len(prices)) if prices else 0.0
    return out


def session_feats(evs: Sequence[Event]) -> Dict[str, float]:
    ss = sessions(evs)
    if not ss:
        return {f"sess_{k}": 0.0 for k in [
            "n", "len_mean", "len_p50", "len_max", "dur_min_mean", "dur_min_p50",
            "depth_clk_mean", "depth_cnv_mean", "bounce_rate", "last_len",
            "last_dur_min", "last_has_cnv", "last_has_clk", "cnv_sess_rate",
            "clk_sess_rate", "len_std", "gap_between_sess_min_mean",
        ]}
    lens = [len(s) for s in ss]
    durs = [(s[-1][2] - s[0][2]) / 60.0 for s in ss]
    clk_d = [sum(1 for e in s if e[1] == A_CLK) for s in ss]
    cnv_d = [sum(1 for e in s if e[1] == A_CNV) for s in ss]
    bounce = sum(1 for s in ss if len(s) <= 1)
    cnv_sess = sum(1 for d in cnv_d if d > 0)
    clk_sess = sum(1 for d in clk_d if d > 0)
    last = ss[-1]
    gaps = []
    for i in range(1, len(ss)):
        gaps.append((ss[i][0][2] - ss[i - 1][-1][2]) / 60.0)
    return {
        "sess_n": float(len(ss)),
        "sess_len_mean": float(np.mean(lens)),
        "sess_len_p50": float(np.median(lens)),
        "sess_len_max": float(np.max(lens)),
        "sess_len_std": float(np.std(lens)),
        "sess_dur_min_mean": float(np.mean(durs)),
        "sess_dur_min_p50": float(np.median(durs)),
        "sess_depth_clk_mean": float(np.mean(clk_d)),
        "sess_depth_cnv_mean": float(np.mean(cnv_d)),
        "sess_bounce_rate": rate(bounce, len(ss)),
        "sess_last_len": float(len(last)),
        "sess_last_dur_min": float((last[-1][2] - last[0][2]) / 60.0),
        "sess_last_has_cnv": 1.0 if any(e[1] == A_CNV for e in last) else 0.0,
        "sess_last_has_clk": 1.0 if any(e[1] == A_CLK for e in last) else 0.0,
        "sess_cnv_sess_rate": rate(cnv_sess, len(ss)),
        "sess_clk_sess_rate": rate(clk_sess, len(ss)),
        "sess_gap_between_sess_min_mean": float(np.mean(gaps)) if gaps else -1.0,
    }


def diversity(evs: Sequence[Event], t_end: int) -> Dict[str, float]:
    if not evs:
        return {"hist_len": 0.0, "active_days": 0.0}
    days = {e[2] // 86400 for e in evs}
    clk_items: Dict[int, int] = defaultdict(int)
    exp_cnt: Dict[int, int] = defaultdict(int)
    clk_cnt: Dict[int, int] = defaultdict(int)
    cnv_cnt: Dict[int, int] = defaultdict(int)
    for iid, act, _, _ in evs:
        if act == A_CLK:
            clk_items[iid] += 1
            clk_cnt[iid] += 1
        elif act == A_EXP:
            exp_cnt[iid] += 1
        elif act == A_CNV:
            cnv_cnt[iid] += 1
    tot = sum(clk_items.values())
    top1 = max(clk_items.values()) / tot if tot else 0.0
    top3 = sum(sorted(clk_items.values(), reverse=True)[:3]) / tot if tot else 0.0

    def since(act: int) -> float:
        ts = [e[2] for e in evs if e[1] == act]
        return float((t_end - max(ts)) / 3600.0) if ts else -1.0

    gaps = [evs[i][2] - evs[i - 1][2] for i in range(1, len(evs))]
    # hour-of-day / weekday histograms on clicks
    hod = [0.0] * 6  # 4h bins
    dow = [0.0] * 7
    for _, act, ts, _ in evs:
        if act != A_CLK:
            continue
        # rough UTC buckets
        hod[(ts % 86400) // 14400] += 1.0
        dow[((ts // 86400) + 4) % 7] += 1.0  # epoch offset-ish
    n_clk = sum(hod) or 1.0
    out = {
        "hist_len": float(len(evs)),
        "span_day": float((evs[-1][2] - evs[0][2]) / 86400.0),
        "active_days": float(len(days)),
        "item_entropy_clk": entropy(clk_items),
        "item_entropy_cnv": entropy(cnv_cnt),
        "top1_item_clk_share": float(top1),
        "top3_item_clk_share": float(top3),
        "hours_since_last_exp": since(A_EXP),
        "hours_since_last_clk": since(A_CLK),
        "hours_since_last_cnv": since(A_CNV),
        "reexpose_rate": rate(sum(1 for c in exp_cnt.values() if c >= 2), len(exp_cnt)),
        "repeat_clk_rate": rate(sum(1 for c in clk_cnt.values() if c >= 2), len(clk_cnt)),
        "repeat_cnv_rate": rate(sum(1 for c in cnv_cnt.values() if c >= 2), len(cnv_cnt)),
        "gap_sec_p50": float(np.median(gaps)) if gaps else -1.0,
        "gap_sec_mean": float(np.mean(gaps)) if gaps else -1.0,
        "gap_sec_std": float(np.std(gaps)) if gaps else -1.0,
        "n_uniq_clk_item": float(len(clk_cnt)),
        "n_uniq_cnv_item": float(len(cnv_cnt)),
        "n_uniq_exp_item": float(len(exp_cnt)),
    }
    for i, v in enumerate(hod):
        out[f"clk_hod_bin{i}_share"] = v / n_clk
    for i, v in enumerate(dow):
        out[f"clk_dow_{i}_share"] = v / n_clk
    return out


def ui_summary(evs: Sequence[Event]) -> Dict[str, float]:
    by: Dict[int, List[Event]] = defaultdict(list)
    for e in evs:
        by[e[0]].append(e)
    pair_ctr, pair_cvr, pair_exp, clk2cnv = [], [], [], []
    multitouch = 0
    only_exp = 0
    for lst in by.values():
        n_exp = sum(1 for e in lst if e[1] == A_EXP)
        n_clk = sum(1 for e in lst if e[1] == A_CLK)
        n_cnv = sum(1 for e in lst if e[1] == A_CNV)
        pair_exp.append(n_exp)
        pair_ctr.append(rate(n_clk, n_exp))
        pair_cvr.append(rate(n_cnv, n_clk))
        if n_exp >= 2 and n_cnv >= 1:
            multitouch += 1
        if n_exp > 0 and n_clk == 0 and n_cnv == 0:
            only_exp += 1
        last_c = None
        for e in lst:
            if e[1] == A_CLK:
                last_c = e[2]
            elif e[1] == A_CNV and last_c is not None:
                clk2cnv.append((e[2] - last_c) / 60.0)
    return {
        "ui_n_items_touched": float(len(by)),
        "ui_exp_per_item_mean": float(np.mean(pair_exp)) if pair_exp else 0.0,
        "ui_exp_per_item_p90": float(np.percentile(pair_exp, 90)) if pair_exp else 0.0,
        "ui_pair_ctr_mean": float(np.mean(pair_ctr)) if pair_ctr else 0.0,
        "ui_pair_cvr_mean": float(np.mean(pair_cvr)) if pair_cvr else 0.0,
        "ui_multitouch_cnv_items": float(multitouch),
        "ui_only_exp_items": float(only_exp),
        "ui_only_exp_share": rate(only_exp, len(by)),
        "ui_pair_clk2cnv_min_p50": float(np.median(clk2cnv)) if clk2cnv else -1.0,
        "ui_pair_clk2cnv_min_mean": float(np.mean(clk2cnv)) if clk2cnv else -1.0,
    }


def transitions(evs: Sequence[Event]) -> Dict[str, float]:
    # order-1
    t1: Dict[Tuple[int, int], int] = defaultdict(int)
    for a, b in zip(evs, evs[1:]):
        t1[(a[1], b[1])] += 1
    tot1 = sum(t1.values()) + EPS
    out = {}
    for a in range(3):
        for b in range(3):
            out[f"trans_{A_NAME[a]}_to_{A_NAME[b]}"] = t1[(a, b)] / tot1
    # order-2 (selected: *→clk→cnv, exp→clk→*, etc.)
    t2: Dict[Tuple[int, int, int], int] = defaultdict(int)
    for a, b, c in zip(evs, evs[1:], evs[2:]):
        t2[(a[1], b[1], c[1])] += 1
    tot2 = sum(t2.values()) + EPS
    keys = [
        (A_EXP, A_CLK, A_CNV),
        (A_EXP, A_EXP, A_CLK),
        (A_CLK, A_CLK, A_CNV),
        (A_CLK, A_EXP, A_CLK),
        (A_EXP, A_CLK, A_CLK),
        (A_CNV, A_EXP, A_CLK),
        (A_CLK, A_CNV, A_EXP),
        (A_EXP, A_CNV, A_EXP),
    ]
    for a, b, c in keys:
        out[f"trans2_{A_NAME[a]}_{A_NAME[b]}_{A_NAME[c]}"] = t2[(a, b, c)] / tot2
    return out


def cross_interactions(feats: Dict[str, float]) -> Dict[str, float]:
    """Cheap multiplicative crosses to push toward 1000 dims."""
    keys = [
        "life_ctr", "life_cvr", "life_ctcvr",
        "7d_ctr", "7d_cvr", "1d_ctr", "1d_cvr",
        "pay_cnt", "arpu_sum_proxy", "arpu_mean_proxy",
        "deal_cheap_cnv_share", "sess_bounce_rate", "sess_n",
        "active_days", "hist_len", "item_entropy_clk",
        "hours_since_last_clk", "hours_since_last_cnv",
        "attr_clk2cnv_min_p50", "attr_anyclk2cnv_min_p50",
        "dec_hl1d_dec_ctr", "dec_hl1d_dec_cvr",
    ]
    present = [k for k in keys if k in feats]
    out: Dict[str, float] = {}
    for i, a in enumerate(present):
        va = feats[a]
        out[f"sq_{a}"] = float(va * va)
        out[f"log1p_abs_{a}"] = float(math.log1p(abs(va)))
        for b in present[i + 1 :]:
            out[f"x_{a}__{b}"] = float(va * feats[b])
    return out


def split_prefix_suffix(
    evs: Sequence[Event], *, prefix_frac: float = 0.8
) -> Tuple[List[Event], List[Event]]:
    """Causal split: features from prefix, label from suffix (future)."""
    if len(evs) < 5:
        return list(evs), []
    cut = max(1, min(len(evs) - 1, int(len(evs) * prefix_frac)))
    # also enforce time order cut at quantile timestamp
    t0, t1 = evs[0][2], evs[-1][2]
    t_cut = t0 + int((t1 - t0) * prefix_frac)
    prefix = [e for e in evs if e[2] <= t_cut]
    suffix = [e for e in evs if e[2] > t_cut]
    if len(prefix) < 1:
        prefix, suffix = list(evs[:cut]), list(evs[cut:])
    return prefix, suffix


def build_user_features(evs: Sequence[Event], price_p50: float) -> Dict[str, float]:
    if not evs:
        return {"hist_len": 0.0, "active_days": 0.0, "pay_user": 0.0, "_seq_t_end": 0.0}
    t_end = evs[-1][2]
    feats: Dict[str, float] = {"_seq_t_end": float(t_end)}

    life = funnel(evs, t_end, 10**12)
    for k, v in life.items():
        feats[f"life_{k}"] = v
    for wname, wsec in WINS.items():
        for k, v in funnel(evs, t_end, wsec).items():
            feats[f"{wname}_{k}"] = v

    # trends
    for short, long in [("1d", "7d"), ("7d", "30d"), ("1h", "1d"), ("3d", "14d")]:
        for m in ("ctr", "cvr", "ctcvr", "n_clk", "n_cnv"):
            feats[f"trend_{m}_{short}_minus_{long}"] = feats[f"{short}_{m}"] - feats[f"{long}_{m}"]
            den = feats[f"{long}_{m}"]
            feats[f"trend_{m}_{short}_over_{long}"] = rate(feats[f"{short}_{m}"], den if den != 0 else EPS)

    for dname, hl in DECAYS.items():
        for k, v in decay_counts(evs, t_end, hl).items():
            feats[f"dec_{dname}_{k}"] = v

    feats.update(session_feats(evs))
    feats.update(diversity(evs, t_end))
    feats.update(attribution(evs))
    mon = monetization(evs, price_p50)
    ad = feats.get("active_days", 0.0)
    mon["arpu_per_active_day_proxy"] = mon["arpu_sum_proxy"] / ad if ad > 0 else 0.0
    mon["deal_hesitation_exp_before_cnv"] = feats.get("attr_exp_before_cnv_mean", 0.0)
    feats.update(mon)
    feats.update(ui_summary(evs))
    feats.update(transitions(evs))
    feats.update(cross_interactions(feats))
    return feats


def enrich_user_table(feats: Dict[str, float], urec: Optional[pd.Series]) -> Dict[str, float]:
    if urec is None:
        return feats
    for col in ["103", "104", "105", "109"]:
        if col not in urec.index:
            continue
        v = urec[col]
        try:
            fv = float(v) if pd.notna(v) else -1.0
            if not math.isfinite(fv):
                fv = -1.0
        except Exception:
            fv = -1.0
        feats[f"u_side_{col}"] = fv
        feats[f"u_side_{col}_known"] = 0.0 if fv < 0 else 1.0
    for col in ["106", "107", "108", "110"]:
        if col not in urec.index:
            continue
        v = urec[col]
        if v is None or (isinstance(v, float) and math.isnan(v)):
            feats[f"u_side_{col}_len"] = 0.0
            feats[f"u_side_{col}_ent"] = 0.0
            continue
        try:
            lst = list(v)
        except Exception:
            lst = []
        feats[f"u_side_{col}_len"] = float(len(lst))
        cnt: Dict[Any, int] = defaultdict(int)
        for x in lst:
            cnt[x] += 1
        feats[f"u_side_{col}_ent"] = entropy(cnt)
    return feats


def compare_feature_lists(selected: Sequence[str], prev_json: Path) -> Dict[str, Any]:
    """Overlap vs the previous 128-feat F-score board."""
    if not prev_json.is_file():
        return {"prev_path": str(prev_json), "n_prev": 0, "n_new": len(selected)}
    prev = json.loads(prev_json.read_text(encoding="utf-8"))
    old_names = [x["name"] if isinstance(x, dict) else str(x) for x in prev]
    new_s, old_s = set(selected), set(old_names)
    return {
        "prev_path": str(prev_json),
        "n_prev": len(old_names),
        "n_new": len(selected),
        "n_overlap": len(new_s & old_s),
        "n_only_new": len(new_s - old_s),
        "n_only_prev": len(old_s - new_s),
        "only_new": [n for n in selected if n not in old_s],
        "only_prev": [n for n in old_names if n not in new_s],
        "overlap": [n for n in selected if n in old_s],
    }


def pad_or_trim(df: pd.DataFrame, target_dim: int, exclude: Sequence[str]) -> pd.DataFrame:
    """Ensure roughly target_dim feature columns (pad zeros / trim low-var)."""
    feat_cols = [c for c in df.columns if c not in exclude]
    if len(feat_cols) > target_dim:
        # keep highest variance
        var = df[feat_cols].var(axis=0).sort_values(ascending=False)
        keep = list(var.index[:target_dim])
        return df[list(exclude) + keep]
    if len(feat_cols) < target_dim:
        n_pad = target_dim - len(feat_cols)
        pad = pd.DataFrame(
            0.0,
            index=df.index,
            columns=[f"pad_zero_{i}" for i in range(n_pad)],
        )
        df = pd.concat([df, pad], axis=1)
    return df


def select_and_train(
    X: np.ndarray,
    y: np.ndarray,
    names: List[str],
    select_k: int,
    method: str,
) -> Dict[str, Any]:
    # clean
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    # variance filter
    vt = VarianceThreshold(threshold=1e-8)
    Xv = vt.fit_transform(X)
    kept = [n for n, m in zip(names, vt.get_support()) if m]
    print(f"  after variance filter: {len(kept)}")

    k = min(select_k, len(kept), max(int(Xv.shape[0] * 0.5), 8))
    if method == "mi":
        selector = SelectKBest(mutual_info_classif, k=k)
    else:
        selector = SelectKBest(f_classif, k=k)
    Xs = selector.fit_transform(Xv, y)
    scores = selector.scores_
    support = selector.get_support()
    selected = [n for n, s in zip(kept, support) if s]
    score_map = {n: float(sc) for n, sc, s in zip(kept, scores, support) if s and np.isfinite(sc)}
    ranked = sorted(score_map.items(), key=lambda kv: -kv[1])

    Xtr, Xte, ytr, yte = train_test_split(Xs, y, test_size=0.25, random_state=42, stratify=y)

    models = {
        "hgb": HistGradientBoostingClassifier(max_depth=4, max_iter=120, learning_rate=0.08, random_state=42),
        "logreg": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=400, C=0.5, solver="lbfgs")),
        ]),
    }
    results = {}
    for name, model in models.items():
        t0 = time.time()
        model.fit(Xtr, ytr)
        proba = model.predict_proba(Xte)[:, 1]
        pred = (proba >= 0.5).astype(int)
        results[name] = {
            "auc": float(roc_auc_score(yte, proba)),
            "ap": float(average_precision_score(yte, proba)),
            "acc": float(accuracy_score(yte, pred)),
            "sec": float(time.time() - t0),
            "n_selected": len(selected),
        }
        print(f"  {name}: AUC={results[name]['auc']:.4f} AP={results[name]['ap']:.4f} Acc={results[name]['acc']:.4f}")

    return {
        "selected_features": [{"name": n, "score": score_map.get(n, 0.0)} for n, _ in ranked],
        "metrics": results,
        "n_features_in": len(names),
        "n_after_variance": len(kept),
        "n_selected": len(selected),
        "label_pos_rate": float(y.mean()),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/tencent_subset"))
    ap.add_argument("--max-users", type=int, default=8000)
    ap.add_argument("--target-dim", type=int, default=1000)
    ap.add_argument("--select-k", type=int, default=128)
    ap.add_argument("--select-method", choices=["f", "mi"], default="f")
    ap.add_argument(
        "--label",
        choices=["future_cnv", "pay_user", "has_cnv_7d"],
        default="future_cnv",
        help="future_cnv = causal: prefix feats → suffix has conversion",
    )
    ap.add_argument("--prefix-frac", type=float, default=0.75)
    ap.add_argument("--out", type=Path, default=Path("results/tencent_gr_fs"))
    args = ap.parse_args()

    seq_paths = parquets(args.root / "seq")
    if not seq_paths:
        raise SystemExit(f"no seq under {args.root}/seq")

    print("loading price proxy ...")
    price = load_price(args.root / "item_feat")
    pvals = np.asarray(list(price.values()), dtype=float) if price else np.asarray([])
    g_p50 = float(np.median(pvals)) if pvals.size else 0.0
    print(f"  priced={len(price)} p50={g_p50:.1f}")

    rows: List[Dict[str, float]] = []
    uids: List[int] = []
    print(f"featurizing <= {args.max_users} users → ~{args.target_dim} dims ...")
    for uid, seq in iter_users(seq_paths, args.max_users):
        evs = parse_seq(seq, price)
        if args.label == "future_cnv":
            prefix, suffix = split_prefix_suffix(evs, prefix_frac=args.prefix_frac)
            feats = build_user_features(prefix, g_p50)
            feats["_label_future_cnv"] = 1.0 if any(e[1] == A_CNV for e in suffix) else 0.0
            feats["_label_pay_user"] = 1.0 if any(e[1] == A_CNV for e in evs) else 0.0
        else:
            feats = build_user_features(evs, g_p50)
            feats["_label_pay_user"] = feats.get("pay_user", 0.0)
            feats["_label_has_cnv_7d"] = 1.0 if feats.get("7d_n_cnv", 0.0) > 0 else 0.0
            feats["_label_future_cnv"] = feats["_label_pay_user"]
        feats["user_id"] = float(uid)
        rows.append(feats)
        uids.append(uid)
        if len(rows) % 2000 == 0:
            print(f"  {len(rows)} ...")

    df = pd.DataFrame(rows)
    print("joining user_feat ...")
    udf = load_users(args.root / "user_feat", set(uids)).set_index("user_id")
    enriched = []
    for _, r in df.iterrows():
        uid = int(r["user_id"])
        urec = udf.loc[uid] if uid in udf.index else None
        d = {k: float(v) for k, v in r.items()}
        d = enrich_user_table(d, urec)
        enriched.append(d)
    df = pd.DataFrame(enriched)

    label_map = {
        "future_cnv": "_label_future_cnv",
        "pay_user": "_label_pay_user",
        "has_cnv_7d": "_label_has_cnv_7d",
    }
    label_col = label_map[args.label]

    exclude = [
        "user_id",
        "_label_pay_user",
        "_label_has_cnv_7d",
        "_label_future_cnv",
        "_seq_t_end",
    ]
    df = pad_or_trim(df, args.target_dim, exclude)
    feat_cols = [c for c in df.columns if c not in exclude]
    print(f"feature_dim={len(feat_cols)}")

    y = df[label_col].to_numpy().astype(int)
    # soft leak filter still useful for non-causal labels
    leak_substrings = []
    if args.label != "future_cnv":
        leak_substrings = ["pay_cnt", "pay_user", "n_cnv", "cnv_share", "ctcvr", "arpu_sum"]
    feat_cols_use = [
        c for c in feat_cols
        if not any(s in c for s in leak_substrings)
    ]
    print(f"after leak drop: {len(feat_cols_use)} (dropped {len(feat_cols)-len(feat_cols_use)})")

    X = df[feat_cols_use].to_numpy(dtype=float)
    report = select_and_train(X, y, feat_cols_use, args.select_k, args.select_method)

    args.out.mkdir(parents=True, exist_ok=True)
    selected_names = [x["name"] for x in report["selected_features"]]
    keep_extra = [
        c
        for c in ("user_id", "_seq_t_end", "_label_future_cnv", "_label_pay_user", "life_ctcvr")
        if c in df.columns
    ]
    sel_cols = []
    seen = set()
    for c in keep_extra + selected_names:
        if c in df.columns and c not in seen:
            seen.add(c)
            sel_cols.append(c)
    df[sel_cols].to_parquet(args.out / "user_feats_selected.parquet", index=False)
    df[["user_id"] + feat_cols_use[: min(200, len(feat_cols_use))]].to_parquet(
        args.out / "user_feats_preview.parquet", index=False
    )
    (args.out / "selected_features.json").write_text(
        json.dumps(report["selected_features"], indent=2), encoding="utf-8"
    )
    (args.out / "FEATURES_150.txt").write_text(
        "\n".join(f"{i:3d}  {n}" for i, n in enumerate(selected_names, 1)) + "\n",
        encoding="utf-8",
    )
    prev_path = Path("results/tencent_gr_fs/selected_features.json")
    compare = compare_feature_lists(selected_names, prev_path)
    (args.out / "FEATURE_COMPARE_150_vs_128.json").write_text(
        json.dumps(compare, indent=2), encoding="utf-8"
    )
    summary = {
        "n_users": int(len(df)),
        "n_features_generated": len(feat_cols),
        "n_features_used": len(feat_cols_use),
        "label": args.label,
        "select_k": args.select_k,
        "select_method": args.select_method,
        "metrics": report["metrics"],
        "label_pos_rate": report["label_pos_rate"],
        "top20": report["selected_features"][:20],
        "compare_vs_128": {
            "n_overlap": compare.get("n_overlap"),
            "n_only_new": compare.get("n_only_new"),
            "n_only_prev": compare.get("n_only_prev"),
        },
    }
    (args.out / "train_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    overlap_n = compare.get("n_overlap", 0)
    only_new = compare.get("only_new") or []
    only_prev = compare.get("only_prev") or []
    feat_md = "\n".join(f"{i}. `{n}`" for i, n in enumerate(selected_names, 1))
    md = f"""# TencentGR auto-1000 → feature-select → train

- users: **{len(df)}**
- generated dims: **{len(feat_cols)}** (target {args.target_dim})
- used after leak-drop: **{len(feat_cols_use)}**
- selected: **{report['n_selected']}** via `{args.select_method}`
- label: `{args.label}` (pos rate={report['label_pos_rate']:.3f})
- vs previous 128: overlap **{overlap_n}**, only-new **{len(only_new)}**, dropped **{len(only_prev)}**

## Metrics

```json
{json.dumps(report['metrics'], indent=2)}
```

## Top-20 selected

```json
{json.dumps(report['selected_features'][:20], indent=2)}
```

## 150-feature list

{feat_md}

## vs 128 board

- overlap: {overlap_n}
- new (not in 128): {", ".join(f"`{n}`" for n in only_new) or "(none)"}
- dropped from 128: {", ".join(f"`{n}`" for n in only_prev) or "(none)"}

## Recipe

1. **Generate**: combinatorial windows × funnel × decay × session × attribution × ARPU/deal × Markov × TOD × crosses ≈ 1000.
2. **Select**: VarianceThreshold → SelectKBest(F / MI) → top-k.
3. **Train**: HistGradientBoosting + LogisticRegression holdout AUC/AP.

Leakage note: for `pay_user`, raw `pay_cnt` / `life_n_cnv` / `arpu_sum` are dropped before selection.
Label `future_cnv` uses prefix features and suffix conversion, so conversion-count features are not a direct leak.
"""
    (args.out / "FS_TRAIN_REPORT.md").write_text(md, encoding="utf-8")
    print(json.dumps(summary, indent=2)[:2500])
    print("wrote", args.out)


if __name__ == "__main__":
    main()
