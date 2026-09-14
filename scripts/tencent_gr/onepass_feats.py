#!/usr/bin/env python3
"""One-pass TencentGR user features — read the loop, that's the business.

Event: (item_id, act, ts, price)
  act 0=曝光  1=点击  2=转化
  一条用户序列按时间走一遍：漏斗 / 衰减 / 场次 / 归因 / 转移 同时积。
  过完再算比率和 cross（乘，不再扫序列）。

  python3 scripts/tencent_gr/onepass_feats.py
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

EXP, CLK, CNV = 0, 1, 2
NAME = {0: "exp", 1: "clk", 2: "cnv"}
# 业务窗：近 3/7/14/30 天还在不在逛、点、买
WINS = {"3d": 3 * 86400, "7d": 7 * 86400, "14d": 14 * 86400, "30d": 30 * 86400}
# 半衰期：越近的点击/转化权重越大（1天掉一半 / 3天 / 7天）
HLS = {"hl1d": 86400, "hl3d": 3 * 86400, "hl7d": 7 * 86400}
SESS_GAP = 30 * 60  # 30min 没动 = 这场逛完了
EPS = 1e-6
Event = Tuple[int, int, int, Optional[float]]


def rate(n: float, d: float) -> float:
    return float(n) / float(d) if d > 0 else 0.0


def entropy(cnt: Dict[int, int]) -> float:
    """点得越分散（什么都点）熵越大；只点一两件越小。"""
    tot = float(sum(cnt.values()))
    if tot <= 0:
        return 0.0
    h = 0.0
    for c in cnt.values():
        if c <= 0:
            continue
        p = c / tot
        h -= p * math.log(p + EPS)
    return h


def onepass(evs: List[Event]) -> Dict[str, float]:
    if not evs:
        return {"hist_len": 0.0}
    evs = sorted(evs, key=lambda e: e[2])
    t_end = evs[-1][2]

    # —— 漏斗：各窗里曝光/点击/转化条数（谁还在买）
    n = {w: {EXP: 0, CLK: 0, CNV: 0, "tot": 0} for w in ["life", *WINS]}
    # —— 衰减：同样是点/买，但近的更重（活跃度，不是生硬截窗）
    dec = {h: {EXP: 0.0, CLK: 0.0, CNV: 0.0} for h in HLS}
    lam = {h: math.log(2.0) / sec for h, sec in HLS.items()}
    # —— 场次：一场连续逛；蹦出 = 看一眼就走
    sess_n = bounce = clk_sess = cnv_sess = 0
    depth_clk: List[int] = []
    depth_cnv: List[int] = []
    cur_len = cur_clk = cur_cnv = 0
    prev_ts: Optional[int] = None

    def close_sess() -> None:
        nonlocal sess_n, bounce, clk_sess, cnv_sess, cur_len, cur_clk, cur_cnv
        if cur_len <= 0:
            return
        sess_n += 1
        if cur_len <= 1:
            bounce += 1  # 看一条就走
        if cur_clk:
            clk_sess += 1  # 这场有没有动手点
        if cur_cnv:
            cnv_sess += 1  # 这场有没有成交
        depth_clk.append(cur_clk)
        depth_cnv.append(cur_cnv)
        cur_len = cur_clk = cur_cnv = 0

    # —— 归因：上次任意点击 → 这笔成交隔了多久（决策时长；空=-1）
    last_any_clk: Optional[int] = None
    last_clk_item: Dict[int, int] = {}
    any2cnv: List[float] = []
    item2cnv: List[float] = []
    # —— 转移：上一步 act → 这一步（逛完点 / 买完又逛）
    t1 = defaultdict(int)
    t2 = defaultdict(int)
    prev_act = prev2_act = None
    # —— 多样性：点过哪些商品
    clk_iid: Dict[int, int] = defaultdict(int)
    days = set()

    for iid, act, ts, _price in evs:
        age = max(t_end - ts, 0)
        days.add(ts // 86400)

        # 1. 漏斗窗：这条事件算进哪些「近 X 天」
        n["life"][act] += 1
        n["life"]["tot"] += 1
        for w, sec in WINS.items():
            if age < sec:  # (t_end-W, t_end]
                n[w][act] += 1
                n[w]["tot"] += 1

        # 2. 衰减：同样一条点击，昨天的比上周的值钱
        for h, l in lam.items():
            dec[h][act] += math.exp(-l * age)

        # 3. 场次：超过 30min 没动作 = 上一场结束
        if prev_ts is not None and ts - prev_ts > SESS_GAP:
            close_sess()
        cur_len += 1
        cur_clk += int(act == CLK)
        cur_cnv += int(act == CNV)
        prev_ts = ts

        # 4. 归因：点击记下时间；转化时回头找上次点
        #    买后点击的对偶（pending ← CNV，CLK 闭上）见 onepass_post.py
        if act == CLK:
            last_any_clk = ts
            last_clk_item[iid] = ts
            clk_iid[iid] += 1
        elif act == CNV:
            if last_any_clk is not None:
                any2cnv.append((ts - last_any_clk) / 60.0)  # 全局：逛完多久下单
            if iid in last_clk_item:
                item2cnv.append((ts - last_clk_item[iid]) / 60.0)  # 同商品：看这件多久下单

        # 5. Markov：上一条 act 到这一条
        if prev_act is not None:
            t1[(prev_act, act)] += 1
        if prev2_act is not None:
            t2[(prev2_act, prev_act, act)] += 1
        prev2_act, prev_act = prev_act, act

    close_sess()

    def funnel(w: str) -> Dict[str, float]:
        e, c, v, tot = n[w][EXP], n[w][CLK], n[w][CNV], n[w]["tot"]
        return {
            "n_exp": float(e),
            "n_clk": float(c),
            "n_cnv": float(v),
            "ctr": rate(c, e),    # 看见会不会点
            "cvr": rate(v, c),    # 点了会不会买
            "ctcvr": rate(v, e),  # 看见会不会买
            "clk_share": rate(c, tot),
            "cnv_share": rate(v, tot),
        }

    life, d3, d7, d14, d30 = map(funnel, ["life", "3d", "7d", "14d", "30d"])
    pay_cnt = life["n_cnv"]
    hist_len = float(len(evs))
    active_days = float(len(days))
    h_clk = entropy(clk_iid)
    bounce_rate = rate(bounce, sess_n)
    any_p50 = float(np.median(any2cnv)) if any2cnv else -1.0
    any_std = float(np.std(any2cnv)) if any2cnv else -1.0
    item_p50 = float(np.median(item2cnv)) if item2cnv else -1.0
    pairs = sum(t1.values()) + EPS
    triples = sum(t2.values()) + EPS

    feats: Dict[str, float] = {
        # 漏斗：买在行为里占多大（不是绝对次数，避免重度用户刷榜）
        "life_cnv_share": life["cnv_share"],
        "life_ctr": life["ctr"],
        "life_cvr": life["cvr"],
        "life_ctcvr": life["ctcvr"],
        "pay_cnt": pay_cnt,
        "3d_clk_share": d3["clk_share"],
        "7d_ctr": d7["ctr"],
        "7d_cvr": d7["cvr"],
        "14d_clk_share": d14["clk_share"],
        "30d_clk_share": d30["clk_share"],
        "30d_cnv_share": d30["cnv_share"],
        "3d_n_cnv": d3["n_cnv"],
        "14d_n_cnv": d14["n_cnv"],
        "3d_cvr": d3["cvr"],
        "14d_cvr": d14["cvr"],
        "trend_n_cnv_3d_minus_14d": d3["n_cnv"] - d14["n_cnv"],  # 近3天买变少会 <0
        "trend_cvr_3d_minus_14d": d3["cvr"] - d14["cvr"],
        # 衰减：还在不在点/买（半衰期，比硬截 7d 更滑）
        "dec_hl1d_dec_clk": dec["hl1d"][CLK],
        "dec_hl3d_dec_clk": dec["hl3d"][CLK],
        "dec_hl7d_dec_clk": dec["hl7d"][CLK],
        "dec_hl3d_dec_cnv": dec["hl3d"][CNV],
        "dec_hl7d_dec_cnv": dec["hl7d"][CNV],
        # 场次：逛得碎不碎、一场里买不买
        "sess_n": float(sess_n),
        "sess_bounce_rate": bounce_rate,
        "sess_depth_cnv_mean": float(np.mean(depth_cnv)) if depth_cnv else 0.0,
        "sess_cnv_sess_rate": rate(cnv_sess, sess_n),
        # 归因：决策时长散不散（标准差大 = 有时秒下有时纠结很久）
        "attr_anyclk2cnv_min_p50": any_p50,
        "attr_anyclk2cnv_min_std": any_std,
        "attr_clk2cnv_min_p50": item_p50,
        "hist_len": hist_len,
        "active_days": active_days,
        "item_entropy_clk": h_clk,
        # 转移：买完是接着逛还是走了；曝光后点不点
        "trans_cnv_to_exp": t1[(CNV, EXP)] / pairs,
        "trans_exp_to_exp": t1[(EXP, EXP)] / pairs,
        "trans_exp_to_clk": t1[(EXP, CLK)] / pairs,
        "trans2_exp_cnv_exp": t2[(EXP, CNV, EXP)] / triples,  # 看见→买→又看
    }

    # cross：强度 × 犹豫 / 强度 × 广度（不再扫序列）
    feats["x_life_ctcvr__sess_bounce_rate"] = life["ctcvr"] * bounce_rate
    feats["x_life_cvr__sess_bounce_rate"] = life["cvr"] * bounce_rate
    feats["x_life_ctr__active_days"] = life["ctr"] * active_days
    feats["x_life_ctcvr__active_days"] = life["ctcvr"] * active_days
    feats["x_life_ctr__sess_n"] = life["ctr"] * sess_n
    feats["x_life_ctcvr__sess_n"] = life["ctcvr"] * sess_n
    feats["x_pay_cnt__active_days"] = pay_cnt * active_days
    feats["x_life_cvr__7d_ctr"] = life["cvr"] * d7["ctr"]
    feats["x_life_ctr__attr_anyclk2cnv_min_p50"] = life["ctr"] * any_p50
    feats["x_sess_n__item_entropy_clk"] = sess_n * h_clk
    feats["x_hist_len__item_entropy_clk"] = hist_len * h_clk
    feats["x_life_ctcvr__pay_cnt"] = life["ctcvr"] * pay_cnt
    feats["log1p_abs_life_ctcvr"] = math.log1p(abs(life["ctcvr"]))
    return feats


def _demo() -> None:
    # 一个用户：看→点→隔很久再买；第二天又看一眼就走（bounce）
    t0 = 1_700_000_000
    evs: List[Event] = [
        (11, EXP, t0, 99.0),
        (11, CLK, t0 + 20, 99.0),
        (11, CNV, t0 + 3600, 99.0),          # 点完 60min 才买
        (22, EXP, t0 + 86400 + 10, 50.0),    # 次日另一件，看一眼
    ]
    f = onepass(evs)
    want = [
        "life_cnv_share",
        "sess_bounce_rate",
        "sess_n",
        "attr_anyclk2cnv_min_p50",
        "trans_exp_to_clk",
        "x_life_ctcvr__sess_bounce_rate",
    ]
    for k in want:
        print(f"{k:36s} {f[k]:.4f}")


if __name__ == "__main__":
    _demo()
