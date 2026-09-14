#!/usr/bin/env python3
"""行为序列 → 软归因 prototype。没有 GT。

一笔转化 = query（item, price, cnv_ts）。
history = cnv_ts 之前的曝光/点击（最后 K 条）。
规则 last-touch / 时间衰减 / 显式权重 softmax：把 GMV=price 分到触点上。

业务：
  同品权重大  → SKU 漏斗在分账
  任意点击权重大 → 人在场在分账
  这批数据同品触点极稀，last-item 会经常分不到，last-any 才有数。

  python3 scripts/tencent_gr/seq_attr_proto.py
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

EXP, CLK, CNV = 0, 1, 2
NAME = {0: "exp", 1: "clk", 2: "cnv"}
Event = Tuple[int, int, int, Optional[float]]
K = 32
HL_SEC = 3 * 86400.0  # 衰减半衰期 3d


def _logp(p: Optional[float]) -> float:
    if p is None or not math.isfinite(p) or p < 0:
        return 0.0
    return math.log1p(p)


def history_before(evs: List[Event], cnv_ts: int, cnv_item: int) -> List[Dict]:
    """严格早于转化的曝光/点击，截最后 K 条。"""
    evs = sorted(evs, key=lambda e: (e[2], e[1]))
    rows = []
    for iid, act, ts, price in evs:
        if ts >= cnv_ts:
            break
        if act not in (EXP, CLK):
            continue
        rows.append(
            {
                "item_id": iid,
                "act": act,
                "ts": ts,
                "price": price,
                "same": int(iid == cnv_item),
                "dt": cnv_ts - ts,
                "log_price": _logp(price),
            }
        )
    return rows[-K:]


def rule_last(hist: List[Dict], *, same: bool, act: Optional[int]) -> Optional[int]:
    """从后往前第一条满足 same/act 的下标。"""
    for i in range(len(hist) - 1, -1, -1):
        h = hist[i]
        if act is not None and h["act"] != act:
            continue
        if same and not h["same"]:
            continue
        return i
    return None


def one_hot(n: int, i: Optional[int]) -> np.ndarray:
    a = np.zeros(n, dtype=float)
    if i is not None and n:
        a[i] = 1.0
    return a


def time_decay(hist: List[Dict], hl: float = HL_SEC) -> np.ndarray:
    if not hist:
        return np.zeros(0)
    lam = math.log(2.0) / hl
    w = np.array([math.exp(-lam * h["dt"]) for h in hist], dtype=float)
    s = w.sum()
    return w / s if s > 0 else w


def soft_attr(
    hist: List[Dict],
    *,
    w_clk: float = 1.2,
    w_same: float = 2.0,
    w_dt: float = -0.35,
    w_price: float = 0.1,
    q_price: float = 0.0,
) -> np.ndarray:
    """score = w_clk*1{clk} + w_same*same + w_dt*log1p(dt/60) + w_price*|logp-q|。无训练。"""
    if not hist:
        return np.zeros(0)
    s = []
    for h in hist:
        sc = 0.0
        sc += w_clk * float(h["act"] == CLK)
        sc += w_same * h["same"]
        sc += w_dt * math.log1p(h["dt"] / 60.0)
        sc += w_price * abs(h["log_price"] - q_price)
        s.append(sc)
    z = np.array(s, dtype=float)
    z = z - z.max()
    e = np.exp(z)
    return e / e.sum()


def credit(alpha: np.ndarray, gmv: float, hist: List[Dict]) -> Dict[str, float]:
    gmv = float(gmv or 0.0)
    c = alpha * gmv if len(alpha) else np.zeros(0)
    out = {
        "gmv": gmv,
        "mass_same": float((alpha * np.array([h["same"] for h in hist])).sum()) if hist else 0.0,
        "mass_clk": float((alpha * np.array([h["act"] == CLK for h in hist], dtype=float)).sum()) if hist else 0.0,
        "mass_exp": float((alpha * np.array([h["act"] == EXP for h in hist], dtype=float)).sum()) if hist else 0.0,
        "credit_same": float((c * np.array([h["same"] for h in hist])).sum()) if hist else 0.0,
        "empty_hist": float(len(hist) == 0),
    }
    return out


def attribute_cnv(evs: List[Event], cnv_item: int, cnv_ts: int, price: Optional[float]) -> Dict[str, Dict]:
    hist = history_before(evs, cnv_ts, cnv_item)
    n = len(hist)
    gmv = float(price or 0.0)
    q = _logp(price)
    plans = {
        "last_clk_item": one_hot(n, rule_last(hist, same=True, act=CLK)),
        "last_clk_any": one_hot(n, rule_last(hist, same=False, act=CLK)),
        "last_exp_item": one_hot(n, rule_last(hist, same=True, act=EXP)),
        "last_any": one_hot(n, rule_last(hist, same=False, act=None)),
        "decay_3d": time_decay(hist),
        "soft": soft_attr(hist, q_price=q),
    }
    return {name: credit(a, gmv, hist) for name, a in plans.items()} | {"n_hist": n}


def _demo() -> None:
    t0 = 1_700_000_000
    evs: List[Event] = [
        (99, EXP, t0, 10.0),
        (99, CLK, t0 + 30, 10.0),
        (11, EXP, t0 + 100, 99.0),
        (11, CLK, t0 + 120, 99.0),
        (11, CNV, t0 + 3600, 99.0),
        (22, EXP, t0 + 86400, 50.0),
        (22, CLK, t0 + 86400 + 30, 50.0),
        (22, CNV, t0 + 86400 + 40, 50.0),
    ]
    # 第一单：同品有点，last_clk_item 能分到 99
    a = attribute_cnv(evs, 11, t0 + 3600, 99.0)
    print("cnv item=11 gmv=99  n_hist", a["n_hist"])
    for k, v in a.items():
        if k == "n_hist":
            continue
        print(f"  {k:16s} mass_same={v['mass_same']:.2f} mass_clk={v['mass_clk']:.2f} credit_same={v['credit_same']:.1f}")
    # 第二单：同品有点（22），但若把 item 改成「序列里没出现的 SKU」= 空路径
    b = attribute_cnv(evs, 77, t0 + 86400 + 40, 50.0)
    print("cnv item=77（空同品） n_hist", b["n_hist"], "empty", b["last_clk_item"]["empty_hist"] or b["last_clk_item"]["mass_same"] == 0)
    print("  last_clk_item mass_same", b["last_clk_item"]["mass_same"], "last_clk_any mass_clk", b["last_clk_any"]["mass_clk"])
    assert a["last_clk_item"]["mass_same"] == 1.0
    assert b["last_clk_item"]["mass_same"] == 0.0
    assert b["last_clk_any"]["mass_clk"] == 1.0


if __name__ == "__main__":
    _demo()
