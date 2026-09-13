#!/usr/bin/env python3
"""Automated 100+ sequence features for TAAC TencentGR-10M.

Funnel labels:
  action_type 0 = expose, 1 = click, 2 = conversion

Price / ARPU proxy: item encrypted field ``115`` (sparse numeric ~0–1000).
Coupon flags are not public → deal-sensitivity is proxied.

Example:
  python3 scripts/tencent_gr/auto_seq_features.py \\
    --root data/tencent_subset --max-users 20000 --out results/tencent_gr
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

A_EXP, A_CLK, A_CNV = 0, 1, 2
A_NAME = {0: "exp", 1: "clk", 2: "cnv"}
WINS = {"1d": 86400, "3d": 3 * 86400, "7d": 7 * 86400, "30d": 30 * 86400}
SESS_GAP = 30 * 60
EPS = 1e-6
Event = Tuple[int, int, int, Optional[float]]  # item, action, ts, price


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
        }
    a = np.asarray(xs, dtype=float)
    return {
        f"{prefix}_n": float(len(xs)),
        f"{prefix}_min": float(a.min()),
        f"{prefix}_p50": float(np.median(a)),
        f"{prefix}_mean": float(a.mean()),
        f"{prefix}_p90": float(np.percentile(a, 90)),
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
    return {
        "n_tot": float(len(sub)),
        "n_exp": float(n_exp),
        "n_clk": float(n_clk),
        "n_cnv": float(n_cnv),
        "n_uniq_item": float(len({e[0] for e in sub})),
        "ctr": rate(n_clk, n_exp),
        "cvr": rate(n_cnv, n_clk),
        "ctcvr": rate(n_cnv, n_exp),
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
                clk2cnv.append((ts - last_clk[iid]) / 60.0)
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
    return {
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
    }


def session_feats(evs: Sequence[Event]) -> Dict[str, float]:
    ss = sessions(evs)
    empty = {
        "sess_n": 0.0,
        "sess_len_mean": 0.0,
        "sess_len_p50": 0.0,
        "sess_len_max": 0.0,
        "sess_dur_min_mean": 0.0,
        "sess_depth_clk_mean": 0.0,
        "sess_bounce_rate": 0.0,
        "sess_last_len": 0.0,
        "sess_last_dur_min": 0.0,
        "sess_last_has_cnv": 0.0,
    }
    if not ss:
        return empty
    lens = [len(s) for s in ss]
    durs = [(s[-1][2] - s[0][2]) / 60.0 for s in ss]
    clk_d = [sum(1 for e in s if e[1] == A_CLK) for s in ss]
    bounce = sum(1 for s in ss if len(s) <= 1)
    last = ss[-1]
    return {
        "sess_n": float(len(ss)),
        "sess_len_mean": float(np.mean(lens)),
        "sess_len_p50": float(np.median(lens)),
        "sess_len_max": float(np.max(lens)),
        "sess_dur_min_mean": float(np.mean(durs)),
        "sess_depth_clk_mean": float(np.mean(clk_d)),
        "sess_bounce_rate": rate(bounce, len(ss)),
        "sess_last_len": float(len(last)),
        "sess_last_dur_min": float((last[-1][2] - last[0][2]) / 60.0),
        "sess_last_has_cnv": 1.0 if any(e[1] == A_CNV for e in last) else 0.0,
    }


def diversity(evs: Sequence[Event], t_end: int) -> Dict[str, float]:
    if not evs:
        return {"hist_len": 0.0, "active_days": 0.0}
    days = {e[2] // 86400 for e in evs}
    clk_items: Dict[int, int] = defaultdict(int)
    exp_cnt: Dict[int, int] = defaultdict(int)
    clk_cnt: Dict[int, int] = defaultdict(int)
    for iid, act, _, _ in evs:
        if act == A_CLK:
            clk_items[iid] += 1
            clk_cnt[iid] += 1
        elif act == A_EXP:
            exp_cnt[iid] += 1
    tot = sum(clk_items.values())
    top1 = max(clk_items.values()) / tot if tot else 0.0

    def since(act: int) -> float:
        ts = [e[2] for e in evs if e[1] == act]
        return float((t_end - max(ts)) / 3600.0) if ts else -1.0

    gaps = [evs[i][2] - evs[i - 1][2] for i in range(1, len(evs))]
    return {
        "hist_len": float(len(evs)),
        "span_day": float((evs[-1][2] - evs[0][2]) / 86400.0),
        "active_days": float(len(days)),
        "item_entropy_clk": entropy(clk_items),
        "top1_item_clk_share": float(top1),
        "hours_since_last_exp": since(A_EXP),
        "hours_since_last_clk": since(A_CLK),
        "hours_since_last_cnv": since(A_CNV),
        "reexpose_rate": rate(sum(1 for c in exp_cnt.values() if c >= 2), len(exp_cnt)),
        "repeat_clk_rate": rate(sum(1 for c in clk_cnt.values() if c >= 2), len(clk_cnt)),
        "gap_sec_p50": float(np.median(gaps)) if gaps else -1.0,
        "gap_sec_mean": float(np.mean(gaps)) if gaps else -1.0,
    }


def ui_summary(evs: Sequence[Event]) -> Dict[str, float]:
    by: Dict[int, List[Event]] = defaultdict(list)
    for e in evs:
        by[e[0]].append(e)
    pair_ctr, pair_cvr, pair_exp, clk2cnv = [], [], [], []
    multitouch = 0
    for lst in by.values():
        n_exp = sum(1 for e in lst if e[1] == A_EXP)
        n_clk = sum(1 for e in lst if e[1] == A_CLK)
        n_cnv = sum(1 for e in lst if e[1] == A_CNV)
        pair_exp.append(n_exp)
        pair_ctr.append(rate(n_clk, n_exp))
        pair_cvr.append(rate(n_cnv, n_clk))
        if n_exp >= 2 and n_cnv >= 1:
            multitouch += 1
        last_c = None
        for e in lst:
            if e[1] == A_CLK:
                last_c = e[2]
            elif e[1] == A_CNV and last_c is not None:
                clk2cnv.append((e[2] - last_c) / 60.0)
    return {
        "ui_n_items_touched": float(len(by)),
        "ui_exp_per_item_mean": float(np.mean(pair_exp)) if pair_exp else 0.0,
        "ui_pair_ctr_mean": float(np.mean(pair_ctr)) if pair_ctr else 0.0,
        "ui_pair_cvr_mean": float(np.mean(pair_cvr)) if pair_cvr else 0.0,
        "ui_multitouch_cnv_items": float(multitouch),
        "ui_pair_clk2cnv_min_p50": float(np.median(clk2cnv)) if clk2cnv else -1.0,
        "ui_pair_clk2cnv_min_mean": float(np.mean(clk2cnv)) if clk2cnv else -1.0,
    }


def build_user_features(evs: Sequence[Event], price_p50: float) -> Dict[str, float]:
    if not evs:
        return {"hist_len": 0.0, "active_days": 0.0}
    t_end = evs[-1][2]
    feats: Dict[str, float] = {}

    for k, v in funnel(evs, t_end, 10**12).items():
        feats[f"life_{k}"] = v
    for wname, wsec in WINS.items():
        for k, v in funnel(evs, t_end, wsec).items():
            feats[f"{wname}_{k}"] = v

    feats["trend_ctr_1d_minus_7d"] = feats["1d_ctr"] - feats["7d_ctr"]
    feats["trend_cvr_1d_minus_7d"] = feats["1d_cvr"] - feats["7d_cvr"]
    feats["trend_cnv_7d_vs_30d_daily"] = feats["7d_n_cnv"] / 7.0 - feats["30d_n_cnv"] / 30.0

    feats.update(session_feats(evs))
    feats.update(diversity(evs, t_end))
    feats.update(attribution(evs))

    mon = monetization(evs, price_p50)
    ad = feats.get("active_days", 0.0)
    mon["arpu_per_active_day_proxy"] = mon["arpu_sum_proxy"] / ad if ad > 0 else 0.0
    mon["deal_hesitation_exp_before_cnv"] = feats.get("attr_exp_before_cnv_mean", 0.0)
    feats.update(mon)
    feats.update(ui_summary(evs))

    trans: Dict[Tuple[int, int], int] = defaultdict(int)
    for a, b in zip(evs, evs[1:]):
        trans[(a[1], b[1])] += 1
    tot = sum(trans.values()) + EPS
    for a in range(3):
        for b in range(3):
            feats[f"trans_{A_NAME[a]}_to_{A_NAME[b]}"] = trans[(a, b)] / tot
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
        if v is None or (isinstance(v, float) and isinstance(v, float) and math.isnan(v)):
            feats[f"u_side_{col}_len"] = 0.0
            feats[f"u_side_{col}_ent"] = 0.0
            continue
        try:
            lst = list(v) if v is not None else []
        except Exception:
            lst = []
        feats[f"u_side_{col}_len"] = float(len(lst))
        cnt: Dict[Any, int] = defaultdict(int)
        for x in lst:
            cnt[x] += 1
        feats[f"u_side_{col}_ent"] = entropy(cnt)
    return feats


JUSTIFY = {
    "pay_": "付费次数/是否付费用户：conversion(action=2) 强度。",
    "arpu_": "ARPU proxy：转化上 item.115 的 sum/mean/p50；日均 / active_days。",
    "deal_": "优惠敏感：低价成交占比、相对自身更便宜、成交前曝光犹豫、浏览-成交价差。",
    "attr_clk2cnv": "同 item 点击→转化分钟数（归因时长）。",
    "attr_exp2clk": "曝光→点击分钟数。",
    "attr_": "漏斗归因时延与无点击转化率。",
    "sess_": "30min session：当下注意力。",
    "1d_": "近1天漏斗=即时意图。",
    "3d_": "近3天漏斗。",
    "7d_": "近7天漏斗=短期偏好。",
    "30d_": "近30天漏斗=稳定画像。",
    "life_": "全序列漏斗。",
    "u_side_": "官方加密用户侧：桶+list长度/熵；missingness 入模。",
    "ui_": "U×I 汇总与 pair 归因。",
    "trans_": "曝/点/转 bigram 路径偏好。",
    "trend_": "短窗相对长窗趋势。",
}


def catalog(cols: Sequence[str]) -> List[Dict[str, str]]:
    rows = []
    for k in cols:
        just = "seq-derived automated funnel/session stat"
        for pref, text in JUSTIFY.items():
            if k.startswith(pref):
                just = text
                break
        rows.append({"name": k, "family": k.split("_")[0], "justification": just})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/tencent_subset"))
    ap.add_argument("--max-users", type=int, default=20000)
    ap.add_argument("--out", type=Path, default=Path("results/tencent_gr"))
    args = ap.parse_args()

    seq_paths = parquets(args.root / "seq")
    if not seq_paths:
        raise SystemExit(f"no seq parquet under {args.root}/seq")

    print("loading item price proxy (115)...")
    price = load_price(args.root / "item_feat")
    pvals = np.asarray(list(price.values()), dtype=float) if price else np.asarray([])
    g_p50 = float(np.median(pvals)) if pvals.size else 0.0
    print(f"  priced_items={len(price)} global_p50={g_p50:.2f}")

    rows: List[Dict[str, float]] = []
    uids: List[int] = []
    print(f"featurizing <= {args.max_users} users ...")
    for uid, seq in iter_users(seq_paths, args.max_users):
        feats = build_user_features(parse_seq(seq, price), g_p50)
        feats["user_id"] = float(uid)
        rows.append(feats)
        uids.append(uid)
        if len(rows) % 5000 == 0:
            print(f"  {len(rows)} ...")

    feat_df = pd.DataFrame(rows)
    print("joining user_feat ...")
    udf = load_users(args.root / "user_feat", set(uids)).set_index("user_id")
    enriched = []
    for _, r in feat_df.iterrows():
        uid = int(r["user_id"])
        urec = udf.loc[uid] if uid in udf.index else None
        d = {k: float(v) for k, v in r.items() if k != "user_id"}
        d = enrich_user_table(d, urec)
        d["user_id"] = uid
        enriched.append(d)
    out_df = pd.DataFrame(enriched)
    feat_cols = [c for c in out_df.columns if c != "user_id"]
    print(f"feature_dim={len(feat_cols)}")

    args.out.mkdir(parents=True, exist_ok=True)
    out_path = args.out / "user_seq_features.parquet"
    out_df.to_parquet(out_path, index=False)
    (args.out / "feature_catalog.json").write_text(
        json.dumps(catalog(feat_cols), indent=2, ensure_ascii=False), encoding="utf-8"
    )

    key = [
        "pay_cnt", "arpu_sum_proxy", "arpu_mean_proxy", "arpu_per_active_day_proxy",
        "deal_cheap_cnv_share", "deal_hesitation_exp_before_cnv", "deal_vs_self_price_share",
        "attr_clk2cnv_min_p50", "attr_clk2cnv_min_mean", "attr_exp2clk_min_p50",
        "1d_ctr", "7d_cvr", "life_ctcvr", "sess_n", "sess_bounce_rate",
    ]
    summary = {
        "n_users": int(len(out_df)),
        "n_features": len(feat_cols),
        "global_price_p50_proxy": g_p50,
        "pay_user_rate": float((out_df["pay_cnt"] > 0).mean()) if "pay_cnt" in out_df else None,
        "key_stats": {},
    }
    for k in key:
        if k not in out_df.columns:
            continue
        s = out_df[k].replace(-1, np.nan)
        summary["key_stats"][k] = {
            "mean": float(np.nanmean(s)),
            "p50": float(np.nanmedian(s)),
            "p90": float(np.nanpercentile(s.dropna(), 90)) if s.notna().any() else None,
            "nonzero_rate": float((out_df[k] > 0).mean()),
        }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    md = f"""# TencentGR-10M automated sequence features

- users: **{len(out_df)}**
- feature dim: **{len(feat_cols)}**
- funnel: expose=0 / click=1 / conversion=2
- price proxy: item `115` (global p50={g_p50:.1f})

## ARPU / 付费次数 / 优惠敏感度

| 业务概念 | 特征 | 公式 | 说明 |
|---|---|---|---|
| 付费次数 | `pay_cnt` | `#conversion` | action=2 作 pay intensity |
| ARPU | `arpu_sum/mean/p50_proxy` | 转化上 `item.115` | 客单排序强度，非字面人民币 |
| 日均 ARPU | `arpu_per_active_day_proxy` | `arpu_sum / active_days` | 去掉「活得久」混杂 |
| 优惠敏感 | `deal_cheap_cnv_share` | 转化价≤全局p50 占比 | 爱买相对便宜货 |
| 优惠敏感 | `deal_vs_self_price_share` | 低于自身成交中位价占比 | 等更便宜再买 |
| 优惠敏感 | `deal_hesitation_exp_before_cnv` | 转化前同item曝光次数 | 比价/等券 |
| 点击归因 | `attr_clk2cnv_min_*` | 同item上次点击→转化分钟 | 你要的归因时长 |

## 为什么 user 侧要进模

1. `103–110` 是官方用户先验，丢掉白扔信号。
2. 标量桶做人口/等级分桶，利于 CTR/CVR 校准。
3. list → 长度+熵 = 兴趣广度，无需明文。
4. `*_known` 缺失本身可预测（新/低活）。
5. 侧信息=静态长期，seq=动态行为，正交互补。

## Key stats

```json
{json.dumps(summary["key_stats"], indent=2)}
```

Catalog: `{args.out / "feature_catalog.json"}` ({len(feat_cols)} dims).
"""
    (args.out / "FEATURE_REPORT.md").write_text(md, encoding="utf-8")
    docs = Path("docs/tencent_gr")
    docs.mkdir(parents=True, exist_ok=True)
    (docs / "TencentGR_auto_seq_features.md").write_text(md, encoding="utf-8")
    print(json.dumps(summary, indent=2)[:2500])
    print("wrote", out_path)


if __name__ == "__main__":
    main()
