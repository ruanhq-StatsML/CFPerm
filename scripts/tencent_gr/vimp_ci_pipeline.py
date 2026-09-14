#!/usr/bin/env python3
"""Y 可换的 VIMP + bootstrap CI。解耦 user-mixture 和 CTR。

R 不当检验。Insight = 各 Y 的 VIMP 序 + 用户簇 bootstrap 分位区间。

口径
----
事件先炸成 (user_id, item_id, act, ts)。每人时间中位切开：
  左窗 → X（量 / share / 衰减 / 场）
  右窗 → Y（可换）
时钟 W = 1{用户 t_end > 全局中位}。RF-domain 不看 Y，是 P(W|X)=谁来了。
PO-VIMP 看 φ=(Y-μ)(W-e)，是「早/晚对这个 Y 的图」。

Y 注册（右窗；n_exp=0 的 CTR 丢掉，不当 0）：
  n_clk   点击次数（量，易和曝光次数缠）
  ctr     n_clk/n_exp，至少 MIN_EXP 次曝光（每次曝光后点的比例）
  n_cnv   成交次数
  cvr     n_cnv/n_clk，至少 1 次点击
  any_cnv 1{右窗有成交}

  python3 scripts/tencent_gr/vimp_ci_pipeline.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from block_tables import tab_decay, tab_funnel, tab_session  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
EV = ROOT / "results/tencent_gr_fs150/tables/ev.parquet"
OUT = ROOT / "results/tencent_gr_fs150"
SEED = 0
B = 25
TREES = 25
DEPTH = 5
FOLDS = 3
MIN_EXP = 3
SPLIT = 0.5
LO = 0.025
HI = 0.975
# impurity 总 ≥0；CI 下界过这条才算「这把尺子选中」
VIMP_FLOOR = 0.03

X_KEEP = [
    "hist_len",
    "active_days",
    "life_n_exp",
    "life_n_clk",
    "pay_cnt",
    "life_cnv_share",
    "life_clk_share",
    "life_ctr",
    "dec_hl7d_dec_clk",
    "dec_hl7d_dec_cnv",
    "sess_bounce_rate",
    "sess_depth_cnv_mean",
    "sess_cnv_sess_rate",
    "ui_only_exp_share",
]

LOGIC = {
    "hist_len": "tenure/量：看了多久（和删失缠）",
    "active_days": "tenure：来了几天",
    "life_n_exp": "量：左窗曝光次数",
    "life_n_clk": "量：左窗点击次数",
    "pay_cnt": "量：左窗成交次数（不是人）",
    "life_cnv_share": "结构：轨迹里买占多大",
    "life_clk_share": "结构：轨迹里点占多大",
    "life_ctr": "左窗 CTR（过去点率 → 右窗，不是泄漏）",
    "dec_hl7d_dec_clk": "近期点击热度",
    "dec_hl7d_dec_cnv": "近期成交热度",
    "sess_bounce_rate": "场碎",
    "sess_depth_cnv_mean": "场内成交深度",
    "sess_cnv_sess_rate": "多少场真正成交",
    "ui_only_exp_share": "只曝不点的货占比（队列）",
}


def tab_ui_only(ev: pd.DataFrame) -> pd.DataFrame:
    mx = ev.groupby(["user_id", "item_id"], sort=False)["act"].max()
    share = mx.eq(0).groupby(level=0).mean().rename("ui_only_exp_share")
    return share.reset_index()


def assemble(ev: pd.DataFrame) -> pd.DataFrame:
    """左 X、右 Y、时钟 t_end。"""
    mm = ev.groupby("user_id")["ts"].agg(t0="min", t_end="max").reset_index()
    mm["t_cut"] = mm["t0"] + (mm["t_end"] - mm["t0"]) * SPLIT
    ev = ev.merge(mm, on="user_id", how="left")
    left = ev.loc[ev["ts"] <= ev["t_cut"]].drop(columns=["t0", "t_end", "t_cut"])
    right = ev.loc[ev["ts"] > ev["t_cut"]]
    x = tab_funnel(left).merge(tab_decay(left), on="user_id")
    x = x.merge(tab_session(left), on="user_id", how="left")
    x = x.merge(tab_ui_only(left), on="user_id", how="left")
    x["ui_only_exp_share"] = x["ui_only_exp_share"].fillna(1.0)
    for c in ("sess_bounce_rate", "sess_depth_cnv_mean", "sess_cnv_sess_rate"):
        if c in x.columns:
            x[c] = x[c].fillna(0.0)
    y = right.groupby("user_id").agg(
        y_n_exp=("act", lambda a: int((a == 0).sum())),
        y_n_clk=("act", lambda a: int((a == 1).sum())),
        y_n_cnv=("act", lambda a: int((a == 2).sum())),
    ).reset_index()
    y["y_ctr"] = np.where(y["y_n_exp"] >= MIN_EXP, y["y_n_clk"] / y["y_n_exp"], np.nan)
    y["y_cvr"] = np.where(y["y_n_clk"] >= 1, y["y_n_cnv"] / y["y_n_clk"], np.nan)
    y["y_any_cnv"] = (y["y_n_cnv"] > 0).astype(float)
    out = mm.merge(x, on="user_id", how="inner").merge(y, on="user_id", how="left")
    for c in ["y_n_exp", "y_n_clk", "y_n_cnv", "y_any_cnv"]:
        out[c] = out[c].fillna(0.0)
    return out


TARGETS = [
    {"name": "n_clk", "col": "y_n_clk", "kind": "count", "need": None, "ask": "右窗点了几次（量）"},
    {"name": "ctr", "col": "y_ctr", "kind": "rate", "need": "y_ctr", "ask": "右窗每次曝光后点的比例"},
    {"name": "n_cnv", "col": "y_n_cnv", "kind": "count", "need": None, "ask": "右窗买了几次（量）"},
    {"name": "cvr", "col": "y_cvr", "kind": "rate", "need": "y_cvr", "ask": "右窗点了之后买的比例"},
    {"name": "any_cnv", "col": "y_any_cnv", "kind": "binary", "need": None, "ask": "右窗有没有成交"},
]


def _xy(df: pd.DataFrame, ycol: str, need: str | None):
    sub = df if need is None else df.loc[df[need].notna()]
    names = [c for c in X_KEEP if c in sub.columns]
    X = sub[names].fillna(0.0).to_numpy(np.float64)
    y = pd.to_numeric(sub[ycol], errors="coerce").fillna(0.0).to_numpy(np.float64)
    med = float(np.median(df["t_end"].to_numpy(np.float64)))
    w = (sub["t_end"].to_numpy(np.float64) > med).astype(int)
    return sub, names, X, y, w


def rf_domain_vimp(X, W, *, seed: int):
    clf = RandomForestClassifier(
        n_estimators=TREES, max_depth=DEPTH, min_samples_leaf=8, random_state=seed, n_jobs=1
    )
    n = len(W)
    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    rng.shuffle(idx)
    cut = max(1, n * 3 // 4)
    tr, te = idx[:cut], idx[cut:]
    clf.fit(X[tr], W[tr])
    auc = float("nan")
    if len(np.unique(W[te])) > 1:
        auc = float(roc_auc_score(W[te], clf.predict_proba(X[te])[:, 1]))
    # 全样本 VIMP（CI 另用 bootstrap）
    clf.fit(X, W)
    return auc, clf.feature_importances_.astype(float)


def po_vimp(X, Y, W, *, seed: int):
    if len(np.unique(W)) < 2:
        return np.zeros(X.shape[1]), 0.0
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    cv = StratifiedKFold(n_splits=min(FOLDS, int(W.sum()), int((1 - W).sum()) or 1), shuffle=True, random_state=seed)
    try:
        splits = list(cv.split(X, W))
    except ValueError:
        return np.zeros(X.shape[1]), 0.0
    for fold, (tr, te) in enumerate(splits):
        m = RandomForestRegressor(
            n_estimators=TREES, max_depth=DEPTH, min_samples_leaf=8, random_state=seed + fold, n_jobs=1
        )
        e = RandomForestClassifier(
            n_estimators=TREES, max_depth=DEPTH, min_samples_leaf=8, random_state=seed + 40 + fold, n_jobs=1
        )
        m.fit(X[tr], Y[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, 0.05, 0.95)
    po = (Y - m_hat) * (W.astype(float) - e_hat)
    tau = RandomForestRegressor(
        n_estimators=TREES, max_depth=DEPTH, min_samples_leaf=8, random_state=seed + 7, n_jobs=1
    )
    tau.fit(X, po)
    r = float(np.mean(tau.predict(X) ** 2))
    return tau.feature_importances_.astype(float), r


def boot_ci(fn, X, y, w, *, seed: int) -> np.ndarray:
    """用户簇：对行下标有放回。返回 (B, p)。"""
    rng = np.random.default_rng(seed)
    n, p = X.shape
    out = np.zeros((B, p))
    k = 0
    tries = 0
    while k < B and tries < B * 4:
        tries += 1
        ix = rng.integers(0, n, size=n)
        ww = w[ix]
        if ww.min() == ww.max():
            continue
        try:
            v = fn(X[ix], y[ix], ww, seed=seed + tries)
        except Exception:
            continue
        out[k] = v
        k += 1
    return out[:k]


def summarize(names, point, boots) -> list[dict]:
    if len(boots) == 0:
        boots = point.reshape(1, -1)
    lo = np.quantile(boots, LO, axis=0)
    hi = np.quantile(boots, HI, axis=0)
    rows = []
    for i, n in enumerate(names):
        rows.append(
            {
                "name": n,
                "logic": LOGIC.get(n, ""),
                "vimp": float(point[i]),
                "lo": float(lo[i]),
                "hi": float(hi[i]),
                "selected": int(lo[i] >= VIMP_FLOOR),
            }
        )
    rows.sort(key=lambda r: -r["vimp"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


def insights(mix_rows: list[dict], by_y: dict) -> list[str]:
    mix_sel = {r["name"] for r in mix_rows if r["selected"]}
    ctr_sel = {r["name"] for r in by_y.get("ctr", {}).get("rows", []) if r["selected"]}
    nclk_sel = {r["name"] for r in by_y.get("n_clk", {}).get("rows", []) if r["selected"]}
    cnv_sel = {r["name"] for r in by_y.get("any_cnv", {}).get("rows", []) if r["selected"]}
    cards = []
    only_mix = sorted(mix_sel - ctr_sel)
    only_ctr = sorted(ctr_sel - mix_sel)
    count_not_rate = sorted(nclk_sel - ctr_sel)
    cards.append(
        f"mixture∩CTR 选中 = {sorted(mix_sel & ctr_sel) or '∅'}；"
        f"只在 mixture = {only_mix or '∅'}；只在 CTR = {only_ctr or '∅'}。"
    )
    if only_mix:
        cards.append(
            "只被 RF-domain 选中的列在回答「后来的人是谁」（队列/tenure/量），"
            "不是「给定这些人之后怎么点」。解耦的就是这一刀。"
        )
    if "pay_cnt" in mix_sel and "pay_cnt" not in ctr_sel:
        cards.append("`pay_cnt`：mixture 选中、CTR 未选中 → 次数代表人/谁来了，不是点击图。")
    if "pay_cnt" not in mix_sel:
        cards.append("`pay_cnt` 在 mixture 下界也没过线：这根时钟上人身上的量更可能在 `hist_len`/`life_n_exp`/`ui_only_exp_share`。")
    if "life_cnv_share" in (ctr_sel | cnv_sel) and "pay_cnt" not in (ctr_sel | cnv_sel):
        cards.append("`life_cnv_share` 进了点/买的图，`pay_cnt` 没有 → 结构不是次数。")
    if "dec_hl7d_dec_clk" in ctr_sel:
        cards.append("`dec_hl7d_dec_clk` 进 CTR 图：近期点击热度，不是 lifetime 单量。")
    if count_not_rate:
        cards.append(
            f"点次数选中但 CTR 未选中 {count_not_rate}：量跟着曝光走，换成每次曝光后的比例就掉了。"
        )
    if "ui_only_exp_share" in mix_sel and "ui_only_exp_share" not in ctr_sel:
        cards.append("`ui_only_exp_share` 只在 mixture：队列（只曝不点谁来了），不是 CTR 图。")
    elif "ui_only_exp_share" in mix_sel and "ui_only_exp_share" in ctr_sel:
        cards.append("`ui_only_exp_share` 两板都过线：后来的人只曝不点，右窗点率也还咬着这列——队列和 CTR 没完全拆开。")
    if "dec_hl7d_dec_clk" in ctr_sel and "dec_hl7d_dec_cnv" in (cnv_sel | {r["name"] for r in by_y.get("n_cnv", {}).get("rows", []) if r["selected"]}):
        cards.append("CTR 选中点击衰减、成交 Y 选中成交衰减：点图和买图不是同一套热度。")
    if "life_cnv_share" in mix_sel and "life_cnv_share" not in ctr_sel:
        cards.append("`life_cnv_share` 只在 mixture：轨迹里买占多大是谁来了，不是右窗怎么点。")
    return cards


def fmt_table(rows: list[dict]) -> str:
    lines = ["| rank | feat | 逻辑 | VIMP | 95% CI | 选中 |", "|---|---|---|---:|---|:---:|"]
    for r in rows:
        mark = "yes" if r["selected"] else ""
        lines.append(
            f"| {r['rank']} | `{r['name']}` | {r['logic']} | {r['vimp']:.3f} | "
            f"[{r['lo']:.3f}, {r['hi']:.3f}] | {mark} |"
        )
    return "\n".join(lines)


def main() -> None:
    ev = pd.read_parquet(EV)
    print(f"ev {len(ev)} users {ev.user_id.nunique()}", flush=True)
    df = assemble(ev)
    print(f"assembled n={len(df)} ctr_ok={int(df.y_ctr.notna().sum())} cvr_ok={int(df.y_cvr.notna().sum())}", flush=True)

    # mixture：所有人，Y 不用
    names = [c for c in X_KEEP if c in df.columns]
    med = float(np.median(df["t_end"].to_numpy(np.float64)))
    X_all = df[names].fillna(0.0).to_numpy(np.float64)
    w_all = (df["t_end"].to_numpy(np.float64) > med).astype(int)
    auc, v_mix = rf_domain_vimp(X_all, w_all, seed=SEED)
    print(f"RF-domain AUC={auc:.3f}", flush=True)

    def mix_fn(X, y, w, seed):
        # y ignored
        _, v = rf_domain_vimp(X, w, seed=seed)
        return v

    print(f"bootstrap mixture B={B} ...", flush=True)
    mix_boot = boot_ci(mix_fn, X_all, w_all, w_all, seed=SEED + 1)
    mix_rows = summarize(names, v_mix, mix_boot)

    by_y = {}
    for t in TARGETS:
        sub, nm, X, y, w = _xy(df, t["col"], t["need"])
        print(f"Y={t['name']} n={len(sub)} pos/mean={float(y.mean()):.4f} W1={w.mean():.3f}", flush=True)
        v, r = po_vimp(X, y, w, seed=SEED + 11)
        print(f"  PO-risk(desc)={r:.6f}  (不是检验)", flush=True)

        def fn(X_, y_, w_, seed, _nm=nm):
            vv, _ = po_vimp(X_, y_, w_, seed=seed)
            return vv

        boots = boot_ci(fn, X, y, w, seed=SEED + 100 + TARGETS.index(t) * 17)
        rows = summarize(nm, v, boots)
        by_y[t["name"]] = {"ask": t["ask"], "kind": t["kind"], "n": int(len(sub)), "mean": float(y.mean()), "po_risk_desc": r, "rows": rows}

    cards = insights(mix_rows, by_y)
    payload = {
        "split": SPLIT,
        "min_exp": MIN_EXP,
        "B": B,
        "vimp_floor": VIMP_FLOOR,
        "rf_domain_auc": auc,
        "mixture": mix_rows,
        "targets": {k: {kk: vv for kk, vv in rec.items() if kk != "rows"} | {"top": rec["rows"][:8]} for k, rec in by_y.items()},
        "cards": cards,
    }
    # keep full rows in a slimmer dump
    payload["targets_full"] = {k: rec["rows"] for k, rec in by_y.items()}
    (OUT / "VIMP_CI.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md = [
        "# VIMP + bootstrap CI：换 Y，解耦 mixture / CTR",
        "",
        "R 不当检验。用户簇 bootstrap 给的是 **VIMP 的分位区间**。",
        f"左窗 X、右窗 Y（prefix 内时间 {SPLIT:.0%} 切开）。W=用户 t_end 中位。CTR 分母 `n_exp≥{MIN_EXP}`，缺测丢掉。",
        f"RF-domain AUC **{auc:.3f}**（P(X) 谁来了）。选中规则：CI 下界 ≥ {VIMP_FLOOR}。B={B}。",
        "",
        "## 自动 insights",
        "",
    ]
    for c in cards:
        md.append(f"- {c}")
    md += ["", "## mixture（RF-domain，与 Y 无关）", "", fmt_table(mix_rows), ""]
    for t in TARGETS:
        rec = by_y[t["name"]]
        md += [
            f"## Y = `{t['name']}` — {t['ask']}",
            "",
            f"n={rec['n']} mean={rec['mean']:.4f}  PO-risk(描述)={rec['po_risk_desc']:.6f}",
            "",
            fmt_table(rec["rows"]),
            "",
        ]
    md += [
        "## 标准化机制",
        "",
        "同一套 (左 X, 时钟 W, 用户 bootstrap)。只换右窗 Y，再跑 PO-VIMP+CI。",
        "和 mixture 板求交/差：只在 mixture = 人/队列；只在 CTR = 点击图；点次数有、CTR 无 = 量和曝光缠在一起。",
        "`pay_cnt` vs `life_cnv_share` vs `dec_hl7d_*` 三列对照着看：量 / 结构 / 近期热度。",
        "",
    ]
    (OUT / "VIMP_CI.md").write_text("\n".join(md), encoding="utf-8")
    print("cards:")
    for c in cards:
        print(" -", c)
    print("wrote", OUT / "VIMP_CI.md")


if __name__ == "__main__":
    main()
