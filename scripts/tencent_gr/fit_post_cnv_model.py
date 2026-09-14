#!/usr/bin/env python3
"""转化粒预报：满窗 y_post_clk_1d。X 停在 cnv_ts。按用户切。

SKU 漏斗每天看构成率：占比 <1% 则 SKU 支路关掉（这批 ~0.11%）。
不是归因 GT。AP/AUC 只评「跟满窗后还会不会点」。

  python3 scripts/tencent_gr/fit_post_cnv_model.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, log_loss, roc_auc_score
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cross_feats import add_cross  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
POST = ROOT / "results/tencent_gr_fs150/tables/post.parquet"
OUT = ROOT / "results/tencent_gr_fs150/MODEL.md"
SKU_GATE_MIN = 0.01  # 同品 last-clk 占比低于这个，SKU 支路不进主模型

# 支路。全部 cnv_ts 已知。不准进：n_after / lift / dt_next / y_post_*
BRANCHES = {
    "heat": ["n_clk_before_7d", "n_clk_before_1d", "n_clk_before_1h"],
    "empty": ["empty_any", "lag_empty_any"],
    "any_path": ["dt_any_min"],
    "sess": ["sess_clk_before", "sess_pos"],
    "lag": ["lag_post_clk_1d_rate", "n_prior_cnv", "log1p_price"],
    "sku": ["wo_prior_clk", "dt_item_min", "item_within_5m", "n_clk_same_before"],
}
ABLATION = [
    ("heat", ["heat"]),
    ("+empty", ["heat", "empty"]),
    ("+any_path", ["heat", "empty", "any_path"]),
    ("+sess", ["heat", "empty", "any_path", "sess"]),
    ("+lag", ["heat", "empty", "any_path", "sess", "lag"]),
    ("+cross", ["heat", "empty", "any_path", "sess", "lag", "cross"]),
    ("+sku", ["heat", "empty", "any_path", "sess", "lag", "cross", "sku"]),
]


def daily_mix(post: pd.DataFrame) -> dict:
    """每天迭代一版：先打构成，再决定 SKU 闸。"""
    sku_clk = float((~post["wo_prior_clk"].fillna(True).astype(bool)).mean())
    empty = float(post["empty_any"].mean()) if "empty_any" in post.columns else float("nan")
    sess0 = float((post["sess_clk_before"].fillna(0) == 0).mean())
    return {
        "n_cnv": int(len(post)),
        "sku_clk_share": sku_clk,
        "empty_any": empty,
        "sess_clk0": sess0,
        "sku_gate": int(sku_clk >= SKU_GATE_MIN),
    }


def prep(post: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    x = pd.DataFrame(index=post.index)
    for c in cols:
        s = post[c] if c in post.columns else pd.Series(np.nan, index=post.index)
        if c in ("wo_prior_clk", "item_within_5m", "item_within_1h", "empty_any"):
            x[c] = s.fillna(False).astype(float)
        elif c in ("dt_any_min", "dt_item_min"):
            x[c + "_miss"] = s.isna().astype(float)
            x[c] = s.fillna(-1.0).astype(float)
        else:
            x[c] = pd.to_numeric(s, errors="coerce")
            # HGB 能吃 NaN；线性模型后面 fill
            x[c] = x[c]
    return x


def user_split(post: pd.DataFrame, rng: np.random.Generator, frac: float = 0.7):
    u = post["user_id"].drop_duplicates().to_numpy().copy()
    rng.shuffle(u)
    n = max(1, int(len(u) * frac))
    tr_u, te_u = set(u[:n]), set(u[n:])
    return post["user_id"].isin(tr_u), post["user_id"].isin(te_u)


def metrics(y, p) -> dict:
    y = np.asarray(y).astype(int)
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    out = {"pos": float(y.mean()), "n": int(len(y))}
    if y.min() == y.max():
        out.update(auc=float("nan"), ap=float("nan"), logloss=float("nan"))
        return out
    out["auc"] = float(roc_auc_score(y, p))
    out["ap"] = float(average_precision_score(y, p))
    out["logloss"] = float(log_loss(y, p))
    return out


def fit_pack(xtr, ytr, xte, yte) -> dict:
    xtr_f = xtr.fillna(0.0)
    xte_f = xte.fillna(0.0)
    pack = {}
    log = Pipeline(
        [("sc", StandardScaler()), ("clf", LogisticRegression(max_iter=400, C=0.5, solver="lbfgs"))]
    )
    log.fit(xtr_f, ytr)
    pack["logreg"] = metrics(yte, log.predict_proba(xte_f)[:, 1])
    clf = log.named_steps["clf"]
    pack["log_coef"] = dict(zip(xtr_f.columns.astype(str), [float(z) for z in clf.coef_.ravel()]))
    hgb = HistGradientBoostingClassifier(
        max_depth=3, max_iter=80, learning_rate=0.08, random_state=0
    )
    hgb.fit(xtr, ytr)  # 留 NaN
    pack["hgb"] = metrics(yte, hgb.predict_proba(xte)[:, 1])
    mlp = Pipeline(
        [
            ("sc", StandardScaler()),
            (
                "clf",
                MLPClassifier(
                    hidden_layer_sizes=(16, 8),
                    activation="relu",
                    max_iter=200,
                    random_state=0,
                    early_stopping=True,
                    validation_fraction=0.15,
                ),
            ),
        ]
    )
    mlp.fit(xtr_f, ytr)
    pack["mlp"] = metrics(yte, mlp.predict_proba(xte_f)[:, 1])
    return pack


def architecture_text(mix: dict, sku_on: bool) -> str:
    gate = "ON" if sku_on else "OFF (构成率<1%，只监控)"
    return f"""
```
post 转化粒                    daily mix: sku_clk={mix['sku_clk_share']:.4f}  empty={mix['empty_any']:.3f}  sess0={mix['sess_clk0']:.3f}
        |                      SKU_GATE={gate}
   满窗滤 y_post_clk_1d
   按 user_id 70/30
        |
  ┌─────────┬──────────┬───────────┬──────────┬────────────┬─────────────────┐
  │ heat    │ empty    │ any_path  │ sess     │ lag        │ sku × GATE      │
  │ 7d/1d/1h│ empty_any│ dt_any    │ pos      │ post_clk_1d│ wo_prior_clk    │
  │ n_clk   │ lag_empty│ (+miss)   │ clk_before│ n_prior    │ dt_item/within  │
  └─────────┴──────────┴───────────┴──────────┴────────────┴─────────────────┘
        | concat + User×Ctx 交叉（SKU 交叉 × GATE）
        +-- LogReg  (scale → 线性)           系数可读
        +-- HGB     (depth=3, 80 iter)       吃 NaN
        +-- MLP     (16 → 8 → 1, ReLU)       浅层非线性
        |
  Ŷ = P(满窗后 1d 内任意点击 | cnv_ts 已知)
  不是 P(哪次点击导致购买)
```
每天：dump 中间表 → 打构成 → SKU 占比跨过 {SKU_GATE_MIN:.0%} 才打开 sku 支路 → 按同一协议重训。
"""


def write_md(path: Path, mix: dict, sku_on: bool, rows: list, coef: dict) -> None:
    lines = [
        "# 转化粒模型 v1",
        "",
        "Y = 满窗 `y_post_clk_1d`。X 停在 `cnv_ts`。按用户切。SKU 支路看每天构成率。",
        architecture_text(mix, sku_on),
        "## daily mix（这批 prefix）",
        "",
        f"- n_cnv={mix['n_cnv']}  sku_clk_share=**{mix['sku_clk_share']:.4f}**  empty_any={mix['empty_any']:.3f}  sess0={mix['sess_clk0']:.3f}",
        f"- SKU_GATE={'开' if sku_on else '关'}（阈值 {SKU_GATE_MIN:.0%}）",
        "",
        "## 满窗 1d 预报（user 70/30）",
        "",
        "| 消融 | LogReg AUC | HGB AUC | MLP AUC | HGB AP |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, pack in rows:
        lg, hg, mp = pack["logreg"], pack["hgb"], pack["mlp"]
        lines.append(
            f"| {name} | {lg['auc']:.3f} | {hg['auc']:.3f} | {mp['auc']:.3f} | {hg['ap']:.3f} |"
        )
    if coef:
        lines += ["", "LogReg 系数（+cross，标准化后）：", "", "| feat | coef |", "|---|---:|"]
        for k, v in sorted(coef.items(), key=lambda kv: -abs(kv[1]))[:16]:
            lines.append(f"| `{k}` | {v:+.3f} |")
    lines += [
        "",
        "AP/AUC 只说明满窗后还会不会点预报得怎样。SKU 漏斗闭没闭合看 daily mix，不看这张表。",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    print("wrote", path)


def main() -> None:
    post = pd.read_parquet(POST)
    mix = daily_mix(post)
    sku_on = bool(mix["sku_gate"])
    print("daily mix", json.dumps(mix, indent=2))
    print(architecture_text(mix, sku_on))

    ycol = "y_post_clk_1d"
    sub = post.loc[post[ycol].notna()].copy()
    y = sub[ycol].astype(int)
    tr, te = user_split(sub, np.random.default_rng(0))
    print(f"满窗 n={len(sub)} pos={float(y.mean()):.3f} tr={int(tr.sum())} te={int(te.sum())}")

    rows = []
    last_coef = {}
    for name, brs in ABLATION:
        if "sku" in brs and not sku_on:
            if name == "+sku":
                print("skip +sku (gate off)")
            continue
        cols = [c for b in brs if b != "cross" for c in BRANCHES[b]]
        x = prep(sub, cols)
        if "cross" in brs:
            x = add_cross(x, sku_on=sku_on and "sku" in brs)
        pack = fit_pack(x.loc[tr], y.loc[tr], x.loc[te], y.loc[te])
        rows.append((name, pack))
        last_coef = pack.get("log_coef") or last_coef
        print(name, {k: pack[k] for k in ("logreg", "hgb", "mlp")})

    write_md(OUT, mix, sku_on, rows, last_coef)
    (OUT.with_suffix(".json")).write_text(
        json.dumps({"mix": mix, "sku_on": sku_on, "ablation": {n: {k: p[k] for k in ("logreg", "hgb", "mlp")} for n, p in rows}}, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
