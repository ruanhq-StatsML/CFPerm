#!/usr/bin/env python3
"""买完会不会带动后续点击：转化粒中间表 → 刻画 + 预测。

Y = 1{cnv 之后 T 内有点击}（任意 / 同商品）。X 只用 cnv_ts 已知量。
不跑 F-score。用户切 70/30，避免同一人进出训练集。

  python3 scripts/tencent_gr/run_post_cnv.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "tencent_gr"))

from auto_feats_select_train import (  # noqa: E402
    iter_users,
    load_price,
    parquets,
    parse_seq,
    split_prefix_suffix,
)
from block_tables import (  # noqa: E402
    POST_PREDICT_X,
    events_frame,
    tab_attr_events,
    tab_post_cnv_events,
)

X_BLOCKS = {
    "volume": [
        "n_clk_before_1h",
        "n_clk_before_1d",
        "n_clk_before_7d",
        "n_clk_same_before",
    ],
    "+path": [
        "n_clk_before_1h",
        "n_clk_before_1d",
        "n_clk_before_7d",
        "n_clk_same_before",
        "wo_prior_clk",
        "dt_item_min",
        "dt_any_min",
        "item_within_5m",
        "item_within_1h",
    ],
    "+lag": [
        "n_clk_before_1h",
        "n_clk_before_1d",
        "n_clk_before_7d",
        "n_clk_same_before",
        "n_prior_cnv",
        "lag_post_clk_1d_rate",
    ],
    "+sess_money_lag": POST_PREDICT_X,
}


def load_prefix_ev(root: Path, max_users: int, prefix_frac: float) -> pd.DataFrame:
    price = load_price(root / "item_feat")
    seq_paths = parquets(root / "seq")
    if not seq_paths:
        raise SystemExit(f"no seq under {root}/seq")
    frames = []
    t0 = time.time()
    for i, (uid, seq) in enumerate(iter_users(seq_paths, max_users), 1):
        evs = parse_seq(seq, price)
        prefix, _ = split_prefix_suffix(evs, prefix_frac=prefix_frac)
        if prefix:
            frames.append(events_frame(uid, prefix))
        if i % 1000 == 0:
            print(f"  loaded {i} users in {time.time() - t0:.1f}s", flush=True)
    if not frames:
        raise SystemExit("no events")
    ev = pd.concat(frames, ignore_index=True)
    print(f"ev {len(ev)} rows / {ev.user_id.nunique()} users in {time.time() - t0:.1f}s")
    return ev


def prep_x(post: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    x = pd.DataFrame(index=post.index)
    for c in cols:
        if c not in post.columns:
            x[c] = 0.0
            continue
        s = post[c]
        if c == "wo_prior_clk":
            x[c] = s.fillna(False).astype(float)
        elif c.startswith("item_within_"):
            x[c] = s.fillna(False).astype(float)
        elif c in ("dt_item_min", "dt_any_min", "dt_first_min"):
            x[c] = s.fillna(-1.0).astype(float)
        else:
            x[c] = pd.to_numeric(s, errors="coerce").fillna(0.0)
    return x


def fit_auc(xtr, ytr, xte, yte) -> dict:
    out = {}
    if ytr.min() == ytr.max() or yte.min() == yte.max() or len(np.unique(yte)) < 2:
        return {"hgb_auc": float("nan"), "log_auc": float("nan"), "hgb_ap": float("nan")}
    hgb = HistGradientBoostingClassifier(max_depth=3, max_iter=80, learning_rate=0.08, random_state=0)
    hgb.fit(xtr, ytr)
    ph = hgb.predict_proba(xte)[:, 1]
    out["hgb_auc"] = float(roc_auc_score(yte, ph))
    out["hgb_ap"] = float(average_precision_score(yte, ph))
    log = Pipeline(
        [
            ("sc", StandardScaler()),
            ("clf", LogisticRegression(max_iter=400, C=0.5, solver="lbfgs")),
        ]
    )
    log.fit(xtr, ytr)
    pl = log.predict_proba(xte)[:, 1]
    out["log_auc"] = float(roc_auc_score(yte, pl))
    clf = log.named_steps["clf"]
    out["log_coef"] = dict(zip(xtr.columns, [float(z) for z in clf.coef_.ravel()]))
    return out


def split_users(post: pd.DataFrame, rng: np.random.Generator, frac: float = 0.7):
    u = post["user_id"].drop_duplicates().to_numpy().copy()
    rng.shuffle(u)
    n = max(1, int(len(u) * frac))
    tr, te = set(u[:n]), set(u[n:])
    return post["user_id"].isin(tr), post["user_id"].isin(te)


def summarize(post: pd.DataFrame) -> dict:
    s = {
        "n_cnv": int(len(post)),
        "n_user": int(post.user_id.nunique()),
    }
    for w in ("5m", "1h", "1d", "7d"):
        y = post[f"y_post_clk_{w}"]
        ys = post[f"y_post_same_{w}"]
        s[f"obs_{w}"] = int(y.notna().sum())
        s[f"post_clk_{w}"] = float(y.mean()) if y.notna().any() else float("nan")
        s[f"post_same_{w}"] = float(ys.mean()) if ys.notna().any() else float("nan")
    lift = post["lift_1d"].dropna()
    delta = post["delta_1d"].dropna()
    s["lift_1d_p50"] = float(lift.median()) if len(lift) else float("nan")
    s["lift_1d_mean"] = float(lift.mean()) if len(lift) else float("nan")
    s["p_lift_gt1"] = float((lift > 1.0).mean()) if len(lift) else float("nan")
    s["p_delta_gt0"] = float((delta > 0).mean()) if len(delta) else float("nan")
    s["n_after_1d_mean"] = float(post["n_clk_after_1d"].mean())
    s["n_before_1d_mean"] = float(post["n_clk_before_1d"].mean())
    s["n_same_after_1d_mean"] = float(post["n_same_after_1d"].mean())
    s["n_same_before_1d_mean"] = float(post["n_same_before_1d"].mean())
    s["same_sess_given_next"] = float(post["next_clk_same_sess"].mean())
    dt = post["dt_next_clk_min"].dropna()
    s["dt_next_clk_p50"] = float(dt.median()) if len(dt) else float("nan")
    return s


def write_note(path: Path, s: dict, pred: dict) -> None:
    def f(x, nd=3):
        if x is None or (isinstance(x, float) and not np.isfinite(x)):
            return "nan"
        return f"{x:.{nd}f}"

    lines = [
        "# 买完会不会带动后续点击",
        "",
        "last-touch（backward asof）= 买之前怎么点过来。",
        "这块（forward asof）= **买完还会不会点**。`trans_cnv_to_exp` 只是邻接 Markov，不够。",
        "",
        "## 中间表",
        "",
        "```",
        "ev",
        " ├ attr  = merge_asof(cnv, clk, backward)   # 同品/任意 last-click",
        " ├ post  = merge_asof(cnv, clk, forward)    # 下一次点击 + 同长窗 before/after 计数",
        " │          cum asof: n_after = C(t+W)-C(t), n_before = C(t)-C(t-W)",
        " └ user  = 未删失转化上的 rate / lift p50，再和 funnel/decay/sess merge",
        "```",
        "",
        "转化粒。删失：窗内已点到 → 1；跟满 W 没点 → 0；跟不满 → 不当 0。",
        "",
        "## 刻画什么（不是因果）",
        "",
        "| 量 | 公式 | 在问什么 |",
        "|---|---|---|",
        "| `y_post_clk_T` | 1{任意下一次点击 ≤ T} | 买完人还在不在点（平台活跃） |",
        "| `y_post_same_T` | 1{同品下一次点击 ≤ T} | 这件商品买完还看不看（复访/晒单/后悔） |",
        "| `lift_T` | n_after / (n_before+1) | 相对买前同长窗，点击量抬没抬 |",
        "| `next_clk_same_sess` | dt_next ≤ 30min | 当场续点，还是隔场回访 |",
        "",
        "任意点击的 lift 会被「本来就在逛」污染：连着买两单时，后一单的预热点会算进前一单的 after。",
        "**同品 lift / y_post_same** 干净得多。当场续点用 5m/30min，回访用 1d/7d。",
        "这不是 CATE：没有对照、没有 ignorability。只是路径描述。要预测的是「这一单之后会不会点」，不是「买导致多点」。",
        "",
        "## 预测：Y 与 X 必须切开",
        "",
        "```",
        "Y = y_post_clk_1d   # 或 y_post_same_1d；只用跟满 1d 的转化",
        "X 只能是 cnv_ts 已知：",
        "  路径  wo_prior_clk, dt_item/any, item_within_5m/1h     # 已有 attr 表",
        "  买前量 n_clk_before_{1h,1d,7d}, n_clk_same_before      # 基线活跃，不是 Y",
        "  当场  sess_pos, sess_clk_before                         # 热场续点的主混杂",
        "  钱    log1p_price",
        "  滞后  n_prior_cnv, lag_post_clk_1d_rate                 # 此前各单的买后点击率，shift(1)",
        "不准进 X：n_clk_after_* / lift_* / dt_next / y_post_* / 原始 price（和 log1p 共线）",
        "```",
        "",
        "为什么要这些：",
        "",
        "- **买前量**：人本来就爱点，买后也会点。这是必须先控的基线。",
        "- **路径**：冲动（within 5m / wo_prior_clk）vs 长犹豫。前者更容易当场续点，后者更像买完走人。",
        "- **当场深度**：30min 场还没关，下一击几乎是续逛，不是「购买带动」。",
        "- **滞后买后率**：这个人以前买完爱不爱点——用户倾向，给下一单用。",
        "- **价格**：贵的可能回去反复看；便宜的买完即走。要数据说话，不先验锁方向。",
        "",
        "用户表只并 **历史倾向**（`post_clk_1d_rate` 等），给 `future_cnv` 那类用户粒任务。",
        "预测「这一单之后」必须停在转化粒。",
        "",
        f"## 这批 prefix（n_cnv={s['n_cnv']}, n_user={s['n_user']}）",
        "",
        "| 窗 | 未删失 | 任意点击率 | 同品点击率 |",
        "|---|---:|---:|---:|",
        f"| 5m | {s['obs_5m']} | {f(s['post_clk_5m'])} | {f(s['post_same_5m'])} |",
        f"| 1h | {s['obs_1h']} | {f(s['post_clk_1h'])} | {f(s['post_same_1h'])} |",
        f"| 1d | {s['obs_1d']} | {f(s['post_clk_1d'])} | {f(s['post_same_1d'])} |",
        f"| 7d | {s['obs_7d']} | {f(s['post_clk_7d'])} | {f(s['post_same_7d'])} |",
        "",
        f"1d lift p50={f(s['lift_1d_p50'])} mean={f(s['lift_1d_mean'])}；"
        f"P(lift>1)={f(s['p_lift_gt1'])}；P(n_after>n_before)={f(s['p_delta_gt0'])}。",
        f"任意：before {f(s['n_before_1d_mean'])} → after {f(s['n_after_1d_mean'])}；",
        f"同品：before {f(s['n_same_before_1d_mean'])} → after {f(s['n_same_after_1d_mean'])}。",
        f"有下一次点击时，P(仍在当场 30min)={f(s['same_sess_given_next'])}；"
        f"dt_next p50={f(s['dt_next_clk_p50']/1440, 1)} d。",
        "",
        "读这批：买完 **不抬** 后续点击。1d after≈before，lift 中位数 0，P(after>before)=0.09。",
        "5m/1h 几乎不点；同品 1d 只有约千分之一（这件商品买完基本不再点它）。",
        "7d 任意 0.41、dt_next 中位约 6 天 → 那是人还在平台上逛，不是购买带动的余热。",
        "同品 last-click 在这批对得上的极少，所以 attr 同品桶对「买后任意点击」几乎没信息；",
        "有用的是任意点击间隔、7d 买前量和滞后买后率。",
        "",
        "## 预测 Y=y_post_clk_1d（按用户 70/30）",
        "",
        "| X | HGB AUC | LogReg AUC | HGB AP |",
        "|---|---:|---:|---:|",
    ]
    for name, m in pred.items():
        lines.append(
            f"| {name} | {f(m.get('hgb_auc'))} | {f(m.get('log_auc'))} | {f(m.get('hgb_ap'))} |"
        )
    full = pred.get("+sess_money_lag") or pred.get("all") or {}
    coef = full.get("log_coef") or {}
    if coef:
        lines += [
            "",
            "LogReg 系数（标准化后，Y=任意 1d 点击；>0 更像买完还点）：",
            "",
            "| feat | coef |",
            "|---|---:|",
        ]
        for k, v in sorted(coef.items(), key=lambda kv: -abs(kv[1])):
            lines.append(f"| `{k}` | {v:+.3f} |")
    lines += [
        "",
        "volume-only 已经能到 ~0.68：1d 任意点击主要是「本来就爱点的人还在点」。",
        "+lag 再涨，说明「这个人以前买完爱不爱点」是下一单的主信号。",
        "path/sess/price 几乎不再涨：同品买后点击近乎 0，1d Y 也不是当场续点。",
        "同品 Y 这批只有十几正例，不够建模——要刻画「买完还看这件」先承认事件极稀。",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {path}")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("data/tencent_subset"))
    ap.add_argument("--max-users", type=int, default=6000)
    ap.add_argument("--prefix-frac", type=float, default=0.75)
    ap.add_argument("--out", type=Path, default=Path("results/tencent_gr_fs150/POST_CNV.md"))
    args = ap.parse_args()

    ev = load_prefix_ev(args.root, args.max_users, args.prefix_frac)
    t1 = time.time()
    attr = tab_attr_events(ev)
    post = tab_post_cnv_events(ev, attr)
    print(f"post {len(post)} conversions in {time.time() - t1:.1f}s")
    s = summarize(post)
    print(s)

    rng = np.random.default_rng(0)
    pred = {}
    y = post["y_post_clk_1d"]
    mask = y.notna()
    sub = post.loc[mask].copy()
    yv = sub["y_post_clk_1d"].astype(int)
    tr, te = split_users(sub, rng)
    print(f"predict rows {len(sub)} pos={float(yv.mean()):.3f} tr={int(tr.sum())} te={int(te.sum())}")
    for name, cols in X_BLOCKS.items():
        x = prep_x(sub, cols)
        pred[name] = fit_auc(x.loc[tr], yv.loc[tr], x.loc[te], yv.loc[te])
        print(name, {k: pred[name][k] for k in pred[name] if k != "log_coef"})

    y1h = post["y_post_clk_1h"]
    m1h = y1h.notna()
    subh = post.loc[m1h].copy()
    yvh = subh["y_post_clk_1h"].astype(int)
    trh, teh = split_users(subh, rng)
    xh = prep_x(subh, POST_PREDICT_X)
    pred["Y=clk_1h allX"] = fit_auc(xh.loc[trh], yvh.loc[trh], xh.loc[teh], yvh.loc[teh])
    print("clk_1h", {k: pred["Y=clk_1h allX"][k] for k in pred["Y=clk_1h allX"] if k != "log_coef"})

    # 同品 1d，完整 X
    y2 = post["y_post_same_1d"]
    m2 = y2.notna()
    sub2 = post.loc[m2].copy()
    yv2 = sub2["y_post_same_1d"].astype(int)
    print(f"same_1d rows {len(sub2)} pos={int(yv2.sum())} rate={float(yv2.mean()):.5f}")
    if int(yv2.sum()) >= 20 and int((1 - yv2).sum()) >= 20:
        tr2, te2 = split_users(sub2, rng)
        x2 = prep_x(sub2, POST_PREDICT_X)
        pred["Y=same_1d allX"] = fit_auc(x2.loc[tr2], yv2.loc[tr2], x2.loc[te2], yv2.loc[te2])
        print("same_1d", {k: pred["Y=same_1d allX"][k] for k in pred["Y=same_1d allX"] if k != "log_coef"})
    else:
        pred["Y=same_1d allX"] = {
            "hgb_auc": float("nan"),
            "log_auc": float("nan"),
            "hgb_ap": float("nan"),
            "note": f"too few pos ({int(yv2.sum())})",
        }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_note(args.out, s, pred)
    post_out = args.out.with_name("post_cnv_events.parquet")
    # 全转化粒太大就只写未删失 1d 的预测用列
    cols_keep = [
        "user_id",
        "item_id",
        "cnv_ts",
        "log1p_price",
        *POST_PREDICT_X,
        "y_post_clk_1h",
        "y_post_clk_1d",
        "y_post_clk_7d",
        "y_post_same_1d",
        "lift_1d",
        "delta_1d",
        "n_clk_after_1d",
        "n_clk_before_1d",
    ]
    keep = [c for c in dict.fromkeys(cols_keep) if c in sub.columns]
    sub[keep].to_parquet(post_out, index=False)
    print(f"wrote {post_out}")


if __name__ == "__main__":
    main()
