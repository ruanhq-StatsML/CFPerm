#!/usr/bin/env python3
"""W1 vs W2 compare → feature selection (落地向).

两段时间窗（≥1 个月间隔）各自 ``feature_engineer(t_start,t_end)``，
再 **直接比较 W1 与 W2**，用差异做 feature selection —— 选出驱动时段漂移的
图谱特征。不用 network。

Scores (per feature):
  - |SMD|  standardized mean difference (W2−W1)
  - domain-AUC  how well the feature alone predicts period (W)
  - KS statistic
  - multivariate: HGB period-classifier + permutation importance

  PYTHONPATH=. python3 scripts/tencent_gr/run_w1w2_feature_select.py \\
    --root data/tencent_subset --max-users 20000 --gap-days 30
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from time_window_feats import (  # noqa: E402
    DAY,
    feature_engineer,
    fsds_feature_columns,
    propose_two_windows,
    scan_time_range,
)

ROOT = Path(__file__).resolve().parents[2]


def _matrix(df: pd.DataFrame, cols: Sequence[str]) -> np.ndarray:
    return df.reindex(columns=list(cols)).replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(
        np.float64
    )


def univariate_period_scores(
    X1: np.ndarray, X2: np.ndarray, cols: Sequence[str]
) -> pd.DataFrame:
    """Compare each feature between W1 and W2."""
    rows = []
    n1, n2 = len(X1), len(X2)
    for j, c in enumerate(cols):
        a, b = X1[:, j], X2[:, j]
        m1, m2 = float(np.mean(a)), float(np.mean(b))
        s1, s2 = float(np.std(a)), float(np.std(b))
        pooled = np.sqrt(0.5 * (s1 * s1 + s2 * s2)) + 1e-12
        smd = (m2 - m1) / pooled
        # domain AUC: can this feature alone tell W2 from W1?
        y = np.concatenate([np.zeros(n1), np.ones(n2)])
        x = np.concatenate([a, b])
        try:
            auc = float(roc_auc_score(y, x))
        except ValueError:
            auc = 0.5
        # fold to [0.5, 1]
        dom_auc = max(auc, 1.0 - auc)
        try:
            ks = float(ks_2samp(a, b, method="asymp").statistic)
        except TypeError:
            ks = float(ks_2samp(a, b).statistic)
        rows.append(
            {
                "feature": c,
                "mean_W1": m1,
                "mean_W2": m2,
                "std_W1": s1,
                "std_W2": s2,
                "smd_W2_minus_W1": smd,
                "abs_smd": abs(smd),
                "domain_auc": dom_auc,
                "ks": ks,
            }
        )
    return pd.DataFrame(rows)


def multivariate_period_importance(
    X1: np.ndarray,
    X2: np.ndarray,
    cols: Sequence[str],
    *,
    seed: int,
    max_rows: int,
) -> pd.DataFrame:
    """HGB classifies period W; permutation importance = shift drivers."""
    rng = np.random.default_rng(seed)
    n1 = min(len(X1), max_rows)
    n2 = min(len(X2), max_rows)
    i1 = rng.choice(len(X1), size=n1, replace=False)
    i2 = rng.choice(len(X2), size=n2, replace=False)
    X = np.vstack([X1[i1], X2[i2]])
    W = np.concatenate([np.zeros(n1, dtype=int), np.ones(n2, dtype=int)])
    Xtr, Xte, Wtr, Wte = train_test_split(
        X, W, test_size=0.25, random_state=seed, stratify=W
    )
    clf = HistGradientBoostingClassifier(
        max_depth=5, learning_rate=0.08, max_iter=100, random_state=seed
    )
    t0 = time.time()
    clf.fit(Xtr, Wtr)
    proba = clf.predict_proba(Xte)[:, 1]
    auc = float(roc_auc_score(Wte, proba))
    pi = permutation_importance(
        clf, Xte, Wte, n_repeats=8, random_state=seed, scoring="roc_auc"
    )
    out = pd.DataFrame(
        {
            "feature": list(cols),
            "perm_importance": pi.importances_mean,
            "perm_importance_std": pi.importances_std,
        }
    )
    out.attrs["period_auc"] = auc
    out.attrs["sec"] = float(time.time() - t0)
    out.attrs["n_train"] = int(len(Wtr))
    out.attrs["n_test"] = int(len(Wte))
    return out


def w1w2_candidate_columns(cols: Sequence[str]) -> List[str]:
    """Drop within-window ranks — not comparable across W1/W2 for landing FS."""
    return [c for c in cols if not c.endswith("_rank") and "rank" not in c]


def select_features(uni: pd.DataFrame, multi: pd.DataFrame, *, top_k: int) -> pd.DataFrame:
    """Combine scores → rank for landing feature selection."""
    m = uni.merge(multi, on="feature", how="left")
    m["perm_importance"] = m["perm_importance"].fillna(0.0)
    # rank-average of abs_smd, domain_auc, ks, perm_importance
    for col in ("abs_smd", "domain_auc", "ks", "perm_importance"):
        m[f"r_{col}"] = m[col].rank(ascending=False, method="average")
    m["rank_score"] = m[[f"r_{c}" for c in ("abs_smd", "domain_auc", "ks", "perm_importance")]].mean(
        axis=1
    )
    m = m.sort_values("rank_score").reset_index(drop=True)
    m["selected"] = (m.index < top_k).astype(int)
    m["rank"] = np.arange(1, len(m) + 1)
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=20000)
    ap.add_argument("--window-days", type=int, default=45)
    ap.add_argument("--gap-days", type=int, default=30)
    ap.add_argument("--localize-k", type=int, default=300, help="0 = use all edges")
    ap.add_argument("--top-k", type=int, default=20)
    ap.add_argument("--co-window", type=int, default=8)
    ap.add_argument("--max-rows-multi", type=int, default=25000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_w1w2_fs",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("scan timeline ...", flush=True)
    t_min, t_max, _ = scan_time_range(args.root, max_users=args.max_users)
    w1, w2 = propose_two_windows(
        t_min, t_max, window_days=args.window_days, gap_days=args.gap_days
    )
    gap_d = (w2.t_start - w1.t_end) / DAY
    print(
        f"span={(t_max-t_min)/DAY:.1f}d  W1={w1.n_days:.0f}d  gap={gap_d:.0f}d  W2={w2.n_days:.0f}d",
        flush=True,
    )
    if gap_d < args.gap_days - 1e-6:
        raise SystemExit(f"gap {gap_d:.2f}d < {args.gap_days}d")

    print("FE W1 ...", flush=True)
    p1 = feature_engineer(
        args.root,
        w1.t_start,
        w1.t_end,
        max_users=args.max_users,
        co_window=args.co_window,
        window_name=w1.name,
        terminal_action=1,
    )
    print("FE W2 ...", flush=True)
    p2 = feature_engineer(
        args.root,
        w2.t_start,
        w2.t_end,
        max_users=args.max_users,
        co_window=args.co_window,
        window_name=w2.name,
        terminal_action=1,
    )
    g1, g2 = p1["grid"], p2["grid"]
    print(
        f"W1 edges={len(g1)} pos={p1['meta']['pos_rate']:.4f} | "
        f"W2 edges={len(g2)} pos={p2['meta']['pos_rate']:.4f}",
        flush=True,
    )

    # optional localization subset from W1 shares (landing: focus on attributed items)
    if args.localize_k > 0 and len(p1["item_df"]):
        loc = (
            p1["item_df"]
            .sort_values("share_linear", ascending=False)
            .head(args.localize_k)["item_id"]
            .astype(int)
            .tolist()
        )
        loc_set = set(loc)
        g1 = g1[g1["item_id"].isin(loc_set)].copy()
        g2 = g2[g2["item_id"].isin(loc_set)].copy()
        print(f"localize subset k={len(loc)} → W1={len(g1)} W2={len(g2)}", flush=True)
    else:
        loc = []

    raw_cols = [c for c in fsds_feature_columns(g1) if c in g2.columns]
    cols = w1w2_candidate_columns(raw_cols)
    dropped = sorted(set(raw_cols) - set(cols))
    if dropped:
        print(f"drop within-window ranks: {', '.join(dropped)}", flush=True)
    print(f"candidate features n={len(cols)}", flush=True)
    if not cols:
        raise SystemExit("no feature columns")
    X1, X2 = _matrix(g1, cols), _matrix(g2, cols)

    print("univariate W1 vs W2 ...", flush=True)
    uni = univariate_period_scores(X1, X2, cols)
    print("multivariate period classifier ...", flush=True)
    multi = multivariate_period_importance(
        X1, X2, cols, seed=args.seed, max_rows=args.max_rows_multi
    )
    ranked = select_features(uni, multi, top_k=args.top_k)
    selected = ranked[ranked["selected"] == 1]["feature"].tolist()

    # outcome snapshot (not used for selection — just report)
    y_snap = {
        "W1_pos_rate": float(g1["y_convert"].mean()) if len(g1) else float("nan"),
        "W2_pos_rate": float(g2["y_convert"].mean()) if len(g2) else float("nan"),
        "period_clf_auc": float(multi.attrs.get("period_auc", float("nan"))),
    }

    ranked.to_csv(args.out_dir / "w1w2_feature_ranking.csv", index=False)
    pd.DataFrame({"feature": selected, "rank": np.arange(1, len(selected) + 1)}).to_csv(
        args.out_dir / "selected_shift_features.csv", index=False
    )
    if loc:
        pd.DataFrame({"item_id": loc, "rank": np.arange(1, len(loc) + 1)}).to_csv(
            args.out_dir / "localize_subset_W1.csv", index=False
        )

    blob = {
        "protocol": [
            "FE(time_range) independently on W1 and W2 (gap >= 30d)",
            "optional localization subset from W1 share_linear",
            "COMPARE W1 vs W2 → feature selection (SMD / domain-AUC / KS / perm-imp)",
            "selected features = temporal shift drivers (落地)",
        ],
        "timeline": {
            "t_min": t_min,
            "t_max": t_max,
            "span_days": (t_max - t_min) / DAY,
            "W1": w1.to_dict(),
            "W2": w2.to_dict(),
            "gap_days": gap_d,
        },
        "W1_meta": p1["meta"],
        "W2_meta": p2["meta"],
        "n_W1": int(len(g1)),
        "n_W2": int(len(g2)),
        "n_features": len(cols),
        "top_k": args.top_k,
        "selected": selected,
        "y_snapshot": y_snap,
        "multi_sec": multi.attrs.get("sec"),
    }
    (args.out_dir / "summary.json").write_text(json.dumps(blob, indent=2, default=str))

    # plot
    top = ranked.head(min(15, len(ranked)))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.barh(range(len(top)), top["abs_smd"][::-1], color="#4C78A8", label="|SMD|")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["feature"][::-1], fontsize=8)
    ax.set_xlabel("|SMD| (W2 vs W1)")
    ax.set_title("Univariate shift")
    ax = axes[1]
    ax.barh(range(len(top)), top["perm_importance"][::-1], color="#F58518")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["feature"][::-1], fontsize=8)
    ax.set_xlabel("perm importance (period clf)")
    ax.set_title(f"Multivariate (period AUC={y_snap['period_clf_auc']:.3f})")
    fig.suptitle(
        f"W1 vs W2 feature selection  gap={gap_d:.0f}d  top-{args.top_k}",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(args.out_dir / "w1w2_feature_select.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    md = [
        "# W1 vs W2 → feature selection (落地)",
        "",
        "比较两个时间窗的图谱特征，**用差异做 feature selection**。",
        "选出的是时段漂移驱动特征，可直接进监控 / 重训 / 归因。",
        "",
        "## Protocol",
        f"1. `feature_engineer(t_start,t_end)` ×2，gap = **{gap_d:.0f}** days (≥{args.gap_days})",
        f"2. Localization subset（可选）: W1 `share_linear` top-{args.localize_k}",
        "3. **Compare W1 vs W2**: |SMD| + domain-AUC + KS + period-HGB perm-importance",
        f"4. Rank-average → select top-**{args.top_k}** (drop within-window `*_rank`)",
        "",
        "## Windows",
        f"- span **{(t_max-t_min)/DAY:.1f}**d | W1 n={len(g1)} | W2 n={len(g2)} | feats={len(cols)}",
        f"- period classifier AUC = **{y_snap['period_clf_auc']:.3f}** "
        f"(how separable are the two windows)",
        f"- click rate snapshot: W1={y_snap['W1_pos_rate']:.4f} → W2={y_snap['W2_pos_rate']:.4f}",
        "",
        "## Selected shift drivers",
        "| rank | feature | \\|SMD\\| | domain-AUC | KS | perm-imp |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for _, r in ranked[ranked["selected"] == 1].iterrows():
        md.append(
            f"| {int(r['rank'])} | `{r['feature']}` | {r['abs_smd']:.3f} | "
            f"{r['domain_auc']:.3f} | {r['ks']:.3f} | {r['perm_importance']:.4f} |"
        )
    md += [
        "",
        "## Why this lands",
        "- 输入是业务已有的 `(user,item)` 图谱特征 + 时间窗",
        "- 输出是 **W1→W2 漂移特征清单**，可挂告警 / 重训触发 / 运营解释",
        "- 无图算法；防泄漏：两窗独立 FE，选择信号是 period W 而非 peek 未来标签训练",
        "",
    ]
    (args.out_dir / "W1W2_FEATURE_SELECT_REPORT.md").write_text("\n".join(md))
    print("wrote", args.out_dir)
    print("period AUC", y_snap["period_clf_auc"])
    print("selected:", ", ".join(selected[:12]), ("..." if len(selected) > 12 else ""))


if __name__ == "__main__":
    main()
