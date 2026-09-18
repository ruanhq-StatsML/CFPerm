#!/usr/bin/env python3
"""Streamlined: Standardization → subset-level MMD → FSDS.

Landing procedure (nothing else):
  0. StandardScaler fit on W1
  1. FE independently on W1 / W2 (gap ≥ 30d)
  2. Subset localization = item-level RBF-MMD² (standardized X)
  3. FSDS = StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg
     → feature ranking

  PYTHONPATH=. python3 scripts/tencent_gr/run_standardize_mmd_fsds.py \\
    --root data/tencent_subset --max-users 20000 --gap-days 30 --localize-k 200
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

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
from gt_subset_evaluator import evaluate_gt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def _matrix(df: pd.DataFrame, cols: Sequence[str]) -> np.ndarray:
    return (
        df.reindex(columns=list(cols))
        .replace([np.inf, -np.inf], np.nan)
        .fillna(0.0)
        .to_numpy(np.float64)
    )


def fit_standardizer(X: np.ndarray) -> StandardScaler:
    sc = StandardScaler(with_mean=True, with_std=True)
    sc.fit(np.asarray(X, np.float64))
    scale = np.asarray(sc.scale_, dtype=np.float64)
    scale[~np.isfinite(scale) | (scale < 1e-8)] = 1.0
    sc.scale_ = scale
    return sc


def standardize(sc: StandardScaler, X: np.ndarray) -> np.ndarray:
    return sc.transform(np.asarray(X, np.float64))


def w1w2_candidate_columns(cols: Sequence[str]) -> List[str]:
    return [c for c in cols if not c.endswith("_rank") and "rank" not in c]


def rbf_mmd2(
    X0: np.ndarray,
    X1: np.ndarray,
    *,
    max_n: int = 256,
    rng: Optional[np.random.Generator] = None,
) -> float:
    """Unbiased RBF MMD² (median bandwidth). Expects standardized X."""
    rng = rng or np.random.default_rng(0)
    X0 = np.asarray(X0, np.float64)
    X1 = np.asarray(X1, np.float64)
    if len(X0) < 2 or len(X1) < 2:
        return 0.0
    if len(X0) > max_n:
        X0 = X0[rng.choice(len(X0), max_n, replace=False)]
    if len(X1) > max_n:
        X1 = X1[rng.choice(len(X1), max_n, replace=False)]
    Z = np.vstack([X0, X1])
    Zs = Z if len(Z) <= 400 else Z[rng.choice(len(Z), 400, replace=False)]
    d2 = np.sum((Zs[:, None, :] - Zs[None, :, :]) ** 2, axis=-1)
    med = float(np.median(d2[d2 > 0])) if np.any(d2 > 0) else 1.0
    gamma = 1.0 / (2.0 * med + 1e-12)

    def k(A, B):
        return np.exp(-gamma * np.sum((A[:, None, :] - B[None, :, :]) ** 2, axis=-1))

    Kxx, Kyy, Kxy = k(X0, X0), k(X1, X1), k(X0, X1)
    n, m = len(X0), len(X1)
    np.fill_diagonal(Kxx, 0.0)
    np.fill_diagonal(Kyy, 0.0)
    return float(
        Kxx.sum() / (n * (n - 1) + 1e-12)
        + Kyy.sum() / (m * (m - 1) + 1e-12)
        - 2.0 * Kxy.mean()
    )


def score_item_mmd(
    g1: pd.DataFrame,
    g2: pd.DataFrame,
    cols: Sequence[str],
    sc: StandardScaler,
    *,
    seed: int,
    max_cand: int,
    mmd_max_n: int,
    min_edges: int,
) -> pd.DataFrame:
    """Subset localization: rank items by RBF-MMD²(W1, W2) on standardized X."""
    c1 = g1["item_id"].value_counts()
    c2 = g2["item_id"].value_counts()
    common = sorted(set(c1[c1 >= min_edges].index) & set(c2[c2 >= min_edges].index))
    vol = {i: int(c1.get(i, 0) + c2.get(i, 0)) for i in common}
    cand = sorted(common, key=lambda i: -vol[i])[:max_cand]
    print(f"item candidates (≥{min_edges} edges both windows): {len(cand)}", flush=True)

    rng = np.random.default_rng(seed)
    rows = []
    for iid in cand:
        a = standardize(sc, _matrix(g1[g1["item_id"] == iid], cols))
        b = standardize(sc, _matrix(g2[g2["item_id"] == iid], cols))
        rows.append(
            {
                "item_id": int(iid),
                "n_W1": int(len(a)),
                "n_W2": int(len(b)),
                "mmd2": rbf_mmd2(a, b, max_n=mmd_max_n, rng=rng),
            }
        )
    scor = pd.DataFrame(rows).sort_values("mmd2", ascending=False).reset_index(drop=True)
    scor["rank"] = np.arange(1, len(scor) + 1)
    return scor


def run_fsds(
    grid_tr: pd.DataFrame,
    grid_te: pd.DataFrame,
    cols: List[str],
    *,
    select_k: int,
    seed: int,
) -> Dict:
    """FSDS: StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg."""
    Xtr = _matrix(grid_tr, cols)
    ytr = grid_tr["y_convert"].to_numpy(int)
    Xte = _matrix(grid_te, cols)
    yte = grid_te["y_convert"].to_numpy(int)
    if len(np.unique(ytr)) < 2:
        return {"ok": False, "reason": "train needs both classes"}

    k = min(select_k, Xtr.shape[1], max(1, Xtr.shape[0] - 1))
    pipe = Pipeline(
        [
            ("sc", StandardScaler(with_mean=True, with_std=True)),
            ("var", VarianceThreshold(1e-8)),
            ("sel", SelectKBest(f_classif, k=k)),
        ]
    )
    t0 = time.time()
    Xt = pipe.fit_transform(Xtr, ytr)
    var_mask = pipe.named_steps["var"].get_support()
    cols_var = [c for c, m in zip(cols, var_mask) if m]
    sel = pipe.named_steps["sel"]
    selected = [c for c, m in zip(cols_var, sel.get_support()) if m]
    ranking = (
        pd.DataFrame({"feature": cols_var, "f_score": sel.scores_})
        .sort_values("f_score", ascending=False)
        .reset_index(drop=True)
    )
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    ranking["selected"] = ranking["feature"].isin(selected).astype(int)

    out: Dict = {
        "ok": True,
        "n_train": int(len(ytr)),
        "n_test": int(len(yte)),
        "pos_train": float(ytr.mean()),
        "pos_test": float(yte.mean()) if len(yte) else float("nan"),
        "n_in": len(cols),
        "n_selected": len(selected),
        "selected": selected,
        "ranking": ranking,
        "sec_select": float(time.time() - t0),
        "models": {},
        "pipeline": "StandardScaler → VarianceThreshold → SelectKBest → model",
    }
    if len(yte) == 0 or len(np.unique(yte)) < 2:
        out["models"]["note"] = "test missing both classes — ranking only"
        return out

    Xv = pipe.transform(Xte)
    hgb = HistGradientBoostingClassifier(
        max_depth=6, learning_rate=0.08, max_iter=120, random_state=seed, class_weight="balanced"
    )
    t0 = time.time()
    hgb.fit(Xt, ytr)
    ph = hgb.predict_proba(Xv)[:, 1]
    out["models"]["hgb"] = {
        "auc": float(roc_auc_score(yte, ph)),
        "ap": float(average_precision_score(yte, ph)),
        "sec": float(time.time() - t0),
    }
    lr = LogisticRegression(max_iter=400, C=0.5, class_weight="balanced", random_state=seed)
    t0 = time.time()
    lr.fit(Xt, ytr)
    pl = lr.predict_proba(Xv)[:, 1]
    out["models"]["logreg"] = {
        "auc": float(roc_auc_score(yte, pl)),
        "ap": float(average_precision_score(yte, pl)),
        "sec": float(time.time() - t0),
    }
    return out


def plot_results(subset: pd.DataFrame, feat_rank: pd.DataFrame, *, out_path: Path, title: str) -> None:
    top_i = subset.head(15)
    top_f = feat_rank.head(15)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ax = axes[0]
    ax.barh(range(len(top_i)), top_i["mmd2"][::-1], color="#F58518")
    ax.set_yticks(range(len(top_i)))
    ax.set_yticklabels([str(x) for x in top_i["item_id"][::-1]], fontsize=7)
    ax.set_xlabel("MMD² (standardized)")
    ax.set_title("Subset localization · MMD")

    ax = axes[1]
    ax.barh(range(len(top_f)), top_f["f_score"][::-1], color="#54A24B")
    ax.set_yticks(range(len(top_f)))
    ax.set_yticklabels(top_f["feature"][::-1], fontsize=8)
    ax.set_xlabel("FSDS F-score")
    ax.set_title("Feature ranking · FSDS")

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Standardization → subset-MMD → FSDS")
    ap.add_argument("--root", type=Path, default=ROOT / "data" / "tencent_subset")
    ap.add_argument("--max-users", type=int, default=20000)
    ap.add_argument("--window-days", type=int, default=45)
    ap.add_argument("--gap-days", type=int, default=30)
    ap.add_argument("--localize-k", type=int, default=200)
    ap.add_argument("--max-cand", type=int, default=400)
    ap.add_argument("--min-edges", type=int, default=3)
    ap.add_argument("--mmd-max-n", type=int, default=128)
    ap.add_argument("--select-k", type=int, default=15)
    ap.add_argument("--co-window", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gt-items", type=Path, default=None)
    ap.add_argument("--gt-orders", type=Path, default=None)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_standardize_mmd_fsds",
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
        args.root, w1.t_start, w1.t_end,
        max_users=args.max_users, co_window=args.co_window,
        window_name=w1.name, terminal_action=1,
    )
    print("FE W2 ...", flush=True)
    p2 = feature_engineer(
        args.root, w2.t_start, w2.t_end,
        max_users=args.max_users, co_window=args.co_window,
        window_name=w2.name, terminal_action=1,
    )
    g1, g2 = p1["grid"], p2["grid"]
    print(
        f"W1 edges={len(g1)} pos={p1['meta']['pos_rate']:.4f} | "
        f"W2 edges={len(g2)} pos={p2['meta']['pos_rate']:.4f}",
        flush=True,
    )

    raw_cols = [c for c in fsds_feature_columns(g1) if c in g2.columns]
    cols = w1w2_candidate_columns(raw_cols)
    print(f"features n={len(cols)}", flush=True)

    # 1) Standardization
    print("1/3 Standardization (fit on W1) ...", flush=True)
    sc = fit_standardizer(_matrix(g1, cols))

    # 2) Subset-level MMD
    print("2/3 subset-level MMD ...", flush=True)
    scor = score_item_mmd(
        g1, g2, cols, sc,
        seed=args.seed, max_cand=args.max_cand,
        mmd_max_n=args.mmd_max_n, min_edges=args.min_edges,
    )
    k_loc = min(args.localize_k, len(scor))
    scor["selected"] = (scor.index < k_loc).astype(int)
    subset = scor[scor["selected"] == 1].copy()
    loc_items = subset["item_id"].astype(int).tolist()
    loc_set = set(loc_items)
    print(f"localized items k={len(loc_items)}", flush=True)

    gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=args.seed)
    tr_idx, va_idx = next(gss.split(g1, groups=g1["user_id"].to_numpy()))
    g1_tr = g1.iloc[tr_idx]
    g1_va = g1.iloc[va_idx]
    g1_tr = g1_tr[g1_tr["item_id"].isin(loc_set)].copy()
    g1_va = g1_va[g1_va["item_id"].isin(loc_set)].copy()
    g2_loc = g2[g2["item_id"].isin(loc_set)].copy()
    if len(g1_tr) < 40:
        g1_tr = g1[g1["item_id"].isin(loc_set)].copy()
        g1_va = g1_tr.sample(frac=0.25, random_state=args.seed) if len(g1_tr) else g1_tr
    print(
        f"localized edges W1-train={len(g1_tr)} W1-holdout={len(g1_va)} W2={len(g2_loc)}",
        flush=True,
    )

    # 3) FSDS
    print("3/3 FSDS ranking ...", flush=True)
    res_w1 = run_fsds(g1_tr, g1_va, cols, select_k=args.select_k, seed=args.seed)
    res_w2 = run_fsds(g1_tr, g2_loc, cols, select_k=args.select_k, seed=args.seed)
    ranking = res_w1.get("ranking")
    if not isinstance(ranking, pd.DataFrame) or ranking.empty:
        ranking = pd.DataFrame({"feature": cols, "f_score": 0.0, "rank": range(1, len(cols) + 1), "selected": 0})

    scor.to_csv(args.out_dir / "item_mmd_scores.csv", index=False)
    subset.to_csv(args.out_dir / "localized_subset_items.csv", index=False)
    ranking.to_csv(args.out_dir / "fsds_feature_ranking.csv", index=False)

    gt_eval = evaluate_gt(
        loc_items,
        gt_items_path=args.gt_items,
        gt_orders_path=args.gt_orders,
        ks=(50, 100, min(200, len(loc_items)) or 1),
    )
    if gt_eval.get("available"):
        (args.out_dir / "gt_eval.json").write_text(json.dumps(gt_eval, indent=2))
        print("GT eval:", gt_eval.get("items"), flush=True)

    def _strip(res: Dict) -> Dict:
        out = {k: v for k, v in res.items() if k != "ranking"}
        if isinstance(res.get("ranking"), pd.DataFrame):
            out["top_features"] = res["ranking"].head(args.select_k)["feature"].tolist()
        return out

    blob = {
        "procedure": ["Standardization", "subset-level MMD", "FSDS"],
        "timeline": {
            "span_days": (t_max - t_min) / DAY,
            "W1": w1.to_dict(),
            "W2": w2.to_dict(),
            "gap_days": gap_d,
        },
        "W1_meta": p1["meta"],
        "W2_meta": p2["meta"],
        "n_features": len(cols),
        "localize_k": len(loc_items),
        "localized_items_head": loc_items[:30],
        "n_localized_edges": {
            "W1_train": int(len(g1_tr)),
            "W1_holdout": int(len(g1_va)),
            "W2": int(len(g2_loc)),
        },
        "fsds_W1_holdout": _strip(res_w1),
        "fsds_W2_temporal": _strip(res_w2),
        "gt_eval": gt_eval,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(blob, indent=2, default=str))

    plot_results(
        subset,
        ranking,
        out_path=args.out_dir / "standardize_mmd_fsds.png",
        title=f"Standardize → subset-MMD → FSDS  gap={gap_d:.0f}d  k={len(loc_items)}",
    )

    md = [
        "# Streamlined: Standardization → subset-MMD → FSDS",
        "",
        "```",
        "StandardScaler(W1) → item MMD²(W1,W2) → FSDS ranking",
        "```",
        "",
        "## Protocol",
        "1. **Standardization** fit on W1",
        f"2. **Subset-level MMD** → top-**{len(loc_items)}** items (gap={gap_d:.0f}d)",
        f"3. **FSDS**: StandardScaler → var → SelectKBest(k={args.select_k}) → HGB/LogReg",
        "",
        "## Localized subset (by MMD²)",
        "| rank | item_id | MMD² | n_W1 | n_W2 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for _, r in subset.head(12).iterrows():
        md.append(
            f"| {int(r['rank'])} | {int(r['item_id'])} | {r['mmd2']:.4f} | "
            f"{int(r['n_W1'])} | {int(r['n_W2'])} |"
        )
    md += [
        "",
        "## FSDS feature ranking",
        "| rank | feature | F-score | selected |",
        "|---:|---|---:|---:|",
    ]
    for _, r in ranking.head(args.select_k).iterrows():
        md.append(
            f"| {int(r['rank'])} | `{r['feature']}` | {float(r['f_score']):.3f} | "
            f"{int(r['selected'])} |"
        )

    def _auc(tag: str, res: Dict) -> str:
        if not res.get("ok"):
            return f"- {tag}: fail ({res.get('reason')})"
        mods = res.get("models") or {}
        bits = []
        for name in ("hgb", "logreg"):
            if name in mods and "auc" in mods[name]:
                bits.append(f"{name} AUC={mods[name]['auc']:.3f}")
        return f"- {tag}: n={res.get('n_train')}/{res.get('n_test')} | " + (
            " | ".join(bits) if bits else "ranking only"
        )

    md += ["", "## Holdout", _auc("W1 user-holdout", res_w1), _auc("W2 temporal", res_w2), ""]
    if gt_eval.get("available"):
        it = gt_eval.get("items") or {}
        md.append(
            f"## GT\n- P@100={it.get('precision@100', float('nan')):.3f} | "
            f"R@100={it.get('recall@100', float('nan')):.3f}\n"
        )
    (args.out_dir / "STANDARDIZE_MMD_FSDS_REPORT.md").write_text("\n".join(md))
    print("wrote", args.out_dir)
    print("top features:", ", ".join(ranking.head(8)["feature"].astype(str).tolist()))


if __name__ == "__main__":
    main()
