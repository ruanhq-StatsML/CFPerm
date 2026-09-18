#!/usr/bin/env python3
"""W1 vs W2: MMD + PO-risk + conditional-mean subset localization → FSDS ranking.

落地流程（无 network；FSDS 标准流 + 换 dataset）
----------------------------------------------
0. **Standardization**（最上面）：W1 上 fit ``StandardScaler``，全程共用
1. ``feature_engineer`` 独立跑 W1 / W2（gap ≥ 30d）
2. **Subset localization**（item 级）用三个直观分数：
   - conditional-mean ‖μ_W2 − μ_W1‖（standardized space）
   - RBF-MMD²(X|item,W1 ; X|item,W2)
   - PO-risk 聚集：全局 τ̂ 在该 item 上的 mean(τ̂²)
3. Rank-average → top-k subset，可视化
4. **FSDS** 标准流：StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg
   给出特征 ranking（W1-train fit；W2 temporal holdout）

  PYTHONPATH=. python3 scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py \\
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
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestRegressor
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit, StratifiedKFold
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
    """Top-of-pipeline standardization (fit once on W1)."""
    sc = StandardScaler(with_mean=True, with_std=True)
    sc.fit(np.asarray(X, np.float64))
    # avoid zero-scale blowups on near-constant cols
    scale = np.asarray(sc.scale_, dtype=np.float64)
    scale[~np.isfinite(scale) | (scale < 1e-8)] = 1.0
    sc.scale_ = scale
    return sc


def standardize(sc: StandardScaler, X: np.ndarray) -> np.ndarray:
    return sc.transform(np.asarray(X, np.float64))


def w1w2_candidate_columns(cols: Sequence[str]) -> List[str]:
    """Drop within-window ranks — not comparable across W1/W2."""
    return [c for c in cols if not c.endswith("_rank") and "rank" not in c]


def rbf_mmd2(
    X0: np.ndarray,
    X1: np.ndarray,
    *,
    max_n: int = 256,
    rng: Optional[np.random.Generator] = None,
) -> float:
    """Unbiased RBF MMD² (median bandwidth). Expects already-standardized X."""
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


def conditional_mean_l2(X0: np.ndarray, X1: np.ndarray) -> float:
    """‖μ1 − μ0‖₂ in standardized feature space."""
    if len(X0) == 0 or len(X1) == 0:
        return 0.0
    return float(np.linalg.norm(X1.mean(axis=0) - X0.mean(axis=0)))


def po_risk_fit(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int,
    n_trees: int = 40,
) -> Dict:
    """DR-style PO → τ̂(X); risk = mean(τ̂²). Matches repo PO-risk usage."""
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    Yf = Y.astype(float)
    n_splits = 3 if n >= 60 else 2
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=n_trees,
            max_depth=6,
            min_samples_leaf=4,
            random_state=seed + fold,
            n_jobs=1,
        )
        from sklearn.ensemble import RandomForestClassifier

        e = RandomForestClassifier(
            n_estimators=n_trees,
            max_depth=6,
            min_samples_leaf=4,
            random_state=seed + 40 + fold,
            n_jobs=1,
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, 0.05, 0.95)
    po = (Yf - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=max(60, n_trees),
        max_depth=8,
        min_samples_leaf=4,
        random_state=seed + 7,
        n_jobs=1,
    )
    tau.fit(X, po)
    tau_hat = tau.predict(X)
    return {
        "risk": float(np.mean(tau_hat**2)),
        "tau_hat": tau_hat,
        "po": po,
        "vimp": tau.feature_importances_.astype(float),
    }


def feature_conditional_mean(X0: np.ndarray, X1: np.ndarray, cols: Sequence[str]) -> pd.DataFrame:
    """Per-feature |μ_W2 − μ_W1| in standardized space."""
    rows = []
    for j, c in enumerate(cols):
        a, b = X0[:, j], X1[:, j]
        rows.append(
            {
                "feature": c,
                "cmean_abs": abs(float(b.mean() - a.mean())),
                "mean_W1": float(a.mean()),
                "mean_W2": float(b.mean()),
            }
        )
    return pd.DataFrame(rows)


def feature_mmd_loco(
    X0: np.ndarray,
    X1: np.ndarray,
    cols: Sequence[str],
    *,
    seed: int,
    max_n: int = 200,
) -> pd.DataFrame:
    """Leave-one-covariate-out ΔMMD as feature contribution (standardized X)."""
    rng = np.random.default_rng(seed)
    full = rbf_mmd2(X0, X1, max_n=max_n, rng=rng)
    rows = []
    for j, c in enumerate(cols):
        keep = [i for i in range(X0.shape[1]) if i != j]
        if not keep:
            delta = 0.0
        else:
            d = rbf_mmd2(X0[:, keep], X1[:, keep], max_n=max_n, rng=rng)
            delta = full - d
        rows.append({"feature": c, "mmd_loco": float(delta), "mmd_full": full})
    return pd.DataFrame(rows)


def score_item_subset(
    g1: pd.DataFrame,
    g2: pd.DataFrame,
    cols: Sequence[str],
    sc: StandardScaler,
    *,
    seed: int,
    max_cand: int,
    mmd_max_n: int,
    min_edges: int,
) -> Tuple[pd.DataFrame, Dict]:
    """Item-level subset localization via cmean / MMD / PO-risk (on standardized X)."""
    c1 = g1["item_id"].value_counts()
    c2 = g2["item_id"].value_counts()
    common = sorted(set(c1[c1 >= min_edges].index) & set(c2[c2 >= min_edges].index))
    vol = {i: int(c1.get(i, 0) + c2.get(i, 0)) for i in common}
    cand = sorted(common, key=lambda i: -vol[i])[:max_cand]
    print(f"item candidates (both windows, ≥{min_edges} edges): {len(cand)}", flush=True)

    rng = np.random.default_rng(seed)
    n_po = min(8000, len(g1) + len(g2))
    i1 = rng.choice(len(g1), size=min(len(g1), n_po // 2), replace=False)
    i2 = rng.choice(len(g2), size=min(len(g2), n_po // 2), replace=False)
    Xp = standardize(
        sc, np.vstack([_matrix(g1.iloc[i1], cols), _matrix(g2.iloc[i2], cols)])
    )
    Wp = np.concatenate([np.zeros(len(i1), dtype=int), np.ones(len(i2), dtype=int)])
    # continuous Y for stable PO when clicks are rare (already on standardized X)
    Yp = PCA(n_components=1, random_state=seed).fit_transform(Xp).ravel()
    order = np.argsort(Yp, kind="mergesort")
    ranks = np.empty(len(Yp), float)
    ranks[order] = np.linspace(0.0, 1.0, len(Yp))
    print("fit global PO-risk ...", flush=True)
    t0 = time.time()
    po_fit = po_risk_fit(Xp, ranks, Wp, seed=seed)
    print(f"  PO-risk={po_fit['risk']:.6f}  sec={time.time()-t0:.1f}", flush=True)

    items_po = np.concatenate(
        [g1.iloc[i1]["item_id"].to_numpy(), g2.iloc[i2]["item_id"].to_numpy()]
    )
    tau2 = po_fit["tau_hat"] ** 2

    rows = []
    for iid in cand:
        a = standardize(sc, _matrix(g1[g1["item_id"] == iid], cols))
        b = standardize(sc, _matrix(g2[g2["item_id"] == iid], cols))
        cm = conditional_mean_l2(a, b)
        mm = rbf_mmd2(a, b, max_n=mmd_max_n, rng=rng)
        mask = items_po == iid
        po_s = float(tau2[mask].mean()) if mask.any() else 0.0
        rows.append(
            {
                "item_id": int(iid),
                "n_W1": int(len(a)),
                "n_W2": int(len(b)),
                "cmean_l2": cm,
                "mmd2": mm,
                "po_tau2_mean": po_s,
            }
        )
    scor = pd.DataFrame(rows)
    for col in ("cmean_l2", "mmd2", "po_tau2_mean"):
        scor[f"r_{col}"] = scor[col].rank(ascending=False, method="average")
    scor["rank_score"] = scor[["r_cmean_l2", "r_mmd2", "r_po_tau2_mean"]].mean(axis=1)
    scor = scor.sort_values("rank_score").reset_index(drop=True)
    scor["rank"] = np.arange(1, len(scor) + 1)
    meta = {
        "po_risk_global": po_fit["risk"],
        "po_vimp": {c: float(v) for c, v in zip(cols, po_fit["vimp"])},
        "n_cand": len(cand),
        "n_po_rows": int(len(Wp)),
        "standardized": True,
    }
    return scor, meta


def run_fsds(
    grid_tr: pd.DataFrame,
    grid_te: pd.DataFrame,
    cols: List[str],
    *,
    select_k: int,
    seed: int,
) -> Dict:
    """Standard FSDS: StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg."""
    Xtr = _matrix(grid_tr, cols)
    ytr = grid_tr["y_convert"].to_numpy(int)
    Xte = _matrix(grid_te, cols)
    yte = grid_te["y_convert"].to_numpy(int)
    if len(np.unique(ytr)) < 2:
        return {"ok": False, "reason": "train needs both classes"}

    k = min(select_k, Xtr.shape[1], max(1, Xtr.shape[0] - 1))
    # standardization FIRST — then var filter + SelectKBest (FSDS standard)
    pipe_prep = Pipeline(
        [
            ("sc", StandardScaler(with_mean=True, with_std=True)),
            ("var", VarianceThreshold(1e-8)),
            ("sel", SelectKBest(f_classif, k=k)),
        ]
    )
    t0 = time.time()
    Xt = pipe_prep.fit_transform(Xtr, ytr)
    var_mask = pipe_prep.named_steps["var"].get_support()
    cols_var = [c for c, m in zip(cols, var_mask) if m]
    sel = pipe_prep.named_steps["sel"]
    sel_mask = sel.get_support()
    selected = [c for c, m in zip(cols_var, sel_mask) if m]
    scores = sel.scores_
    ranking = (
        pd.DataFrame({"feature": cols_var, "f_score": scores})
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

    Xv = pipe_prep.transform(Xte)
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
    # already standardized upstream — LogReg without nested scaler
    lr = LogisticRegression(
        max_iter=400, C=0.5, class_weight="balanced", random_state=seed
    )
    t0 = time.time()
    lr.fit(Xt, ytr)
    pl = lr.predict_proba(Xv)[:, 1]
    out["models"]["logreg"] = {
        "auc": float(roc_auc_score(yte, pl)),
        "ap": float(average_precision_score(yte, pl)),
        "sec": float(time.time() - t0),
    }
    return out


def plot_results(
    subset: pd.DataFrame,
    feat_rank: pd.DataFrame,
    *,
    out_path: Path,
    title: str,
) -> None:
    top_i = subset.head(15)
    top_f = feat_rank.head(15)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    ax.barh(range(len(top_i)), top_i["cmean_l2"][::-1], color="#4C78A8")
    ax.set_yticks(range(len(top_i)))
    ax.set_yticklabels([str(x) for x in top_i["item_id"][::-1]], fontsize=7)
    ax.set_xlabel("conditional-mean ‖μ₂−μ₁‖")
    ax.set_title("Behavior-shift subset · cmean")

    ax = axes[0, 1]
    ax.barh(range(len(top_i)), top_i["mmd2"][::-1], color="#F58518")
    ax.set_yticks(range(len(top_i)))
    ax.set_yticklabels([str(x) for x in top_i["item_id"][::-1]], fontsize=7)
    ax.set_xlabel("MMD²")
    ax.set_title("Purchase-intent shift · MMD")

    ax = axes[1, 0]
    ax.scatter(subset["cmean_l2"], subset["mmd2"], c=subset["po_tau2_mean"], cmap="viridis", s=18)
    ax.set_xlabel("cmean")
    ax.set_ylabel("MMD²")
    ax.set_title("Items (color = PO-risk τ²)")

    ax = axes[1, 1]
    ax.barh(range(len(top_f)), top_f["f_score"][::-1], color="#54A24B")
    ax.set_yticks(range(len(top_f)))
    ax.set_yticklabels(top_f["feature"][::-1], fontsize=8)
    ax.set_xlabel("FSDS F-score")
    ax.set_title("Attribution feature ranking (FSDS)")

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
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
    ap.add_argument(
        "--gt-items",
        type=Path,
        default=None,
        help="GT item list CSV/parquet/txt (item_id / oid / sku_id …)",
    )
    ap.add_argument(
        "--gt-orders",
        type=Path,
        default=None,
        help="GT order×item CSV/parquet (item_id + optional order_id)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_w1w2_mmd_po_fsds",
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

    raw_cols = [c for c in fsds_feature_columns(g1) if c in g2.columns]
    cols = w1w2_candidate_columns(raw_cols)
    print(f"features n={len(cols)} (dropped ranks={len(raw_cols)-len(cols)})", flush=True)

    # --- 0. Standardization on W1 (top of pipeline; no W2 peek) ---
    print("fit StandardScaler on W1 ...", flush=True)
    sc = fit_standardizer(_matrix(g1, cols))

    # --- subset localization ---
    print("subset localization (cmean + MMD + PO) on standardized X ...", flush=True)
    scor, loc_meta = score_item_subset(
        g1,
        g2,
        cols,
        sc,
        seed=args.seed,
        max_cand=args.max_cand,
        mmd_max_n=args.mmd_max_n,
        min_edges=args.min_edges,
    )
    k_loc = min(args.localize_k, len(scor))
    scor["selected"] = (scor.index < k_loc).astype(int)
    subset = scor[scor["selected"] == 1].copy()
    loc_items = subset["item_id"].astype(int).tolist()
    loc_set = set(loc_items)
    print(f"localized items k={len(loc_items)}", flush=True)

    # user-grouped W1 split; shares never needed — subset from MMD/PO/cmean only
    gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=args.seed)
    u = g1["user_id"].to_numpy()
    tr_idx, va_idx = next(gss.split(g1, groups=u))
    g1_tr_all, g1_va_all = g1.iloc[tr_idx], g1.iloc[va_idx]
    g1_tr = g1_tr_all[g1_tr_all["item_id"].isin(loc_set)].copy()
    g1_va = g1_va_all[g1_va_all["item_id"].isin(loc_set)].copy()
    g2_loc = g2[g2["item_id"].isin(loc_set)].copy()
    print(
        f"localized edges W1-train={len(g1_tr)} W1-holdout={len(g1_va)} W2={len(g2_loc)}",
        flush=True,
    )

    # feature-level shift diagnostics on standardized localized matrices
    X1s = standardize(sc, _matrix(g1_tr if len(g1_tr) else g1, cols))
    X2s = standardize(sc, _matrix(g2_loc if len(g2_loc) else g2, cols))
    rng = np.random.default_rng(args.seed)
    n1 = min(len(X1s), 1500)
    n2 = min(len(X2s), 1500)
    X1b = X1s[rng.choice(len(X1s), n1, replace=False)] if len(X1s) else X1s
    X2b = X2s[rng.choice(len(X2s), n2, replace=False)] if len(X2s) else X2s
    feat_cmean = feature_conditional_mean(X1b, X2b, cols)
    feat_mmd = feature_mmd_loco(X1b, X2b, cols, seed=args.seed, max_n=args.mmd_max_n)
    feat_diag = feat_cmean.merge(feat_mmd, on="feature")
    po_v = loc_meta.get("po_vimp", {})
    feat_diag["po_vimp"] = feat_diag["feature"].map(lambda c: float(po_v.get(c, 0.0)))

    # --- FSDS (StandardScaler first) ---
    print("FSDS ranking on localized subset (StandardScaler → var → SelectKBest) ...", flush=True)
    if len(g1_tr) < 40:
        g1_tr = g1[g1["item_id"].isin(loc_set)].copy()
        g1_va = g1_tr.sample(frac=0.25, random_state=args.seed) if len(g1_tr) else g1_tr
    res_w1 = run_fsds(g1_tr, g1_va, cols, select_k=args.select_k, seed=args.seed)
    res_w2 = run_fsds(g1_tr, g2_loc, cols, select_k=args.select_k, seed=args.seed)

    ranking = res_w1.get("ranking")
    if ranking is None or not isinstance(ranking, pd.DataFrame) or ranking.empty:
        # fallback ranking from shift diagnostics
        ranking = feat_diag.copy()
        for col in ("cmean_abs", "mmd_loco", "po_vimp"):
            ranking[f"r_{col}"] = ranking[col].rank(ascending=False, method="average")
        ranking["rank_score"] = ranking[["r_cmean_abs", "r_mmd_loco", "r_po_vimp"]].mean(axis=1)
        ranking = ranking.sort_values("rank_score").reset_index(drop=True)
        ranking["rank"] = np.arange(1, len(ranking) + 1)
        ranking["f_score"] = ranking["cmean_abs"]
        ranking["selected"] = (ranking.index < args.select_k).astype(int)

    # persist
    scor.to_csv(args.out_dir / "item_subset_scores.csv", index=False)
    subset.to_csv(args.out_dir / "localized_subset_items.csv", index=False)
    ranking.to_csv(args.out_dir / "fsds_feature_ranking.csv", index=False)
    feat_diag.to_csv(args.out_dir / "feature_shift_diagnostics.csv", index=False)

    # optional GT evaluator (orders / items) — plug in when available
    gt_eval = evaluate_gt(
        loc_items,
        gt_items_path=args.gt_items,
        gt_orders_path=args.gt_orders,
        ks=(50, 100, min(200, len(loc_items)) or 1),
    )
    if gt_eval.get("available"):
        (args.out_dir / "gt_eval.json").write_text(json.dumps(gt_eval, indent=2))
        print(
            "GT eval:",
            {k: gt_eval.get("items", {}).get(k) for k in ("n_gt", "precision@100", "recall@100")},
            flush=True,
        )
    else:
        print("GT eval skipped (pass --gt-items / --gt-orders to enable)", flush=True)

    def _strip(res: Dict) -> Dict:
        out = {k: v for k, v in res.items() if k != "ranking"}
        if "ranking" in res and isinstance(res["ranking"], pd.DataFrame):
            out["top_features"] = res["ranking"].head(args.select_k)["feature"].tolist()
        return out

    blob = {
        "protocol": [
            "Standardization: StandardScaler fit on W1 (top of pipeline)",
            "FE(time) independently on W1 and W2 (gap >= 30d)",
            "subset localization = rank-average(cmean, MMD, PO-risk) over items (standardized X)",
            "visualize localized subset (行为/购买欲变动归因)",
            "FSDS: StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg → feature ranking",
            "optional GT evaluator on orders/items → hit/precision/recall@k",
        ],
        "timeline": {
            "span_days": (t_max - t_min) / DAY,
            "W1": w1.to_dict(),
            "W2": w2.to_dict(),
            "gap_days": gap_d,
        },
        "W1_meta": p1["meta"],
        "W2_meta": p2["meta"],
        "n_features": len(cols),
        "standardization": True,
        "localize_k": len(loc_items),
        "localized_items_head": loc_items[:30],
        "po_risk_global": loc_meta["po_risk_global"],
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
        scor,
        ranking,
        out_path=args.out_dir / "w1w2_mmd_po_localize_fsds.png",
        title=f"Behavior / purchase-intent shift  W1↔W2  gap={gap_d:.0f}d  k={len(loc_items)}",
    )

    # report
    md = [
        "# 行为 / 购买欲变动归因：MMD + PO-risk + conditional-mean → FSDS",
        "",
        "直观 concise 链路（**先 Standardization**）：",
        "1. W1/W2 独立 FE → standardized X",
        "2. **Subset localization**：conditional-mean / MMD² / PO-risk",
        "3. **可视化** 漂移商品 subset（用户行为 & 购买欲变动）",
        "4. **FSDS** 标准流 → 归因特征 ranking",
        "5. （可选）GT 订单/商品 list → hit / precision / recall@k",
        "",
        "## Protocol",
        "0. **StandardScaler** fit on W1（pipeline 最上面）",
        f"1. FE ×2，gap = **{gap_d:.0f}**d",
        f"2. Item subset rank-average(cmean, MMD, PO) → top-**{len(loc_items)}**",
        "3. Viz subset",
        f"4. FSDS: **StandardScaler** → var → SelectKBest(k={args.select_k}) → HGB/LogReg",
        "5. GT evaluator（`--gt-items` / `--gt-orders`）",
        "",
        "## Localized subset (head)",
        "| rank | item_id | cmean | MMD² | PO τ² | n_W1 | n_W2 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in subset.head(12).iterrows():
        md.append(
            f"| {int(r['rank'])} | {int(r['item_id'])} | {r['cmean_l2']:.3f} | "
            f"{r['mmd2']:.4f} | {r['po_tau2_mean']:.4f} | {int(r['n_W1'])} | {int(r['n_W2'])} |"
        )

    def _auc_line(tag: str, res: Dict) -> str:
        if not res.get("ok"):
            return f"- {tag}: **fail** ({res.get('reason')})"
        mods = res.get("models") or {}
        bits = []
        for name in ("hgb", "logreg"):
            if name in mods and "auc" in mods[name]:
                bits.append(f"{name} AUC={mods[name]['auc']:.3f}")
        sel = ", ".join(f"`{c}`" for c in (res.get("selected") or [])[:8])
        return f"- {tag}: n={res.get('n_train')}/{res.get('n_test')} | " + (
            " | ".join(bits) if bits else "ranking only"
        ) + (f" | selected: {sel}" if sel else "")

    md += [
        "",
        f"global PO-risk = **{loc_meta['po_risk_global']:.6f}**",
        "",
        "## FSDS feature ranking",
        "| rank | feature | F-score | selected |",
        "|---:|---|---:|---:|",
    ]
    for _, r in ranking.head(args.select_k).iterrows():
        md.append(
            f"| {int(r['rank'])} | `{r['feature']}` | {float(r.get('f_score', 0)):.3f} | "
            f"{int(r.get('selected', 0))} |"
        )
    md += [
        "",
        "## Holdout",
        _auc_line("W1 user-holdout", res_w1),
        _auc_line("W2 temporal", res_w2),
        "",
        "## GT evaluator",
    ]
    if gt_eval.get("available"):
        it = gt_eval.get("items") or {}
        md.append(
            f"- items: n_gt={int(it.get('n_gt', 0))} | "
            f"P@100={it.get('precision@100', float('nan')):.3f} | "
            f"R@100={it.get('recall@100', float('nan')):.3f}"
        )
        if "orders" in gt_eval:
            od = gt_eval["orders"]
            md.append(
                f"- orders: n={int(od.get('n_orders', 0))} | "
                f"coverage={od.get('order_coverage', float('nan')):.3f}"
            )
    else:
        md.append("- skipped（提供 `--gt-items` / `--gt-orders` 即可评 hit/precision/recall）")
    md.append("")
    (args.out_dir / "W1W2_MMD_PO_LOCALIZE_FSDS_REPORT.md").write_text("\n".join(md))
    print("wrote", args.out_dir)
    print("top features:", ", ".join(ranking.head(8)["feature"].astype(str).tolist()))


if __name__ == "__main__":
    main()
