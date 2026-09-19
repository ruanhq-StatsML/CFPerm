#!/usr/bin/env python3
"""Multi-step FSDS iteration on graph features (locked scope).

Scope: 图谱特征 only. No community / ego / graph-local backends.
Baseline FSDS: StandardScaler → VarianceThreshold → SelectKBest → model.
This script iterates *selection steps* and logs discoveries.

  PYTHONPATH=. python3 scripts/tencent_gr/run_fsds_multistep_iterate.py \\
    --iter-tag iter01

Hourly overnight: bump --iter-tag / try new variants; append to FINDINGS.md.
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif, mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit, StratifiedKFold
from sklearn.preprocessing import StandardScaler

_HERE = Path(__file__).resolve().parent
import sys

if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from time_window_feats import fsds_feature_columns  # noqa: E402
from run_standardize_mmd_fsds import run_fsds, w1w2_candidate_columns  # noqa: E402
from po_risk_fsds import (  # noqa: E402
    blend_cmean_po_scores,
    fit_po_on_windows,
)

ROOT = Path(__file__).resolve().parents[2]


def _matrix(df: pd.DataFrame, cols: Sequence[str]) -> np.ndarray:
    return df.loc[:, list(cols)].to_numpy(dtype=float, copy=True)


def _corr_prune(X: np.ndarray, cols: List[str], thr: float) -> List[str]:
    """Greedy keep: scan by column variance descending, drop |corr|>=thr."""
    if X.shape[1] <= 1:
        return list(cols)
    var = X.var(axis=0)
    order = np.argsort(-var)
    keep_idx: List[int] = []
    for j in order:
        ok = True
        for i in keep_idx:
            c = np.corrcoef(X[:, i], X[:, j])[0, 1]
            if np.isfinite(c) and abs(c) >= thr:
                ok = False
                break
        if ok:
            keep_idx.append(int(j))
    return [cols[i] for i in sorted(keep_idx)]


def _cmean_abs_delta(X1: np.ndarray, X2: np.ndarray) -> np.ndarray:
    return np.abs(X2.mean(axis=0) - X1.mean(axis=0))


def _stability_select(
    X: np.ndarray,
    y: np.ndarray,
    cols: List[str],
    *,
    k: int,
    n_splits: int,
    seed: int,
    score_fn,
) -> Tuple[List[str], pd.DataFrame]:
    """Multi-fold SelectKBest; π_j = fraction of folds feature in TopK."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    hits = np.zeros(len(cols), dtype=float)
    score_sum = np.zeros(len(cols), dtype=float)
    n_ok = 0
    for tr, _ in skf.split(X, y):
        if len(np.unique(y[tr])) < 2:
            continue
        kk = min(k, X.shape[1], max(1, len(tr) - 1))
        sel = SelectKBest(score_fn, k=kk)
        sel.fit(X[tr], y[tr])
        hits += sel.get_support().astype(float)
        sc = np.nan_to_num(sel.scores_, nan=0.0)
        score_sum += sc
        n_ok += 1
    if n_ok == 0:
        return cols[:k], pd.DataFrame({"feature": cols, "pi": 0.0, "mean_score": 0.0})
    pi = hits / n_ok
    mean_sc = score_sum / n_ok
    rank = np.argsort(-pi, kind="mergesort")
    # tie-break by mean score
    order = sorted(range(len(cols)), key=lambda j: (-pi[j], -mean_sc[j]))
    selected = [cols[j] for j in order[:k]]
    tab = (
        pd.DataFrame({"feature": cols, "pi": pi, "mean_score": mean_sc})
        .sort_values(["pi", "mean_score"], ascending=False)
        .reset_index(drop=True)
    )
    tab["rank"] = np.arange(1, len(tab) + 1)
    tab["selected"] = tab["feature"].isin(selected).astype(int)
    return selected, tab


@dataclass
class VariantResult:
    name: str
    selected: List[str]
    ranking: pd.DataFrame
    metrics: Dict
    sec: float
    notes: str = ""


def eval_models(
    Xtr: np.ndarray,
    ytr: np.ndarray,
    Xte: np.ndarray,
    yte: np.ndarray,
    *,
    seed: int,
) -> Dict:
    out: Dict = {}
    if len(yte) == 0 or len(np.unique(yte)) < 2 or len(np.unique(ytr)) < 2:
        out["note"] = "need both classes"
        return out
    hgb = HistGradientBoostingClassifier(
        max_depth=6, learning_rate=0.08, max_iter=120, random_state=seed, class_weight="balanced"
    )
    hgb.fit(Xtr, ytr)
    ph = hgb.predict_proba(Xte)[:, 1]
    out["hgb"] = {
        "auc": float(roc_auc_score(yte, ph)),
        "ap": float(average_precision_score(yte, ph)),
    }
    lr = LogisticRegression(max_iter=400, C=0.5, class_weight="balanced", random_state=seed)
    lr.fit(Xtr, ytr)
    pl = lr.predict_proba(Xte)[:, 1]
    out["logreg"] = {
        "auc": float(roc_auc_score(yte, pl)),
        "ap": float(average_precision_score(yte, pl)),
    }
    return out


def run_variant(
    name: str,
    g_tr: pd.DataFrame,
    g_te: pd.DataFrame,
    g_w2_ref: pd.DataFrame,
    cols: List[str],
    *,
    k: int,
    seed: int,
    builder: Callable,
    use_official_fsds: bool = True,
) -> VariantResult:
    """Select via builder; score via official FSDS pipeline when possible.

    Official FSDS (locked stats method): Scaler → Var → SelectKBest → HGB/LogReg,
    fit on W1-train only; holdout / W2 never enter selection.
    """
    t0 = time.time()
    # If builder is pure FSDS on full cols, call official run_fsds end-to-end.
    if use_official_fsds and name.split("|")[0] in ("A_baseline_F", "B_baseline_MI"):
        # still use builder for ranking table consistency; metrics from run_fsds
        selected, ranking, notes = builder(g_tr, g_w2_ref, cols, k=k, seed=seed)
        res = run_fsds(g_tr, g_te, cols, select_k=k, seed=seed)
        metrics = {"hgb": res.get("models", {}).get("hgb", {}), "logreg": res.get("models", {}).get("logreg", {})}
        if res.get("ok") and isinstance(res.get("ranking"), pd.DataFrame):
            ranking = res["ranking"].rename(columns={"f_score": "score"})
            selected = list(res.get("selected") or selected)
        notes = notes + " | metrics=official_run_fsds"
    else:
        selected, ranking, notes = builder(g_tr, g_w2_ref, cols, k=k, seed=seed)
        if use_official_fsds:
            # Fuse: builder proposes final column set (or a shortlist).
            # Official FSDS runs on exactly that set — W2 labels never select.
            pool = list(selected) if selected else list(cols)
            res = run_fsds(g_tr, g_te, pool, select_k=min(k, len(pool)), seed=seed)
            if res.get("ok"):
                selected = list(res.get("selected") or selected)
                if isinstance(res.get("ranking"), pd.DataFrame):
                    ranking = res["ranking"].rename(columns={"f_score": "score"})
                metrics = {
                    "hgb": res.get("models", {}).get("hgb", {}),
                    "logreg": res.get("models", {}).get("logreg", {}),
                }
                notes = notes + " | fused→official_FSDS"
            else:
                sc = StandardScaler()
                Xtr = sc.fit_transform(_matrix(g_tr, selected))
                Xte = sc.transform(_matrix(g_te, selected))
                metrics = eval_models(
                    Xtr,
                    g_tr["y_convert"].to_numpy(int),
                    Xte,
                    g_te["y_convert"].to_numpy(int),
                    seed=seed,
                )
                notes = notes + " | fallback_eval_models"
        else:
            sc = StandardScaler()
            Xtr = sc.fit_transform(_matrix(g_tr, selected))
            Xte = sc.transform(_matrix(g_te, selected))
            metrics = eval_models(
                Xtr,
                g_tr["y_convert"].to_numpy(int),
                Xte,
                g_te["y_convert"].to_numpy(int),
                seed=seed,
            )
    return VariantResult(
        name=name,
        selected=selected,
        ranking=ranking if isinstance(ranking, pd.DataFrame) else pd.DataFrame(),
        metrics=metrics,
        sec=float(time.time() - t0),
        notes=notes,
    )


# ---- builders: (g_tr, g_w2_for_delta, cols, k, seed) -> selected, ranking, notes ----

def build_baseline_f(g_tr, g_w2, cols, *, k, seed):
    X = _matrix(g_tr, cols)
    y = g_tr["y_convert"].to_numpy(int)
    pipe_k = min(k, X.shape[1], max(1, len(y) - 1))
    vt = VarianceThreshold(1e-8)
    Xs = StandardScaler().fit_transform(X)
    Xv = vt.fit_transform(Xs)
    cols_v = [c for c, m in zip(cols, vt.get_support()) if m]
    sel = SelectKBest(f_classif, k=min(pipe_k, len(cols_v)))
    sel.fit(Xv, y)
    selected = [c for c, m in zip(cols_v, sel.get_support()) if m]
    ranking = (
        pd.DataFrame({"feature": cols_v, "score": sel.scores_})
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    ranking["selected"] = ranking["feature"].isin(selected).astype(int)
    return selected, ranking, "Scaler→Var→SelectKBest(F)"


def build_baseline_mi(g_tr, g_w2, cols, *, k, seed):
    X = StandardScaler().fit_transform(_matrix(g_tr, cols))
    y = g_tr["y_convert"].to_numpy(int)
    vt = VarianceThreshold(1e-8)
    Xv = vt.fit_transform(X)
    cols_v = [c for c, m in zip(cols, vt.get_support()) if m]
    sel = SelectKBest(
        lambda a, b: mutual_info_classif(a, b, random_state=seed),
        k=min(k, len(cols_v)),
    )
    sel.fit(Xv, y)
    selected = [c for c, m in zip(cols_v, sel.get_support()) if m]
    ranking = (
        pd.DataFrame({"feature": cols_v, "score": sel.scores_})
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    ranking["selected"] = ranking["feature"].isin(selected).astype(int)
    return selected, ranking, "Scaler→Var→SelectKBest(MI)"


def build_cmean_then_f(g_tr, g_w2, cols, *, k, seed):
    """Step1: keep top by |δ| (W1-train vs W2 ref means); Step2: F on remaining."""
    X1 = _matrix(g_tr, cols)
    X2 = _matrix(g_w2, cols)
    dlt = _cmean_abs_delta(X1, X2)
    pre_n = min(len(cols), max(k * 2, k + 5, int(0.6 * len(cols))))
    pre_idx = np.argsort(-dlt)[:pre_n]
    pre_cols = [cols[i] for i in pre_idx]
    return build_baseline_f(g_tr, g_w2, pre_cols, k=k, seed=seed)[:2] + (
        f"cmean|δ| prefilter→{pre_n} then F→{k}",
    )


def build_f_then_corr_prune(g_tr, g_w2, cols, *, k, seed):
    """Step1: SelectKBest wider; Step2: corr prune to ~k."""
    wide = min(len(cols), max(k * 2, k + 8))
    sel_cols, ranking, _ = build_baseline_f(g_tr, g_w2, cols, k=wide, seed=seed)
    X = StandardScaler().fit_transform(_matrix(g_tr, sel_cols))
    pruned = _corr_prune(X, sel_cols, thr=0.92)
    if len(pruned) > k:
        # keep highest F among pruned
        score_map = dict(zip(ranking["feature"], ranking["score"]))
        pruned = sorted(pruned, key=lambda c: -float(score_map.get(c, 0.0)))[:k]
    ranking = ranking.copy()
    ranking["selected"] = ranking["feature"].isin(pruned).astype(int)
    return pruned, ranking, f"F wide={wide} → corr-prune@0.92 →{len(pruned)}"


def build_stable_pi_f(g_tr, g_w2, cols, *, k, seed):
    X = StandardScaler().fit_transform(_matrix(g_tr, cols))
    y = g_tr["y_convert"].to_numpy(int)
    vt = VarianceThreshold(1e-8)
    Xv = vt.fit_transform(X)
    cols_v = [c for c, m in zip(cols, vt.get_support()) if m]
    selected, tab = _stability_select(
        Xv, y, cols_v, k=k, n_splits=5, seed=seed, score_fn=f_classif
    )
    return selected, tab.rename(columns={"mean_score": "score"}), "5-fold π-stable SelectKBest(F)"


def build_cmean_stable(g_tr, g_w2, cols, *, k, seed):
    """cmean prefilter then π-stable F."""
    X1 = _matrix(g_tr, cols)
    X2 = _matrix(g_w2, cols)
    dlt = _cmean_abs_delta(X1, X2)
    pre_n = min(len(cols), max(k * 2, k + 5, int(0.6 * len(cols))))
    pre_cols = [cols[i] for i in np.argsort(-dlt)[:pre_n]]
    selected, tab, _ = build_stable_pi_f(g_tr, g_w2, pre_cols, k=k, seed=seed)
    return selected, tab, f"cmean pre→{pre_n} + π-stable F→{k}"


def build_hgb_importance_refine(g_tr, g_w2, cols, *, k, seed):
    """Step1: F wide screen; Step2: HGB impurity importance on screen; take top-k."""
    wide = min(len(cols), max(k * 2, k + 8))
    screen, ranking_f, _ = build_baseline_f(g_tr, g_w2, cols, k=wide, seed=seed)
    X = StandardScaler().fit_transform(_matrix(g_tr, screen))
    y = g_tr["y_convert"].to_numpy(int)
    hgb = HistGradientBoostingClassifier(
        max_depth=6, learning_rate=0.08, max_iter=80, random_state=seed, class_weight="balanced"
    )
    hgb.fit(X, y)
    # sklearn HGB has no feature_importances_ until recent — use permutation on train proxy:
    from sklearn.inspection import permutation_importance

    imp = permutation_importance(hgb, X, y, n_repeats=5, random_state=seed, scoring="roc_auc")
    scores = imp.importances_mean
    order = np.argsort(-scores)[:k]
    selected = [screen[i] for i in order]
    tab = (
        pd.DataFrame({"feature": screen, "score": scores})
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
    tab["rank"] = np.arange(1, len(tab) + 1)
    tab["selected"] = tab["feature"].isin(selected).astype(int)
    return selected, tab, f"F screen={wide} → HGB perm-imp →{k}"


def build_cmean_stable_mi(g_tr, g_w2, cols, *, k, seed):
    """Like F_cmean_stable but SelectKBest uses mutual_info (auto_feats path)."""
    X1 = _matrix(g_tr, cols)
    X2 = _matrix(g_w2, cols)
    dlt = _cmean_abs_delta(X1, X2)
    pre_n = min(len(cols), max(k * 2, k + 5, int(0.6 * len(cols))))
    pre_cols = [cols[i] for i in np.argsort(-dlt)[:pre_n]]
    X = StandardScaler().fit_transform(_matrix(g_tr, pre_cols))
    y = g_tr["y_convert"].to_numpy(int)
    vt = VarianceThreshold(1e-8)
    Xv = vt.fit_transform(X)
    cols_v = [c for c, m in zip(pre_cols, vt.get_support()) if m]
    selected, tab = _stability_select(
        Xv,
        y,
        cols_v,
        k=k,
        n_splits=5,
        seed=seed,
        score_fn=lambda a, b: mutual_info_classif(a, b, random_state=seed),
    )
    return (
        selected,
        tab.rename(columns={"mean_score": "score"}),
        f"cmean pre→{pre_n} + π-stable MI→{k}",
    )


def build_combined(g_tr, g_w2, cols, *, k, seed):
    """ empirically combined recipe (iter01–03):

        soft cmean guidance → π-stable F (same as F_cmean_stable) → official FSDS

    Explicitly **omits** hard/soft corr prune and J*-only replace (hurt W2).
    """
    selected, ranking, note_f = build_cmean_stable(g_tr, g_w2, cols, k=k, seed=seed)
    note = f"COMBINED(=cmean+π→FSDS): {note_f}"
    return selected, ranking, note


def build_combined_mi(g_tr, g_w2, cols, *, k, seed):
    selected, ranking, note_f = build_cmean_stable_mi(g_tr, g_w2, cols, k=k, seed=seed)
    return selected, ranking, f"COMBINED-MI(=cmean+π-MI→FSDS): {note_f}"


def build_po_vimp_fsds(g_tr, g_w2, cols, *, k, seed):
    """PO-VIMP ranks features (period W); then official FSDS on top pool."""
    po = fit_po_on_windows(g_tr, g_w2, cols, seed=seed, max_n=6000)
    tab = po["feature_table"].rename(columns={"po_vimp": "score"})
    # tight pool so PO actually shapes the set (not a no-op on 21 feats)
    pre_n = min(len(cols), max(k, k + 3))
    pre_cols = list(tab["feature"].head(pre_n))
    selected, ranking, _ = build_baseline_f(g_tr, g_w2, pre_cols, k=k, seed=seed)
    return (
        selected,
        ranking,
        f"PO-VIMP pre→{pre_n} (risk={po['risk']:.6g}) then FSDS-F→{k}",
    )


def build_combined_po(g_tr, g_w2, cols, *, k, seed):
    """Data-science combo: cmean |δ| ⋈ PO-VIMP → π-stable F → FSDS.

    PO fit once (W=period). Helps rank shift-relevant graph feats without
    claiming ATE. α=0.5 blend; tight guided pool so PO is not a no-op.
    """
    X1 = _matrix(g_tr, cols)
    X2 = _matrix(g_w2, cols)
    dlt = _cmean_abs_delta(X1, X2)
    po = fit_po_on_windows(g_tr, g_w2, cols, seed=seed, max_n=6000)
    blend = blend_cmean_po_scores(cols, dlt, po["vimp"], alpha=0.5)
    pre_n = min(len(cols), max(k, k + 3))
    guided = list(blend["feature"].head(pre_n))
    selected, ranking, note_pi = build_stable_pi_f(g_tr, g_w2, guided, k=k, seed=seed)
    note = (
        f"COMBINED-PO: cmean⋈PO-VIMP→{pre_n} | {note_pi} | "
        f"PO-risk={po['risk']:.6g} | →FSDS"
    )
    ranking = ranking.copy()
    return selected, ranking, note


def build_soft_corr_prune(g_tr, g_w2, cols, *, k, seed):
    """FSDS-wide then soft corr prune @0.98 (less aggressive than 0.92)."""
    wide = min(len(cols), max(k * 2, k + 8))
    sel_cols, ranking, _ = build_baseline_f(g_tr, g_w2, cols, k=wide, seed=seed)
    X = StandardScaler().fit_transform(_matrix(g_tr, sel_cols))
    pruned = _corr_prune(X, sel_cols, thr=0.98)
    if len(pruned) > k:
        score_map = dict(zip(ranking["feature"], ranking["score"]))
        pruned = sorted(pruned, key=lambda c: -float(score_map.get(c, 0.0)))[:k]
    ranking = ranking.copy()
    ranking["selected"] = ranking["feature"].isin(pruned).astype(int)
    return pruned, ranking, f"F wide={wide} → soft-corr@0.98 →{len(pruned)}"


def build_stable_pi_f3(g_tr, g_w2, cols, *, k, seed):
    """π-stable with 3 folds (rare-positive friendly vs 5-fold)."""
    X = StandardScaler().fit_transform(_matrix(g_tr, cols))
    y = g_tr["y_convert"].to_numpy(int)
    vt = VarianceThreshold(1e-8)
    Xv = vt.fit_transform(X)
    cols_v = [c for c, m in zip(cols, vt.get_support()) if m]
    n_pos = int(y.sum())
    n_splits = max(2, min(3, n_pos))
    selected, tab = _stability_select(
        Xv, y, cols_v, k=k, n_splits=n_splits, seed=seed, score_fn=f_classif
    )
    return (
        selected,
        tab.rename(columns={"mean_score": "score"}),
        f"{n_splits}-fold π-stable SelectKBest(F) (n_pos={n_pos})",
    )


def build_delta_share_then_fsds_cols(g_tr, g_w2, cols, *, k, seed):
    """Method fusion: J*=TopShare(|δ|) then keep those cols for official FSDS."""
    X1 = _matrix(g_tr, cols)
    X2 = _matrix(g_w2, cols)
    dlt = _cmean_abs_delta(X1, X2)
    energy = dlt ** 2
    share = energy / (energy.sum() + 1e-12)
    # keep until cumulative share >= 0.8 or at least k feats
    order = np.argsort(-share)
    cum = np.cumsum(share[order])
    n_keep = int(np.searchsorted(cum, 0.80) + 1)
    n_keep = max(k, min(len(cols), n_keep))
    pre_cols = [cols[i] for i in order[:n_keep]]
    # second step still FSDS select to k
    selected, ranking, _ = build_baseline_f(g_tr, g_w2, pre_cols, k=k, seed=seed)
    return selected, ranking, f"J* cumshare≥0.8 →{n_keep} cols then FSDS-F→{k}"


VARIANTS: Dict[str, Callable] = {
    "A_baseline_F": build_baseline_f,
    "B_baseline_MI": build_baseline_mi,
    "C_cmean_then_F": build_cmean_then_f,
    "D_F_corr_prune": build_f_then_corr_prune,
    "E_stable_pi_F": build_stable_pi_f,
    "F_cmean_stable": build_cmean_stable,
    "G_F_then_HGB_perm": build_hgb_importance_refine,
    "H_soft_corr": build_soft_corr_prune,
    "I_stable_pi_F3": build_stable_pi_f3,
    "J_delta_share_FSDS": build_delta_share_then_fsds_cols,
    "Z_combined": build_combined,
    "Z_combined_MI": build_combined_mi,
    "P_po_vimp_FSDS": build_po_vimp_fsds,
    "Z_combined_PO": build_combined_po,
}


def _df_md(df: pd.DataFrame) -> str:
    """Markdown table without tabulate dependency."""
    cols = list(df.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if isinstance(v, float):
                cells.append(f"{v:.4f}" if np.isfinite(v) else "")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def jaccard(a: Sequence[str], b: Sequence[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    return float(len(sa & sb) / len(sa | sb))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--w1-grid",
        type=Path,
        default=ROOT / "results" / "tencent_gr_localize_fsds_time" / "W1_localized_grid.parquet",
    )
    ap.add_argument(
        "--w2-grid",
        type=Path,
        default=ROOT
        / "results"
        / "tencent_gr_localize_fsds_time"
        / "W2_localized_grid_sample.parquet",
    )
    ap.add_argument("--select-k", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--iter-tag", type=str, default="iter01")
    ap.add_argument(
        "--variants",
        type=str,
        default=(
            "A_baseline_F,Z_combined,P_po_vimp_FSDS,Z_combined_PO"
        ),
    )
    ap.add_argument(
        "--fuse-official-fsds",
        action="store_true",
        default=True,
        help="Score via official run_fsds (Scaler→Var→SelectKBest→model)",
    )
    ap.add_argument("--no-fuse-official-fsds", action="store_false", dest="fuse_official_fsds")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "tencent_gr_fsds_iterate",
    )
    args = ap.parse_args()
    out = args.out_dir / args.iter_tag
    out.mkdir(parents=True, exist_ok=True)

    g1 = pd.read_parquet(args.w1_grid)
    g2 = pd.read_parquet(args.w2_grid)
    raw = [c for c in fsds_feature_columns(g1) if c in g2.columns]
    cols = w1w2_candidate_columns(raw)
    print(f"grids W1={len(g1)} W2={len(g2)} feats={len(cols)}", flush=True)

    gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=args.seed)
    tr_idx, va_idx = next(gss.split(g1, groups=g1["user_id"].to_numpy()))
    g_tr, g_va = g1.iloc[tr_idx].copy(), g1.iloc[va_idx].copy()
    print(f"split train={len(g_tr)} va={len(g_va)} pos_tr={g_tr['y_convert'].mean():.4f}", flush=True)

    names = [x.strip() for x in args.variants.split(",") if x.strip()]
    results: List[VariantResult] = []
    for name in names:
        if name not in VARIANTS:
            print(f"skip unknown variant {name}", flush=True)
            continue
        print(f"run {name} ...", flush=True)
        # W1-holdout eval
        r_va = run_variant(
            name + "|W1hold",
            g_tr,
            g_va,
            g2,
            cols,
            k=args.select_k,
            seed=args.seed,
            builder=VARIANTS[name],
            use_official_fsds=args.fuse_official_fsds,
        )
        # W2 temporal: official FSDS transform path — select on train only
        selected = r_va.selected
        if args.fuse_official_fsds:
            res_w2 = run_fsds(g_tr, g2, selected, select_k=min(args.select_k, len(selected)), seed=args.seed)
            m_w2 = {
                "hgb": res_w2.get("models", {}).get("hgb", {}),
                "logreg": res_w2.get("models", {}).get("logreg", {}),
            }
            if not res_w2.get("ok"):
                sc = StandardScaler()
                Xtr = sc.fit_transform(_matrix(g_tr, selected))
                Xw2 = sc.transform(_matrix(g2, selected))
                m_w2 = eval_models(
                    Xtr,
                    g_tr["y_convert"].to_numpy(int),
                    Xw2,
                    g2["y_convert"].to_numpy(int),
                    seed=args.seed,
                )
        else:
            sc = StandardScaler()
            Xtr = sc.fit_transform(_matrix(g_tr, selected))
            Xw2 = sc.transform(_matrix(g2, selected))
            m_w2 = eval_models(
                Xtr,
                g_tr["y_convert"].to_numpy(int),
                Xw2,
                g2["y_convert"].to_numpy(int),
                seed=args.seed,
            )
        r_va.metrics["W2_temporal"] = m_w2
        results.append(r_va)
        r_va.ranking.to_csv(out / f"ranking_{name}.csv", index=False)
        print(
            f"  selected={selected[:5]}... "
            f"W1 hgb={r_va.metrics.get('hgb', {}).get('auc')} "
            f"W2 hgb={m_w2.get('hgb', {}).get('auc')} "
            f"W2 ap={m_w2.get('hgb', {}).get('ap')}  [{r_va.notes}]",
            flush=True,
        )

    # compare to baseline set
    base = next((r for r in results if r.name.startswith("A_baseline_F")), None)
    rows = []
    for r in results:
        h1 = r.metrics.get("hgb", {}) or {}
        l1 = r.metrics.get("logreg", {}) or {}
        w2 = (r.metrics.get("W2_temporal") or {}).get("hgb", {}) or {}
        w2l = (r.metrics.get("W2_temporal") or {}).get("logreg", {}) or {}
        jac = jaccard(base.selected, r.selected) if base else float("nan")
        rows.append(
            {
                "variant": r.name,
                "notes": r.notes,
                "n_selected": len(r.selected),
                "jaccard_vs_baseline": jac,
                "W1_hgb_auc": h1.get("auc"),
                "W1_hgb_ap": h1.get("ap"),
                "W1_logreg_auc": l1.get("auc"),
                "W2_hgb_auc": w2.get("auc"),
                "W2_hgb_ap": w2.get("ap"),
                "W2_logreg_auc": w2l.get("auc"),
                "sec": r.sec,
                "top5": ",".join(r.selected[:5]),
            }
        )
    cmp = pd.DataFrame(rows)
    cmp.to_csv(out / "variant_compare.csv", index=False)

    # discoveries
    lines = [
        f"# FSDS multi-step iteration `{args.iter_tag}`",
        "",
        f"- grids: `{args.w1_grid.name}` / `{args.w2_grid.name}`",
        f"- n_train={len(g_tr)} n_va={len(g_va)} n_w2={len(g2)} feats={len(cols)} k={args.select_k}",
        f"- fuse_official_fsds={args.fuse_official_fsds} (Scaler→Var→SelectKBest→model; W2 never selects)",
        f"- pos_train={float(g_tr['y_convert'].mean()):.6f} (rare-positive regime)",
        "",
        "## Variant table",
        "",
        _df_md(cmp),
        "",
        "## Findings",
        "",
    ]
    if base is not None and len(cmp):
        best_w2 = cmp.sort_values("W2_hgb_auc", ascending=False).iloc[0]
        best_w1 = cmp.sort_values("W1_hgb_auc", ascending=False).iloc[0]
        lines.append(
            f"- Best **W2** HGB AUC: `{best_w2['variant']}` = {best_w2['W2_hgb_auc']:.4f} "
            f"(baseline A={cmp.loc[cmp['variant'].str.startswith('A_baseline_F'), 'W2_hgb_auc'].values[0]:.4f})"
        )
        lines.append(
            f"- Best **W1-hold** HGB AUC: `{best_w1['variant']}` = {best_w1['W1_hgb_auc']:.4f}"
        )
        # redundancy note from corr prune
        drow = cmp[cmp["variant"].str.contains("corr_prune", na=False)]
        if len(drow) and base is not None:
            lines.append(
                f"- Corr-prune Jaccard vs baseline={float(drow.iloc[0]['jaccard_vs_baseline']):.3f} "
                f"(drops near-duplicate share/credit twins when present)"
            )
        # cmean alignment
        crow = cmp[cmp["variant"].str.contains("cmean", na=False)]
        if len(crow):
            lines.append(
                "- cmean-|δ| prefilter aligns FS with period-shift guidance (same δ story as drill doc)"
            )
        lines.append(
            "- Scope locked: graph features only; selection steps only; no graph algorithms."
        )
    (out / "FINDINGS.md").write_text("\n".join(lines) + "\n")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "iter_tag": args.iter_tag,
                "n_train": len(g_tr),
                "n_va": len(g_va),
                "n_w2": len(g2),
                "n_feats": len(cols),
                "compare": rows,
            },
            indent=2,
        )
    )

    # append master log
    master = args.out_dir / "ITERATION_LOG.md"
    prev = master.read_text() if master.exists() else "# FSDS multi-step overnight log\n\n"
    master.write_text(
        prev
        + f"\n## {args.iter_tag} ({time.strftime('%Y-%m-%d %H:%M UTC')})\n\n"
        + _df_md(cmp)
        + "\n\n"
        + f"See `{args.iter_tag}/FINDINGS.md`.\n"
    )
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
