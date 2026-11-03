"""Text multiclass: feature-level vs sample-level weights on FULL eval.

No EWMA. No first-n subset KPI as primary — report overall accuracy and
mean OVA-AUC on the full held-out LIVE set.

Schemes
-------
  uniform    : CE, flat sample weights, raw features
  feature_w  : CE, flat sample weights, X *= (1 + α · vimp)   ← feature-level
  sample_w   : CE, sample_weight from membership, raw features
  both       : feature_w + sample_w

Curves: overall acc / mean OVA-AUC vs novel-sample budget.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from text_feature_sample_weights import load_text_boards_full
from weighted_cold_start import (
    SEED,
    attribute_rf_domain,
    cold_start_partial_fit,
    localize_subset,
    make_sgd,
    weights_from_attribution,
)

LOSS_NAME = "CE (log_loss)"
SCHEMES = ("uniform", "feature_w", "sample_w", "both")


def feature_level_scale(
    X: np.ndarray,
    vimp: np.ndarray,
    *,
    alpha: float = 3.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Feature-level weights: scale each column by (1 + α · vimp_j).

    vimp is L1-normalized RF-domain importance (LIVE vs REF).
    Returns (X_scaled, scale_vector).
    """
    v = np.asarray(vimp, float).ravel()
    v = v / (v.sum() + 1e-12)
    # stretch so mean scale ≈ 1, high-VIMP coords get more mass
    scale = 1.0 + float(alpha) * (v / (v.max() + 1e-12))
    X = np.asarray(X, float) * scale.reshape(1, -1)
    return X, scale


def _ova_auc_mean(
    y_true: np.ndarray,
    proba: np.ndarray,
    classes: np.ndarray,
    novel: Sequence[int],
) -> float:
    y_true = np.asarray(y_true).ravel()
    classes = list(np.asarray(classes).ravel())
    aucs = []
    for c in novel:
        c = int(c)
        if c not in classes:
            continue
        col = classes.index(c)
        binary = (y_true == c).astype(int)
        if binary.min() == binary.max():
            continue
        aucs.append(float(roc_auc_score(binary, proba[:, col])))
    return float(np.mean(aucs)) if aucs else float("nan")


def fit_full_scheme(
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    X_tr: np.ndarray,
    y_tr: np.ndarray,
    *,
    scheme: str,
    seed: int,
    classes: np.ndarray,
    novel: Sequence[int],
    top_k: int = 16,
    lam: float = 2.5,
    feat_alpha: float = 3.0,
) -> Tuple[SGDClassifier, StandardScaler, Optional[np.ndarray]]:
    """Fit CE head under uniform / feature_w / sample_w / both. No EWMA."""
    scheme = str(scheme).lower()
    sc = StandardScaler().fit(X_ref)
    X0 = sc.transform(X_ref)
    X1 = sc.transform(X_tr)
    classes = np.asarray(classes, int)
    clf = make_sgd(seed=seed)
    clf.partial_fit(X0, y_ref, classes=classes)

    scale = None
    w = np.ones(len(y_tr), dtype=float)

    if scheme != "uniform":
        vimp = attribute_rf_domain(X0, X1, seed=seed, n_estimators=80)
        mem, feat = localize_subset(X0, X1, vimp, top_k=top_k)
        # private-live fuse for sample channel
        mem_f = mem.copy()
        if novel:
            mem_f = np.maximum(mem_f, np.isin(y_tr, list(novel)).astype(float))

        if scheme in ("feature_w", "both"):
            X1, scale = feature_level_scale(X1, vimp, alpha=feat_alpha)
            # same geometry on REF + re-init so coefs live in scaled coords
            X0 = X0 * scale.reshape(1, -1)
            clf = make_sgd(seed=seed)
            clf.partial_fit(X0, y_ref, classes=classes)

        if scheme in ("sample_w", "both"):
            wt = weights_from_attribution(
                mem_f, feat, X1.shape[1], lam=lam, mode="fused", temper=True,
            )
            w = wt["sample_weight"]
        elif scheme == "feature_w":
            w = np.ones(len(y_tr), dtype=float)

    w = w * (len(w) / max(w.sum(), 1e-12))
    clf = cold_start_partial_fit(
        clf, X1, y_tr, sample_weight=w, classes=classes, n_passes=4,
    )
    return clf, sc, scale


def predict_full(
    clf: SGDClassifier,
    sc: StandardScaler,
    scale: Optional[np.ndarray],
    X: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    Xs = sc.transform(np.asarray(X, float))
    if scale is not None:
        Xs = Xs * scale.reshape(1, -1)
    return clf.predict(Xs), clf.predict_proba(Xs)


def eval_full_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    proba: np.ndarray,
    classes: np.ndarray,
    *,
    known: Sequence[int],
    novel: Sequence[int],
) -> Dict:
    y_true = np.asarray(y_true).ravel()
    out = {
        "acc": float(accuracy_score(y_true, y_pred)),
        "ova_auc_novel_mean": _ova_auc_mean(y_true, proba, classes, novel),
    }
    if known:
        m = np.isin(y_true, list(known))
        out["known_acc"] = float(np.mean(y_pred[m] == y_true[m])) if m.any() else float("nan")
    if novel:
        m = np.isin(y_true, list(novel))
        out["novel_acc"] = float(np.mean(y_pred[m] == y_true[m])) if m.any() else float("nan")
    return out


def run_board_full_curves(
    board: Dict,
    *,
    seed: int = SEED,
    budgets: Optional[Sequence[int]] = None,
    schemes: Sequence[str] = SCHEMES,
    feat_alpha: float = 3.0,
) -> Dict:
    """Learning curves on FULL eval: vary #novel adapt samples."""
    budgets = list(budgets or (20, 40, 60, 80, 100, 150, 200))
    known = list(board["known_labels"])
    novel = list(board["novel_labels"])
    X_ref = np.asarray(board["X_ref"], float)
    y_ref = np.asarray(board["y_ref"]).ravel()
    X_live = np.asarray(board["X_live"], float)
    y_live = np.asarray(board["y_live"]).ravel()

    rng = np.random.default_rng(int(seed))
    # fixed full holdout on LIVE (stratified-ish)
    idx = rng.permutation(len(y_live))
    n_ev = max(80, len(y_live) // 3)
    ev, pool = idx[:n_ev], idx[n_ev:]
    X_ev, y_ev = X_live[ev], y_live[ev]
    X_pool, y_pool = X_live[pool], y_live[pool]

    known_pool = np.flatnonzero(np.isin(y_pool, known))
    novel_pool = np.flatnonzero(np.isin(y_pool, novel))
    rng.shuffle(known_pool)
    rng.shuffle(novel_pool)

    classes = np.unique(np.concatenate([y_ref, y_live]))
    curves = {sch: [] for sch in schemes}

    for b in budgets:
        b = int(b)
        take_n = novel_pool[: min(b, len(novel_pool))]
        # keep known context ~ 2x novel budget
        take_k = known_pool[: min(max(40, 2 * len(take_n)), len(known_pool))]
        tr = np.concatenate([take_k, take_n]) if len(take_n) else take_k
        if tr.size < 8:
            continue
        X_tr, y_tr = X_pool[tr], y_pool[tr]

        for sch in schemes:
            clf, sc, scale = fit_full_scheme(
                X_ref, y_ref, X_tr, y_tr,
                scheme=sch, seed=seed + b, classes=classes, novel=novel,
                feat_alpha=feat_alpha,
            )
            pred, proba = predict_full(clf, sc, scale, X_ev)
            m = eval_full_metrics(
                y_ev, pred, proba, classes, known=known, novel=novel,
            )
            curves[sch].append({
                "budget": b,
                "n_novel_train": int(len(take_n)),
                "n_known_train": int(len(take_k)),
                **m,
            })

    # endpoint summary = last budget
    endpoint = {}
    for sch in schemes:
        endpoint[sch] = curves[sch][-1] if curves[sch] else {}

    def _delta(metric):
        out = {}
        base = endpoint.get("uniform", {}).get(metric)
        for sch in schemes:
            if sch == "uniform" or base is None or base != base:
                continue
            v = endpoint.get(sch, {}).get(metric)
            out[sch] = float(v - base) if v is not None and v == v else None
        return out

    return {
        "dataset": board["name"],
        "note": board.get("note"),
        "loss": LOSS_NAME,
        "ewma": False,
        "eval": "full_holdout_LIVE",
        "n_eval": int(len(y_ev)),
        "budgets": list(budgets),
        "schemes": list(schemes),
        "curves": curves,
        "endpoint": endpoint,
        "delta_acc_vs_uniform": _delta("acc"),
        "delta_ova_auc_vs_uniform": _delta("ova_auc_novel_mean"),
        "delta_novel_acc_vs_uniform": _delta("novel_acc"),
    }


def run_text_full_weight_curves(
    *,
    seed: int = SEED,
    n_first_board: int = 400,
    budgets: Optional[Sequence[int]] = None,
    feat_alpha: float = 3.0,
) -> Dict:
    """Text-only boards; full-eval curves for feature/sample/both weights."""
    # Full LIVE mass (not first-n subset); cnae9 capped inside loader.
    boards = load_text_boards_full(seed=seed, n_per_novel=n_first_board)
    rows = []
    for b in boards:
        print("[full-curves]", b["name"], "n_live", len(b["y_live"]), flush=True)
        rows.append(
            run_board_full_curves(
                b, seed=seed, budgets=budgets, feat_alpha=feat_alpha,
            )
        )
    # aggregate endpoint deltas
    def _mean_delta(key, sch):
        vals = [
            r[key].get(sch)
            for r in rows
            if r.get(key, {}).get(sch) is not None
        ]
        return float(np.mean(vals)) if vals else None

    return {
        "domain": "text_multiclass_only",
        "tabular_excluded": True,
        "loss": LOSS_NAME,
        "ewma": False,
        "primary_kpi": "full_holdout overall acc + mean novel OVA-AUC",
        "seed": int(seed),
        "boards": rows,
        "mean_delta_acc_feature_w": _mean_delta("delta_acc_vs_uniform", "feature_w"),
        "mean_delta_acc_sample_w": _mean_delta("delta_acc_vs_uniform", "sample_w"),
        "mean_delta_acc_both": _mean_delta("delta_acc_vs_uniform", "both"),
        "mean_delta_ova_feature_w": _mean_delta("delta_ova_auc_vs_uniform", "feature_w"),
        "mean_delta_ova_sample_w": _mean_delta("delta_ova_auc_vs_uniform", "sample_w"),
        "mean_delta_ova_both": _mean_delta("delta_ova_auc_vs_uniform", "both"),
        "takeaway": (
            "Feature-level = X·(1+α·vimp); sample-level = membership w_i; both = product. "
            "Curves on FULL LIVE holdout (not first-n subset). No EWMA."
        ),
    }
