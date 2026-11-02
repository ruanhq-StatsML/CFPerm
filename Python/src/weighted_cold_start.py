"""Weighted cold-start under category shift (non-consistent classes).

Pipeline
--------
  detect (error-rise gate)
  → FSDS/CFPerm-style attribution (RF-domain VIMP) + subset localization
  → boundary check (label-support overlap → comparable | incomparable)
  → upweight localized subsets + selected features
  → sklearn partial_fit cold-start

Weights are *derived* from attribution (ranking / membership), not raw VIMP.
Incomparable batches skip directed adaptation (honest boundary).
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

SEED = 2026


# ---------------------------------------------------------------------------
# Detect → attribute → localize → boundary → weights → partial_fit
# ---------------------------------------------------------------------------


def detect_error_rise(
    clf,
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    X_live: np.ndarray,
    y_live: np.ndarray,
    *,
    thresh: float = 1.15,
) -> Dict:
    """OnlineRFPerm-style gate: LIVE/REF error ratio (labels for eval/gate only)."""
    X_ref, X_live = np.asarray(X_ref, float), np.asarray(X_live, float)
    y_ref, y_live = np.asarray(y_ref).ravel(), np.asarray(y_live).ravel()
    e_ref = 1.0 - float(accuracy_score(y_ref, clf.predict(X_ref)))
    e_live = 1.0 - float(accuracy_score(y_live, clf.predict(X_live)))
    ratio = e_live / max(e_ref, 1e-6)
    return {
        "e_ref": e_ref,
        "e_live": e_live,
        "ratio": float(ratio),
        "trigger": bool(ratio >= float(thresh)),
    }


def attribute_rf_domain(
    X_ref: np.ndarray,
    X_live: np.ndarray,
    *,
    seed: int = SEED,
    n_estimators: int = 80,
) -> np.ndarray:
    """FSDS RF-domain arm: RF predicts batch W∈{0,1} → feature_importances_."""
    X0, X1 = np.asarray(X_ref, float), np.asarray(X_live, float)
    X = np.vstack([X0, X1])
    W = np.array([0] * len(X0) + [1] * len(X1), dtype=int)
    rf = RandomForestClassifier(
        n_estimators=int(n_estimators), max_depth=6, min_samples_leaf=5,
        random_state=int(seed), n_jobs=1,
    )
    rf.fit(X, W)
    v = np.asarray(rf.feature_importances_, float)
    return v / (v.sum() + 1e-12)


def localize_subset(
    X_ref: np.ndarray,
    X_live: np.ndarray,
    vimp: np.ndarray,
    *,
    top_k: int = 8,
) -> Tuple[np.ndarray, np.ndarray]:
    """Z-score distance in top-K VIMP subspace → LIVE membership in [0, 1]."""
    vimp = np.asarray(vimp, float)
    k = int(min(max(1, top_k), len(vimp)))
    feat = np.argsort(-vimp)[:k]
    mu = np.asarray(X_ref, float)[:, feat].mean(0)
    sd = np.asarray(X_ref, float)[:, feat].std(0) + 1e-6
    z = (np.asarray(X_live, float)[:, feat] - mu) / sd
    d = np.linalg.norm(z, axis=1)
    # high distance in attributed subspace → more "new-category-like"
    lo, hi = np.quantile(d, 0.2), np.quantile(d, 0.8)
    mem = np.clip((d - lo) / max(hi - lo, 1e-6), 0.0, 1.0)
    return mem, feat


def boundary_check(
    labels_ref: Sequence,
    labels_live: Sequence,
    *,
    min_jaccard: float = 0.20,
) -> Dict:
    """Label-support overlap. Low Jaccard → incomparable (skip directed adapt)."""
    a = set(np.asarray(labels_ref).ravel().tolist())
    b = set(np.asarray(labels_live).ravel().tolist())
    inter, union = a & b, a | b
    j = len(inter) / max(len(union), 1)
    shared = sorted(inter)
    private_live = sorted(b - a)
    return {
        "jaccard": float(j),
        "n_shared": len(shared),
        "n_private_live": len(private_live),
        "shared": shared,
        "private_live": private_live,
        "comparable": bool(j >= float(min_jaccard)),
        "decision": "comparable" if j >= float(min_jaccard) else "incomparable",
    }


def weights_from_attribution(
    membership: np.ndarray,
    feat_idx: np.ndarray,
    p: int,
    *,
    lam: float = 2.0,
    cap: float = 5.0,
) -> Dict:
    """Sample weights from localized membership; feature mask from top-K."""
    mem = np.asarray(membership, float).ravel()
    w = np.clip(1.0 + float(lam) * mem, 1.0, float(cap))
    w = w * (len(w) / max(w.sum(), 1e-12))  # mean 1
    mask = np.zeros(int(p), dtype=bool)
    mask[np.asarray(feat_idx, int)] = True
    return {"sample_weight": w, "feature_mask": mask, "feat_idx": np.asarray(feat_idx, int)}


def make_sgd(*, seed: int = SEED) -> SGDClassifier:
    return SGDClassifier(
        loss="log_loss", penalty="l2", alpha=1e-4, random_state=int(seed),
        learning_rate="optimal",
    )


def cold_start_partial_fit(
    clf: SGDClassifier,
    X: np.ndarray,
    y: np.ndarray,
    *,
    sample_weight: Optional[np.ndarray] = None,
    feature_mask: Optional[np.ndarray] = None,
    classes: Optional[Sequence[int]] = None,
    n_passes: int = 1,
) -> SGDClassifier:
    """partial_fit with sample_weight; optional feature mask (zero out rest)."""
    X = np.asarray(X, float).copy()
    y = np.asarray(y).ravel()
    # Feature mask is used for localization / reporting; sample weights carry
    # the cold-start signal. Optional mild amplify of attributed coords only.
    if feature_mask is not None:
        m = np.asarray(feature_mask, bool)
        X[:, m] *= 1.15
    # sklearn requires the same `classes` vector on every partial_fit call
    if classes is None:
        classes = getattr(clf, "classes_", None)
        if classes is None:
            classes = np.unique(y)
    classes = np.asarray(classes, dtype=int).ravel()
    sw = None if sample_weight is None else np.asarray(sample_weight, float).ravel()
    for _ in range(max(1, int(n_passes))):
        clf.partial_fit(X, y, classes=classes, sample_weight=sw)
    return clf


def adapt_batch(
    clf: SGDClassifier,
    X_ref: np.ndarray,
    y_ref: np.ndarray,
    X_live: np.ndarray,
    y_live: np.ndarray,
    *,
    seed: int = SEED,
    top_k: int = 8,
    lam: float = 2.0,
    min_jaccard: float = 0.20,
    detect_thresh: float = 1.15,
    force_adapt: bool = False,
) -> Dict:
    """One stream step: detect → attribute → boundary → (maybe) weighted cold-start."""
    y_ref, y_live = np.asarray(y_ref).ravel(), np.asarray(y_live).ravel()
    gate = detect_error_rise(clf, X_ref, y_ref, X_live, y_live, thresh=detect_thresh)
    bound = boundary_check(y_ref, y_live, min_jaccard=min_jaccard)
    vimp = attribute_rf_domain(X_ref, X_live, seed=seed)
    mem, feat = localize_subset(X_ref, X_live, vimp, top_k=top_k)
    # New-category cold-start: fuse subspace membership with private-live labels
    # (pseudo-labels / discovered ids in the real pipeline).
    if bound["private_live"]:
        mem = np.maximum(mem, np.isin(y_live, bound["private_live"]).astype(float))
    wt = weights_from_attribution(mem, feat, X_live.shape[1], lam=lam)

    adapted = False
    reason = "no_trigger"
    if bound["decision"] == "incomparable" and not force_adapt:
        reason = "incomparable"
    elif gate["trigger"] or force_adapt:
        classes = getattr(clf, "classes_", None)
        if classes is None:
            classes = np.unique(np.concatenate([y_ref, y_live]))
        clf = cold_start_partial_fit(
            clf, X_live, y_live,
            sample_weight=wt["sample_weight"],
            feature_mask=wt["feature_mask"],
            classes=classes,
            n_passes=3,
        )
        adapted = True
        reason = "weighted_cold_start" if gate["trigger"] else "forced"
    return {
        "clf": clf,
        "gate": gate,
        "boundary": bound,
        "vimp": vimp,
        "membership": mem,
        "feat_idx": feat,
        "weights": wt,
        "adapted": adapted,
        "reason": reason,
        "vimp_top_mass": float(vimp[feat].sum()),
    }


# ---------------------------------------------------------------------------
# Baselines for ablations
# ---------------------------------------------------------------------------


def adapt_uniform_partial_fit(clf, X_live, y_live, classes):
    return cold_start_partial_fit(clf, X_live, y_live, classes=classes, n_passes=2)


def fit_initial(
    X, y, *, seed: int = SEED, classes: Optional[Sequence[int]] = None,
) -> Tuple[SGDClassifier, StandardScaler]:
    """First partial_fit. Pass full `classes` (incl. future novel) so C-GCD can expand."""
    X, y = np.asarray(X, float), np.asarray(y).ravel()
    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)
    classes = np.array(
        sorted(set(int(c) for c in (classes if classes is not None else y))), dtype=int
    )
    clf = make_sgd(seed=seed)
    clf.partial_fit(Xs, y, classes=classes)
    return clf, sc
