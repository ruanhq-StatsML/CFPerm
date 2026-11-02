"""Weighted cold-start under category shift (non-consistent classes).

Cold-start logic (locked)
-------------------------
  1. Detect error rise (gate).
  2. FSDS RF-domain VIMP → top-K subspace → membership (diagnostic).
  3. Boundary: low label-support Jaccard → incomparable → *no* directed adapt.
  4. Fuse membership with private-live labels (discovered / pseudo ids).
  5. Build weights:
       sample_w  = 1 + λ · membership_fused     (soft)
       subset_w  = boost on hard observation subset (membership ≥ q)
       fused_w   = sample_w × subset_w          ← default cold-start channel
  6. partial_fit with fused_w (+ mild amplify of attributed features).

Pure membership alone is for diagnosis / ablation. Digits/amazon show weights
are not unconditionally helpful — boundary + sparsity still gate the update.
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


def sample_weights(membership: np.ndarray, *, lam: float = 2.5, cap: float = 5.0) -> np.ndarray:
    """Soft per-observation weights from membership (mean-normalized)."""
    mem = np.asarray(membership, float).ravel()
    w = np.clip(1.0 + float(lam) * mem, 1.0, float(cap))
    return w * (len(w) / max(w.sum(), 1e-12))


def subset_weights(
    membership: np.ndarray,
    *,
    q: float = 0.6,
    boost: float = 3.0,
    max_frac: float = 0.45,
) -> Tuple[np.ndarray, np.ndarray]:
    """Hard observation-subset weights: boost the localized high-membership set."""
    mem = np.asarray(membership, float).ravel()
    n = len(mem)
    thr = float(np.quantile(mem, float(q)))
    mask = mem >= thr
    if not mask.any() or mask.mean() > float(max_frac):
        k = max(1, int(round(n * min(1.0 - float(q), float(max_frac)))))
        mask = np.zeros(n, dtype=bool)
        mask[np.argsort(-mem)[:k]] = True
    w = np.ones(n, dtype=float)
    w[mask] = float(boost)
    w = w * (n / max(w.sum(), 1e-12))
    return w, mask


def combine_weights(sample_w: np.ndarray, subset_w: np.ndarray) -> np.ndarray:
    w = np.asarray(sample_w, float) * np.asarray(subset_w, float)
    return w * (len(w) / max(w.sum(), 1e-12))


def weights_from_attribution(
    membership: np.ndarray,
    feat_idx: np.ndarray,
    p: int,
    *,
    lam: float = 2.5,
    cap: float = 5.0,
    subset_q: float = 0.6,
    subset_boost: float = 3.0,
    mode: str = "fused",
) -> Dict:
    """Build cold-start weights from membership.

    mode
      sample   — soft sample_w only
      subset   — hard observation subset_w only
      fused    — sample_w × subset_w  (default cold-start channel)
    """
    mem = np.asarray(membership, float).ravel()
    sw = sample_weights(mem, lam=lam, cap=cap)
    sub_w, sub_mask = subset_weights(mem, q=subset_q, boost=subset_boost)
    mode = str(mode).lower()
    if mode == "sample":
        w = sw
    elif mode == "subset":
        w = sub_w
    else:
        w = combine_weights(sw, sub_w)
        mode = "fused"
    fmask = np.zeros(int(p), dtype=bool)
    fmask[np.asarray(feat_idx, int)] = True
    return {
        "sample_weight": w,
        "sample_w": sw,
        "subset_w": sub_w,
        "subset_mask": sub_mask,
        "feature_mask": fmask,
        "feat_idx": np.asarray(feat_idx, int),
        "weight_mode": mode,
        "subset_frac": float(sub_mask.mean()),
    }


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
    lam: float = 2.5,
    min_jaccard: float = 0.20,
    detect_thresh: float = 1.15,
    force_adapt: bool = False,
    use_private_labels: bool = True,
    weight_mode: str = "fused",
    subset_q: float = 0.6,
    subset_boost: float = 3.0,
) -> Dict:
    """Detect → attribute → boundary → fused sample_w×subset_w cold-start.

    Default weight_mode='fused' (locked cold-start channel).
    Pure membership (use_private_labels=False, weight_mode='sample') is diagnostic.
    """
    y_ref, y_live = np.asarray(y_ref).ravel(), np.asarray(y_live).ravel()
    gate = detect_error_rise(clf, X_ref, y_ref, X_live, y_live, thresh=detect_thresh)
    bound = boundary_check(y_ref, y_live, min_jaccard=min_jaccard)
    vimp = attribute_rf_domain(X_ref, X_live, seed=seed)
    mem_raw, feat = localize_subset(X_ref, X_live, vimp, top_k=top_k)
    mem = mem_raw.copy()
    if use_private_labels and bound["private_live"]:
        mem = np.maximum(mem, np.isin(y_live, bound["private_live"]).astype(float))
    wt = weights_from_attribution(
        mem, feat, X_live.shape[1], lam=lam,
        subset_q=subset_q, subset_boost=subset_boost, mode=weight_mode,
    )

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
        reason = "fused_cold_start" if wt["weight_mode"] == "fused" else "weighted_cold_start"
        if not gate["trigger"]:
            reason = "forced_" + reason
    return {
        "clf": clf,
        "gate": gate,
        "boundary": bound,
        "vimp": vimp,
        "membership": mem,
        "membership_raw": mem_raw,
        "feat_idx": feat,
        "weights": wt,
        "adapted": adapted,
        "reason": reason,
        "vimp_top_mass": float(vimp[feat].sum()),
        "use_private_labels": bool(use_private_labels),
        "weight_mode": wt["weight_mode"],
        "subset_frac": wt["subset_frac"],
        "mean_membership": float(np.mean(mem)),
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
