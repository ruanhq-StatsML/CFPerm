"""Modality-Specific Gap (MSG) via RF-Domain classifiers and CFPerm PO-risk.

g_m^{(t)} = Normalize( AUC_m * VIMP_m + gamma * PO-risk_m )

This is a distance-based statistic, not a causal claim: under the
alternative the usual ignorability conditions for CATE are not met, and
the meta-learner is used as a discrepancy estimator (see CFPerm README).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .chronoberg import MODALITIES, TimeWindowBatch

# Overlap below this triggers the paper's RF-Domain-only fallback.
OVERLAP_FALLBACK = 0.3


@dataclass
class MSGResult:
    modality: str
    auc: float
    vimp: float
    po_risk: float
    gap_raw: float
    overlap: float
    used_fallback: bool
    top_features: Tuple[int, ...] = ()
    feature_vimp: np.ndarray = field(default_factory=lambda: np.zeros(0))


@dataclass
class MSGState:
    gaps: Dict[str, float]
    details: Dict[str, MSGResult]

    def vector(self, modalities: Sequence[str] = MODALITIES) -> np.ndarray:
        return np.array([self.gaps[m] for m in modalities], dtype=np.float64)


def _as_2d(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    return x


def _safe_auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    y_true = np.asarray(y_true).astype(int)
    if np.unique(y_true).size < 2:
        return 0.5
    scores = np.asarray(scores, dtype=np.float64)
    if not np.isfinite(scores).all():
        return 0.5
    try:
        return float(roc_auc_score(y_true, scores))
    except ValueError:
        return 0.5


def rf_domain_auc_vimp(
    X0: np.ndarray,
    Xt: np.ndarray,
    *,
    n_estimators: int = 60,
    max_depth: int = 5,
    seed: int = 0,
    n_splits: int = 3,
) -> Tuple[float, float, np.ndarray, float, Tuple[int, ...]]:
    """Domain classifier separating D0 vs Dt.

    Returns ``(auc, vimp_scalar, per_feature_vimp, overlap, top_features)``.
    ``vimp_scalar`` grows with both separability ``(2*AUC-1)_+`` and how
    concentrated the impurity mass is on a few coordinates.
    """
    X0 = _as_2d(X0)
    Xt = _as_2d(Xt)
    X = np.vstack([X0, Xt])
    y = np.concatenate([np.zeros(len(X0), dtype=int), np.ones(len(Xt), dtype=int)])
    n, p = X.shape
    clf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_leaf=max(1, n // 40),
        random_state=seed,
        n_jobs=1,
    )
    n_splits = max(2, min(n_splits, int(y.sum()), int((1 - y).sum()), n // 4 or 2))
    if n_splits < 2 or min(int(y.sum()), int((1 - y).sum())) < n_splits:
        clf.fit(X, y)
        proba = clf.predict_proba(X)[:, 1]
        auc = _safe_auc(y, proba)
        vimp = np.asarray(clf.feature_importances_, dtype=np.float64)
    else:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        oob = np.zeros(n, dtype=np.float64)
        vimp_acc = np.zeros(p, dtype=np.float64)
        for fold, (tr, te) in enumerate(cv.split(X, y)):
            model = RandomForestClassifier(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_leaf=max(1, n // 40),
                random_state=seed + fold,
                n_jobs=1,
            )
            model.fit(X[tr], y[tr])
            oob[te] = model.predict_proba(X[te])[:, 1]
            vimp_acc += np.asarray(model.feature_importances_, dtype=np.float64)
        auc = _safe_auc(y, oob)
        vimp = vimp_acc / n_splits
        clf.fit(X, y)

    overlap = float(np.clip(1.0 - abs(2.0 * auc - 1.0), 0.0, 1.0))
    shift_mass = float(np.clip(2.0 * auc - 1.0, 0.0, 1.0))
    # Impurity importances sum to 1; peakiness distinguishes a few driving
    # coordinates from a diffuse, uninformative split.
    peakiness = float(np.max(vimp) * p) if p else 1.0
    vimp_scalar = float(shift_mass * (1.0 + peakiness))
    top = tuple(int(i) for i in np.argsort(-vimp)[: min(5, p)])
    return auc, vimp_scalar, vimp, overlap, top


def _crossfit_nuisance(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int,
    n_splits: int,
    clip_e: float,
) -> Tuple[np.ndarray, np.ndarray]:
    n = X.shape[0]
    m_hat = np.zeros(n, dtype=np.float64)
    e_hat = np.zeros(n, dtype=np.float64)
    n_splits = max(2, min(n_splits, int(W.sum()), int((1 - W).sum())))
    if n_splits < 2:
        ridge = Pipeline([("scaler", StandardScaler()), ("model", Ridge(alpha=1.0))])
        logit = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", LogisticRegression(max_iter=200, solver="lbfgs")),
            ]
        )
        ridge.fit(X, Y)
        logit.fit(X, W)
        m_hat[:] = ridge.predict(X)
        e_hat[:] = logit.predict_proba(X)[:, 1]
        return m_hat, np.clip(e_hat, clip_e, 1.0 - clip_e)

    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        y_model = RandomForestRegressor(
            n_estimators=40,
            max_depth=4,
            min_samples_leaf=2,
            random_state=seed + fold,
            n_jobs=1,
        )
        e_model = LogisticRegression(max_iter=300, solver="lbfgs")
        y_model.fit(X[tr], Y[tr])
        e_model.fit(X[tr], W[tr])
        m_hat[te] = y_model.predict(X[te])
        e_hat[te] = e_model.predict_proba(X[te])[:, 1]
    return m_hat, np.clip(e_hat, clip_e, 1.0 - clip_e)


def cfperm_po_risk(
    X0: np.ndarray,
    Y0: np.ndarray,
    Xt: np.ndarray,
    Yt: np.ndarray,
    *,
    seed: int = 0,
    n_splits: int = 3,
    clip_e: float = 0.05,
    clip_wtilde: float = 1e-3,
) -> float:
    """R-risk of a cross-fitted DR / R-learner treating batch as treatment.

    Large values mean the batch indicator still explains leftover outcome
    variation after adjusting for X — a proxy for concept drift. This is
    the online (no permutation) CFPerm statistic used as a clever covariate.
    """
    X0 = _as_2d(X0)
    Xt = _as_2d(Xt)
    X = np.vstack([X0, Xt])
    Y = np.concatenate([np.asarray(Y0, dtype=np.float64).ravel(), np.asarray(Yt, dtype=np.float64).ravel()])
    W = np.concatenate([np.zeros(len(X0), dtype=int), np.ones(len(Xt), dtype=int)])
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    m_hat, e_hat = _crossfit_nuisance(Xs, Y, W, seed=seed, n_splits=n_splits, clip_e=clip_e)
    y_tilde = Y - m_hat
    w_tilde = W.astype(np.float64) - e_hat
    denom = np.where(
        np.abs(w_tilde) < clip_wtilde,
        clip_wtilde * np.where(w_tilde >= 0, 1.0, -1.0),
        w_tilde,
    )
    z = y_tilde / denom
    weights = np.clip(w_tilde ** 2, 1e-6, None)
    tau_model = Ridge(alpha=1.0, random_state=seed)
    tau_model.fit(Xs, z, sample_weight=weights)
    tau_hat = tau_model.predict(Xs)
    # CFPerm uses treatment-effect heterogeneity, not a global intercept
    # shift: a mean change in Y induced by another modality looks like a
    # constant CATE and must not leak into every modality's PO-risk.
    hetero = float(np.var(tau_hat))
    return hetero


def _minmax(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    lo, hi = float(np.min(v)), float(np.max(v))
    if not np.isfinite(v).all() or hi - lo < 1e-12:
        return np.zeros_like(v)
    return (v - lo) / (hi - lo)


def compute_modality_gap(
    X0: np.ndarray,
    Y0: np.ndarray,
    Xt: np.ndarray,
    Yt: np.ndarray,
    *,
    modality: str,
    gamma: float = 1.0,
    seed: int = 0,
    overlap_threshold: float = OVERLAP_FALLBACK,
) -> MSGResult:
    auc, vimp_scalar, vimp_vec, overlap, top = rf_domain_auc_vimp(X0, Xt, seed=seed)
    po = cfperm_po_risk(X0, Y0, Xt, Yt, seed=seed)
    n0 = max(8, len(Y0) // 2)
    po_null = cfperm_po_risk(X0[:n0], Y0[:n0], X0[n0:], Y0[n0:], seed=seed + 1)
    po_lift = float(po / (po_null + 1e-6))
    used_fallback = overlap < overlap_threshold
    cov_term = float(auc * vimp_scalar)
    po_term = 0.0 if used_fallback else float(gamma * np.log1p(max(po_lift - 1.0, 0.0)))
    gap_raw = cov_term + po_term
    return MSGResult(
        modality=modality,
        auc=float(auc),
        vimp=float(vimp_scalar),
        po_risk=float(po_lift),
        gap_raw=float(gap_raw),
        overlap=float(overlap),
        used_fallback=bool(used_fallback),
        top_features=top,
        feature_vimp=vimp_vec,
    )


def compute_msg_state(
    reference: TimeWindowBatch,
    current: TimeWindowBatch,
    *,
    gamma: float = 1.0,
    seed: int = 0,
    modalities: Sequence[str] = MODALITIES,
    overlap_threshold: float = OVERLAP_FALLBACK,
) -> MSGState:
    details: Dict[str, MSGResult] = {}
    raw = []
    for i, m in enumerate(modalities):
        details[m] = compute_modality_gap(
            reference.X[m],
            reference.Y,
            current.X[m],
            current.Y,
            modality=m,
            gamma=gamma,
            seed=seed + 17 * i,
            overlap_threshold=overlap_threshold,
        )
        raw.append(details[m].gap_raw)
    gaps_vec = _minmax(np.asarray(raw, dtype=np.float64))
    # Keep a floor so a fully-stable modality is not hard-zeroed before softmax.
    gaps_vec = 0.05 + 0.95 * gaps_vec
    gaps = {m: float(gaps_vec[i]) for i, m in enumerate(modalities)}
    return MSGState(gaps=gaps, details=details)
