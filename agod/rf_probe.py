"""Shallow RF probe + last-two consecutive-OOS hop gate.

Same probe as OnlineRFPerm (``n_estimators=20``, ``max_depth=4``).
``hop_fires`` is the last-two gate, *not* the Palm–Nagler CI.

sklearn is imported lazily so unit tests of ``hop_fires`` do not need it.
"""
from __future__ import annotations

import numpy as np


class ConstantProbe:
    """Degenerate classifier when a window has a single label."""

    def __init__(self, label: int):
        self.label = int(label)
        self.classes_ = np.asarray([self.label], dtype=int)

    def predict(self, X):
        n = len(np.asarray(X))
        return np.full(n, self.label, dtype=int)

    def predict_proba(self, X):
        n = len(np.asarray(X))
        return np.ones((n, 1), dtype=float)


def fit_online_rf(X, y, *, seed=0, task="acc"):
    """Same shallow RF as the IPTW stream probe."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).ravel()
    if str(task) == "acc":
        y = y.astype(int)
        classes = np.unique(y)
        if classes.size < 2:
            return ConstantProbe(int(classes[0]) if classes.size else 0)
        from sklearn.ensemble import RandomForestClassifier

        clf = RandomForestClassifier(
            n_estimators=20,
            max_depth=4,
            min_samples_leaf=2,
            random_state=int(seed),
            n_jobs=1,
        )
        clf.fit(X, y)
        return clf
    from sklearn.ensemble import RandomForestRegressor

    rf = RandomForestRegressor(
        n_estimators=20,
        max_depth=4,
        min_samples_leaf=3,
        random_state=int(seed),
        n_jobs=1,
    )
    rf.fit(X, y.astype(float))
    return rf


fit_online_probe = fit_online_rf


def predict_p1(model, X) -> np.ndarray:
    """P(Y=1 | X). Missing positive class → 0."""
    X = np.asarray(X, dtype=float)
    if hasattr(model, "predict_proba"):
        proba = np.asarray(model.predict_proba(X), dtype=float)
        classes = [int(c) for c in getattr(model, "classes_", [])]
        if 1 in classes:
            return proba[:, classes.index(1)]
        if classes == [0]:
            return np.zeros(len(X), dtype=float)
        return proba[:, -1]
    return np.asarray(model.predict(X), dtype=float).ravel()


def probe_err(model, X, y, task="mse"):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y).ravel()
    if str(task) == "acc":
        pred = np.asarray(model.predict(X)).astype(int).ravel()
        return 1.0 - float(np.mean(pred == y.astype(int)))
    pred = np.asarray(model.predict(X), dtype=float).ravel()
    return float(np.mean(np.abs(y.astype(float) - pred)))


def brier_score(model, X, y) -> float:
    """Mean squared error of P(Y=1). Smooth analog of the concept-board MSE."""
    p1 = predict_p1(model, X)
    y = np.asarray(y, dtype=float).ravel()
    return float(np.mean((y - p1) ** 2))


score_probe = probe_err


def error_floor(task, n):
    """Minimum reliable OOS denominator.

    Classification on a constant-label stretch has e_prev=0, so
    e_now/e_prev is 10^7-scale and γ is vacuous. Require at least one
    mistake (1/n) and 2% error before a ratio is a hop.
    """
    if str(task) == "acc":
        return max(1.0 / max(int(n), 1), 0.02)
    return 1e-8


def shift_ratio(e_now, e_prev, e_floor=0.0):
    denom = max(float(e_prev), float(e_floor or 0.0), 1e-8)
    return float(e_now) / denom


def hop_fires(e_now, e_prev, gate=1.5, e_floor=0.0):
    """Consecutive OOS gate. First hop / vacuous e_prev → quiet."""
    if e_prev is None:
        return False
    e_prev = float(e_prev)
    floor = float(e_floor or 0.0)
    if e_prev < floor:
        return False
    ratio = shift_ratio(e_now, e_prev, e_floor=floor)
    return bool(np.isfinite(ratio) and ratio >= float(gate))
