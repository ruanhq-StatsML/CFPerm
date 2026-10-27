"""Doubly-robust pseudo-outcome learner for a batch split.

This is a distance / shift estimator on (X, Y, W), not a CATE claim.
W is a batch index (e.g. early vs late timestamp), not a treatment.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler


def _as_1d(x) -> np.ndarray:
    return np.asarray(x).reshape(-1)


def dr_pseudo_outcome(
    y: np.ndarray,
    w: np.ndarray,
    mu0: np.ndarray,
    mu1: np.ndarray,
    e: np.ndarray,
    clip_e: float = 0.01,
) -> np.ndarray:
    e = np.clip(e, clip_e, 1.0 - clip_e)
    w = w.astype(np.float64)
    y = y.astype(np.float64)
    return (mu1 - mu0) + (w / e) * (y - mu1) - ((1.0 - w) / (1.0 - e)) * (y - mu0)


def _fit_mu_split(X: np.ndarray, y: np.ndarray, w: np.ndarray, alpha: float, seed: int):
    mu0 = Ridge(alpha=alpha)
    mu1 = Ridge(alpha=alpha)
    m0 = w == 0
    m1 = w == 1
    if m0.sum() < 2:
        mu0 = None
    else:
        mu0.fit(X[m0], y[m0])
    if m1.sum() < 2:
        mu1 = None
    else:
        mu1.fit(X[m1], y[m1])
    return mu0, mu1


def _predict_mu(model, X: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    if model is None:
        return np.full(X.shape[0], float(np.mean(fallback)), dtype=np.float64)
    return np.asarray(model.predict(X), dtype=np.float64).reshape(-1)


def fit_dr_pseudo_outcome(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    n_splits: int = 3,
    clip_e: float = 0.01,
    ridge_alpha: float = 1.0,
    seed: int = 0,
    max_iter: int = 200,
) -> Dict[str, Any]:
    X = np.asarray(X, dtype=np.float64)
    Y = _as_1d(Y).astype(np.float64)
    W = _as_1d(W).astype(int)
    n, p = X.shape
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)

    mu0_hat = np.zeros(n, dtype=np.float64)
    mu1_hat = np.zeros(n, dtype=np.float64)
    e_hat = np.zeros(n, dtype=np.float64)

    for fold, (tr, te) in enumerate(kf.split(Xs)):
        Xtr, Xte = Xs[tr], Xs[te]
        ytr, wtr = Y[tr], W[tr]
        mu0, mu1 = _fit_mu_split(Xtr, ytr, wtr, ridge_alpha, seed + fold)
        mu0_hat[te] = _predict_mu(mu0, Xte, ytr[wtr == 0] if (wtr == 0).any() else ytr)
        mu1_hat[te] = _predict_mu(mu1, Xte, ytr[wtr == 1] if (wtr == 1).any() else ytr)
        clf = LogisticRegression(max_iter=max_iter, solver="liblinear", random_state=seed + 100 + fold)
        if len(np.unique(wtr)) < 2:
            e_hat[te] = float(np.mean(wtr))
        else:
            clf.fit(Xtr, wtr)
            e_hat[te] = clf.predict_proba(Xte)[:, 1]

    e_hat = np.clip(e_hat, clip_e, 1.0 - clip_e)
    phi = dr_pseudo_outcome(Y, W, mu0_hat, mu1_hat, e_hat, clip_e=clip_e)

    tau_model = Ridge(alpha=ridge_alpha)
    tau_model.fit(Xs, phi)
    tau_hat = np.asarray(tau_model.predict(Xs), dtype=np.float64).reshape(-1)

    out = {
        "n": int(n),
        "p": int(p),
        "n_splits": int(n_splits),
        "clip_e": float(clip_e),
        "ridge_alpha": float(ridge_alpha),
        "seed": int(seed),
        "y_mean": float(Y.mean()),
        "w_mean": float(W.mean()),
        "naive_ate": float(Y[W == 1].mean() - Y[W == 0].mean()) if ((W == 0).any() and (W == 1).any()) else 0.0,
        "phi_mean": float(phi.mean()),
        "po_risk_phi": float(np.mean(phi ** 2)),
        "po_risk_tau": float(np.mean(tau_hat ** 2)),
        "tau_mean": float(tau_hat.mean()),
        "e_mean": float(e_hat.mean()),
        "mu0_mean": float(mu0_hat.mean()),
        "mu1_mean": float(mu1_hat.mean()),
        "note": (
            "DR pseudo-outcome on a batch split. "
            "Not a treatment-effect estimate; ignorability is not claimed. "
            "PO-risk is E[phi^2] and E[tau(X)^2] as a shift distance."
        ),
        "phi": phi,
        "tau_hat": tau_hat,
        "e_hat": e_hat,
        "mu0_hat": mu0_hat,
        "mu1_hat": mu1_hat,
        "scaler_mean": scaler.mean_,
        "scaler_scale": scaler.scale_,
        "tau_coef": np.asarray(tau_model.coef_, dtype=np.float64).reshape(-1),
    }
    return out


def compact_blocks(X: np.ndarray, user_dim: int, emb_dim: int, keep: int = 32) -> np.ndarray:
    """Keep user coords plus the head/tail of each 1056-d mm block (82 then 84)."""
    hist = X[:, user_dim : user_dim + emb_dim]
    tgt = X[:, user_dim + emb_dim :]
    k = min(keep, emb_dim)
    hist_c = np.concatenate([hist[:, :k], hist[:, -min(32, emb_dim) :]], axis=1)
    tgt_c = np.concatenate([tgt[:, :k], tgt[:, -min(32, emb_dim) :]], axis=1)
    return np.concatenate([X[:, :user_dim], hist_c, tgt_c], axis=1)
