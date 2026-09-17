"""Streaming PO-risk: new batch is T=1, reference is T=0.

φ = (Y − μ)(T − e), τ̂(X) ≈ φ, risk = mean(τ̂²).

μ is a linear probe on the table, or the AnyMLP prediction
for that freeze-depth. Read the number; do not bootstrap.
"""
from __future__ import annotations

import numpy as np

CLIP = 1e-3
REF_N = 10_000
# Incoming T=1 block. Smaller than this, PO-risk jitters and you would
# need online-bootstrap on k+1 models — too expensive for the board.
MIN_STREAM_N = 5_000


def zscore(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - X.mean(axis=0)) / sd


def pack_ref_new(X_ref, Y_ref, X_new, Y_new):
    """Concatenate reference (T=0) and the incoming batch (T=1)."""
    X_ref = np.asarray(X_ref, dtype=float)
    X_new = np.asarray(X_new, dtype=float)
    Y = np.concatenate(
        [np.asarray(Y_ref, dtype=float).ravel(), np.asarray(Y_new, dtype=float).ravel()]
    )
    T = np.concatenate([np.zeros(len(X_ref)), np.ones(len(X_new))])
    X = np.vstack([X_ref, X_new])
    return X, Y, T


def po_risk(X, Y, T, mu=None, clip: float = CLIP) -> float:
    """Linear-probe PO-risk. Incoming batch must already be coded T=1."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    T = np.asarray(T, dtype=float).ravel()
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    Xd = np.column_stack([np.ones(len(X)), zscore(X)])
    if mu is None:
        beta, *_ = np.linalg.lstsq(Xd, Y, rcond=None)
        mu = Xd @ beta
    else:
        mu = np.asarray(mu, dtype=float).ravel()
    e_beta, *_ = np.linalg.lstsq(Xd, T, rcond=None)
    e_hat = np.clip(Xd @ e_beta, clip, 1.0 - clip)
    phi = (Y - mu) * (T - e_hat)
    tau, *_ = np.linalg.lstsq(Xd, phi, rcond=None)
    tau_hat = Xd @ tau
    return float(np.mean(tau_hat**2))


def streaming_po_risk(X_ref, Y_ref, X_new, Y_new, mu_fn=None, clip: float = CLIP) -> float:
    """PO-risk on (ref ∪ new) with T=1 on the new batch."""
    X, Y, T = pack_ref_new(X_ref, Y_ref, X_new, Y_new)
    mu = None if mu_fn is None else np.asarray(mu_fn(X), dtype=float).ravel()
    return po_risk(X, Y, T, mu=mu, clip=clip)
