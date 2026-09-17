"""Gradual concept-drift vs gradual covariate-shift DGPs.

Concept: P(X) fixed, P(Y|X) rotates after a labeled onset batch.
Covariate: P(Y|X) fixed, mean of X walks after the same onset.

Labeled onset is a stream-batch index (0-based, after D_ref).
OnlineRFPerm marks WHEN. PO × MSE × MMD vs X_ref marks WHAT.
"""
from __future__ import annotations

import numpy as np

ONSET_BATCH = 2


def _sigmoid(z):
    z = np.clip(np.asarray(z, dtype=float), -30.0, 30.0)
    return 1.0 / (1.0 + np.exp(-z))


def _alpha(n_ref: int, n_new: int, n_batches: int, onset_batch: int) -> np.ndarray:
    n = int(n_ref) + int(n_new) * int(n_batches)
    idx = np.arange(n)
    start = int(n_ref) + int(onset_batch) * int(n_new)
    span = max(n - start - 1, 1)
    a = np.clip((idx - start) / span, 0.0, 1.0)
    a[idx < start] = 0.0
    return a


def make_gradual_concept(
    n_ref: int = 10_000,
    n_new: int = 2_000,
    n_batches: int = 12,
    onset_batch: int = ONSET_BATCH,
    p: int = 8,
    seed: int = 2026,
):
    """P(X)~N(0,I). β rotates from β0 to −β0 after onset_batch."""
    rng = np.random.default_rng(seed)
    n = int(n_ref) + int(n_new) * int(n_batches)
    X = rng.normal(size=(n, p))
    beta0 = np.zeros(p)
    beta0[:4] = 1.2
    beta1 = -beta0
    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    logit = (X @ beta0) * (1.0 - a) + (X @ beta1) * a
    Y = rng.binomial(1, _sigmoid(logit)).astype(float)
    meta = {
        "kind": "concept",
        "onset_batch": int(onset_batch),
        "title": f"DGP gradual concept (β flips after batch {onset_batch})",
        "n_ref": int(n_ref),
        "n_new": int(n_new),
    }
    return X, Y, meta


def make_gradual_covariate(
    n_ref: int = 10_000,
    n_new: int = 2_000,
    n_batches: int = 12,
    onset_batch: int = ONSET_BATCH,
    p: int = 8,
    shift: float = 2.5,
    seed: int = 2026,
):
    """P(Y|X) fixed. First 4 coordinates of X walk by `shift` after onset."""
    rng = np.random.default_rng(seed)
    n = int(n_ref) + int(n_new) * int(n_batches)
    X0 = rng.normal(size=(n, p))
    delta = np.zeros(p)
    delta[:4] = float(shift)
    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    X = X0 + a[:, None] * delta
    beta = np.zeros(p)
    beta[:4] = 1.2
    Y = rng.binomial(1, _sigmoid(X @ beta)).astype(float)
    meta = {
        "kind": "covariate",
        "onset_batch": int(onset_batch),
        "title": f"DGP gradual covariate (μ walks after batch {onset_batch})",
        "n_ref": int(n_ref),
        "n_new": int(n_new),
    }
    return X, Y, meta
