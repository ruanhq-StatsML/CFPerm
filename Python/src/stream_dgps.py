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
    """P(Y|X) is a radial bump on the first 4 coords. Mean of those coords walks.

    Same f before and after. After onset the mass leaves the bump, so a model
    fit on D_ref pays MSE while PO-risk (same f) stays quieter than concept.
    """
    rng = np.random.default_rng(seed)
    n = int(n_ref) + int(n_new) * int(n_batches)
    X0 = rng.normal(size=(n, p))
    delta = np.zeros(p)
    delta[:4] = float(shift)
    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    X = X0 + a[:, None] * delta
    r = np.sqrt(np.sum(X[:, :4] ** 2, axis=1))
    Y = rng.binomial(1, _sigmoid(2.2 - 1.6 * r)).astype(float)
    meta = {
        "kind": "covariate",
        "onset_batch": int(onset_batch),
        "title": f"DGP gradual covariate (bump; μ walks after batch {onset_batch})",
        "n_ref": int(n_ref),
        "n_new": int(n_new),
    }
    return X, Y, meta


TRIMODAL_GROUPS = {
    "video": slice(0, 4),
    "audio": slice(4, 8),
    "text": slice(8, 12),
}


def make_trimodal_stream(
    n_ref: int = 800,
    n_new: int = 200,
    n_batches: int = 6,
    onset_batch: int = ONSET_BATCH,
    kind: str = "concept_video",
    n_slices: int = 4,
    shift: float = 2.2,
    seed: int = 2026,
):
    """Three feature blocks standing in for video / audio / text.

    Y depends on all three blocks. After onset:
      concept_video — flip the video coefficients; P(X) fixed
      covariate_audio — walk the audio mean; same f
      both — both of the above
    Only slice == (n_slices-1) is shifted, so subset tables have a planted peak.
    """
    rng = np.random.default_rng(seed)
    n = int(n_ref) + int(n_new) * int(n_batches)
    p = 12
    X = rng.normal(size=(n, p))
    n_slices = max(2, int(n_slices))
    slice_id = rng.integers(0, n_slices, size=n)
    a = _alpha(n_ref, n_new, n_batches, onset_batch)
    hot = (slice_id == (n_slices - 1)).astype(float)
    w = a * hot

    kind = str(kind)
    if kind in ("covariate_audio", "both"):
        X[:, 4:8] = X[:, 4:8] + w[:, None] * float(shift)

    beta_v = np.array([1.2, 0.4, 0.0, 0.0])
    beta_a = np.array([0.45, 0.0, 0.0, 0.0])
    beta_t = np.array([0.45, 0.0, 0.0, 0.0])
    logit0 = X[:, 0:4] @ beta_v + X[:, 4:8] @ beta_a + X[:, 8:12] @ beta_t
    if kind in ("concept_video", "both"):
        logit = logit0 * (1.0 - w) + (
            X[:, 0:4] @ (-beta_v) + X[:, 4:8] @ beta_a + X[:, 8:12] @ beta_t
        ) * w
    else:
        logit = logit0
    Y = rng.binomial(1, _sigmoid(logit)).astype(float)
    meta = {
        "kind": kind,
        "onset_batch": int(onset_batch),
        "groups": {k: list(range(s.start, s.stop)) for k, s in TRIMODAL_GROUPS.items()},
        "shifted_slice": int(n_slices - 1),
        "n_slices": int(n_slices),
        "title": f"trimodal {kind} (slice {n_slices - 1} after batch {onset_batch})",
        "n_ref": int(n_ref),
        "n_new": int(n_new),
    }
    return X, Y, slice_id, meta
