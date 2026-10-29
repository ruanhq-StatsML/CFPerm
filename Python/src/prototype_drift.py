"""iCaRL NCM prototypes and Yu et al. CVPR 2020 semantic drift compensation.

Prototypes are class means in a modality block (or projector space). SDC
estimates the unknown drift of an old prototype from the observed drift of
current-batch points, then the compensated prototype is what the next batch
would use as an NCM classifier.

Not a new loss. Not wired into TSS.
"""
from __future__ import annotations

import numpy as np

from instance_discrimination import l2_normalize
from msrvtt_continuous_trainer import _split_modalities
from msrvtt_multimodal_attribution import GROUP_NAMES, SEED


def class_prototypes(X, y, n_classes=None):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=int)
    if n_classes is None:
        n_classes = int(y.max()) + 1 if y.size else 0
    mu = np.zeros((int(n_classes), X.shape[1]))
    counts = np.zeros(int(n_classes), dtype=int)
    for c in np.unique(y):
        sl = y == c
        mu[c] = X[sl].mean(axis=0)
        counts[c] = int(sl.sum())
    return mu, counts


def ncm_predict(X, mu):
    Xn, _ = l2_normalize(X)
    Mn, _ = l2_normalize(mu)
    return (Xn @ Mn.T).argmax(axis=1)


def ncm_accuracy(X, y, mu):
    y = np.asarray(y, dtype=int)
    if y.size == 0:
        return float("nan")
    return float(np.mean(ncm_predict(X, mu) == y))


def prototype_cosine(a, b):
    an, _ = l2_normalize(np.asarray(a, dtype=float))
    bn, _ = l2_normalize(np.asarray(b, dtype=float))
    return float(np.mean((an * bn).sum(axis=1)))


def sdc_compensate(mu_old, z_prev, z_curr, sigma=None):
    """Yu CVPR 2020, Eq. of interpolating δ_i onto μ_c.

    z_prev, z_curr are the *same* current-task points before and after the
    observed move (encoder/projector step, or a known location shift).
    """
    mu_old = np.asarray(mu_old, dtype=float)
    z_prev = np.asarray(z_prev, dtype=float)
    z_curr = np.asarray(z_curr, dtype=float)
    delta = z_curr - z_prev
    if sigma is None:
        dif = z_prev - z_prev.mean(axis=0, keepdims=True)
        sigma = float(np.median(np.sqrt((dif * dif).sum(axis=1))))
        sigma = max(sigma, 1e-3)
    out = mu_old.copy()
    for c in range(mu_old.shape[0]):
        d2 = ((z_prev - mu_old[c]) ** 2).sum(axis=1)
        w = np.exp(-d2 / (2.0 * sigma * sigma + 1e-12))
        w = w / (w.sum() + 1e-12)
        out[c] = mu_old[c] + w @ delta
    return out


def location_shift_field(X_prev, X_curr):
    """Constant vector field μ_t − μ_{t-1}. Exact for a pure location shift."""
    return np.asarray(X_curr, dtype=float).mean(axis=0) - np.asarray(X_prev, dtype=float).mean(axis=0)


def run_prototype_stream(stream, seed=SEED):
    """iCaRL NCM from batch t-1 vs SDC-compensated prototypes on batch t."""
    del seed
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=int)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    n_classes = int(y.max()) + 1
    Xs_all = _split_modalities(X)
    history = []
    prev_mu = {}
    prev_X = {}
    for t in range(n_batches):
        idx = np.flatnonzero(batch == t)
        row = {"round": int(t), "n": int(idx.size)}
        for g in GROUP_NAMES:
            Xm = Xs_all[g][idx]
            yt = y[idx]
            mu, _ = class_prototypes(Xm, yt, n_classes=n_classes)
            rec = {
                "ncm_self": ncm_accuracy(Xm, yt, mu),
                "ncm_stale": float("nan"),
                "ncm_sdc": float("nan"),
                "proto_cos_stale": float("nan"),
                "proto_cos_sdc": float("nan"),
            }
            if t > 0:
                rec["ncm_stale"] = ncm_accuracy(Xm, yt, prev_mu[g])
                field = location_shift_field(prev_X[g], Xm)
                z_prev = Xm - field
                mu_sdc = sdc_compensate(prev_mu[g], z_prev, Xm)
                rec["ncm_sdc"] = ncm_accuracy(Xm, yt, mu_sdc)
                rec["proto_cos_stale"] = prototype_cosine(prev_mu[g], mu)
                rec["proto_cos_sdc"] = prototype_cosine(mu_sdc, mu)
            row[g] = rec
            prev_mu[g] = mu
            prev_X[g] = Xm
        history.append(row)
    adapt = history[1:] if len(history) > 1 else history

    def _mean(key):
        return {g: float(np.nanmean([h[g][key] for h in adapt])) for g in GROUP_NAMES}

    return {
        "n_batches": n_batches,
        "mean_ncm_self": {g: float(np.mean([h[g]["ncm_self"] for h in history])) for g in GROUP_NAMES},
        "mean_ncm_stale": _mean("ncm_stale"),
        "mean_ncm_sdc": _mean("ncm_sdc"),
        "mean_proto_cos_stale": _mean("proto_cos_stale"),
        "mean_proto_cos_sdc": _mean("proto_cos_sdc"),
        "history": history,
        "meta": dict(getattr(stream, "meta", {}) or {}),
    }
