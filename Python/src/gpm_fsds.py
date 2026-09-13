"""GPM (Saha et al., ICLR 2021) on FSDS modality blocks, gated by TSS (c, δ).

After batch t the CGS of head m is the leading SVD bases of X^{(m)}.
The next batch's linear-head gradient is projected orthogonal to those
bases *only when* covariate intensity is loud and concept intensity is
quiet — the same thresholds as TSS. Concept-loud heads are not projected:
GPM as written would freeze P(Y|X) in the old span.

This is GPM composed with the locked typed signs, not a new product loop.
"""
from __future__ import annotations

import numpy as np

from msrvtt_continuous_trainer import SeparateHeadProbe, _minibatch, _softmax, _split_modalities
from msrvtt_multimodal_attribution import GROUP_NAMES, SEED
from typed_shift_stepsize import (
    ETA0,
    TAU_C,
    TAU_D,
    _acc,
    concept_intensity,
    covariate_intensity,
)


def representation_bases(X, thresh=0.97, max_k=None):
    """Uncentered SVD of a representation matrix (GPM CGS)."""
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return np.zeros((0, 0)), 0.0, 0
    _, S, Vt = np.linalg.svd(X, full_matrices=False)
    energy = S * S
    tot = float(energy.sum()) + 1e-12
    cdf = np.cumsum(energy) / tot
    k = int(np.searchsorted(cdf, float(thresh)) + 1)
    k = max(1, min(k, Vt.shape[0]))
    if max_k is not None:
        k = min(k, int(max_k))
    explained = float(cdf[k - 1]) if k else 0.0
    return Vt[:k].T.copy(), explained, k


def residual_energy(X, M):
    """Fraction of ||X||_F^2 outside span(M)."""
    X = np.asarray(X, dtype=float)
    den = float((X * X).sum()) + 1e-12
    if M is None or np.asarray(M).size == 0:
        return 1.0
    proj = X @ M @ M.T
    res = X - proj
    return float((res * res).sum() / den)


def gpm_extend(M, X, thresh=0.97, max_k=None):
    """Append new orthogonal bases of the residual representation (GPM)."""
    X = np.asarray(X, dtype=float)
    if M is None or np.asarray(M).size == 0:
        B, _, _ = representation_bases(X, thresh=thresh, max_k=max_k)
        return B
    M = np.asarray(M, dtype=float)
    X_res = X - X @ M @ M.T
    if float((X_res * X_res).sum()) < 1e-10:
        return M
    B, _, _ = representation_bases(X_res, thresh=thresh, max_k=max_k)
    return np.hstack([M, B])


def project_grad(dW, M):
    """ΔW ← (I − MM^⊤) ΔW, GPM orthogonal step on a linear head."""
    if M is None or np.asarray(M).size == 0:
        return dW
    M = np.asarray(M, dtype=float)
    return dW - M @ (M.T @ dW)


def gpm_gate(c, delta, tau_c=TAU_C, tau_d=TAU_D):
    """Project iff covariate loud and concept quiet. Same τ as TSS."""
    return float(c) >= float(tau_c) and float(delta) < float(tau_d)


def _head_grad(probe, Xs, y):
    y = np.asarray(y, dtype=int)
    P = _softmax(probe.logits(Xs))
    n = max(len(y), 1)
    G = P.copy()
    G[np.arange(len(y)), y] -= 1.0
    G /= n
    return G


def gpm_step(probe, Xs, y, lrs, bases, gate, active=None):
    G = _head_grad(probe, Xs, y)
    names = (active,) if active else GROUP_NAMES
    for name in names:
        dW = Xs[name].T @ G
        if gate.get(name, False):
            dW = project_grad(dW, bases.get(name))
        probe.W[name] = probe.W[name] - float(lrs.get(name, 0.0)) * dW


def run_gpm_method(
    stream,
    mode="typed",
    eta0=ETA0,
    steps_per_batch=6,
    warmup_steps=4,
    seed=SEED,
    thresh=0.90,
    max_k=8,
):
    """Linear probe with GPM projection. mode: none | always | typed."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=int)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    n_classes = int(y.max()) + 1
    probe = SeparateHeadProbe(n_classes=n_classes, seed=seed)
    rng = np.random.default_rng(seed)
    lrs = {g: float(eta0) / 3.0 for g in GROUP_NAMES}
    bases = {g: None for g in GROUP_NAMES}
    i0 = np.flatnonzero(batch == 0)
    Xs0, y0 = _split_modalities(X[i0]), y[i0]
    chunk = max(4, i0.size // 2)
    for head in GROUP_NAMES:
        for _ in range(int(max(1, warmup_steps))):
            sl = _minibatch(rng, i0.size, chunk)
            probe.step({g: Xs0[g][sl] for g in GROUP_NAMES}, y0[sl], lrs, active=head)
        bases[head] = gpm_extend(None, Xs0[head], thresh=thresh, max_k=max_k)

    history = []
    c_hat = {g: 0.0 for g in GROUP_NAMES}
    d_hat = {g: 0.0 for g in GROUP_NAMES}

    def snapshot(round_id, phase, idx, gate):
        Xs = _split_modalities(X[idx])
        return {
            "round": int(round_id),
            "phase": phase,
            "mode": mode,
            "acc": _acc(probe.logits(Xs), y[idx]),
            "bwt": _acc(probe.logits(Xs0), y0),
            "c": {g: float(c_hat[g]) for g in GROUP_NAMES},
            "delta": {g: float(d_hat[g]) for g in GROUP_NAMES},
            "gate": {g: bool(gate.get(g, False)) for g in GROUP_NAMES},
            "residual": {g: residual_energy(Xs[g], bases[g]) for g in GROUP_NAMES},
            "rank": {g: int(0 if bases[g] is None else bases[g].shape[1]) for g in GROUP_NAMES},
        }

    history.append(snapshot(0, "warmup", i0, {g: False for g in GROUP_NAMES}))

    for t in range(1, n_batches):
        ip, ic = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        Xsp, Xsc = _split_modalities(X[ip]), _split_modalities(X[ic])
        yp, yc = y[ip], y[ic]
        c_hat, _, _ = covariate_intensity(X[ip], X[ic])
        d_hat, _ = concept_intensity(probe, Xsp, yp, Xsc, yc)
        if mode == "none":
            gate = {g: False for g in GROUP_NAMES}
        elif mode == "always":
            gate = {g: True for g in GROUP_NAMES}
        else:
            gate = {g: gpm_gate(c_hat[g], d_hat[g]) for g in GROUP_NAMES}
        n_steps = max(1, int(steps_per_batch) // 3)
        chunk_t = max(4, ic.size // 2)
        for head in GROUP_NAMES:
            for _ in range(n_steps):
                sl = _minibatch(rng, ic.size, chunk_t)
                gpm_step(
                    probe,
                    {g: Xsc[g][sl] for g in GROUP_NAMES},
                    yc[sl],
                    lrs,
                    bases,
                    gate,
                    active=head,
                )
        history.append(snapshot(t, "adapt", ic, gate))
        for g in GROUP_NAMES:
            bases[g] = gpm_extend(bases[g], Xsc[g], thresh=thresh, max_k=max_k)

    adapt = [h for h in history if h["phase"] == "adapt"]
    split = int(stream.meta.get("concept_at") or max(n_batches // 2, 1))
    post = [h for h in adapt if h["round"] >= split]
    return {
        "mode": mode,
        "n_batches": n_batches,
        "online_acc": float(np.mean([h["acc"] for h in adapt])) if adapt else float("nan"),
        "post_acc": float(np.mean([h["acc"] for h in post])) if post else float("nan"),
        "bwt": history[-1]["bwt"] if history else float("nan"),
        "mean_gate": {g: float(np.mean([h["gate"][g] for h in adapt])) for g in GROUP_NAMES} if adapt else {},
        "mean_residual": {g: float(np.mean([h["residual"][g] for h in adapt])) for g in GROUP_NAMES} if adapt else {},
        "history": history,
        "meta": dict(stream.meta),
    }


def gpm_identification(stream, thresh=0.90, max_k=8):
    """Residual energy of batch t outside CGS of batch t-1, vs FSDS c_m."""
    X = np.asarray(stream.X, dtype=float)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    history = []
    prev = None
    bases = {g: None for g in GROUP_NAMES}
    for t in range(n_batches):
        idx = np.flatnonzero(batch == t)
        Xs = _split_modalities(X[idx])
        row = {"round": int(t)}
        c_hat = {g: 0.0 for g in GROUP_NAMES}
        if prev is not None:
            c_hat, share, _ = covariate_intensity(X[prev], X[idx])
            row["c"] = c_hat
            row["pi"] = share
        else:
            row["c"] = c_hat
            row["pi"] = {g: 1.0 / 3.0 for g in GROUP_NAMES}
        row["residual"] = {}
        row["rank"] = {}
        for g in GROUP_NAMES:
            row["residual"][g] = residual_energy(Xs[g], bases[g])
            bases[g] = gpm_extend(bases[g], Xs[g], thresh=thresh, max_k=max_k)
            row["rank"][g] = int(bases[g].shape[1])
        history.append(row)
        prev = idx
    adapt = history[1:] if len(history) > 1 else history
    return {
        "mean_residual": {g: float(np.mean([h["residual"][g] for h in adapt])) for g in GROUP_NAMES},
        "mean_c": {g: float(np.mean([h["c"][g] for h in adapt])) for g in GROUP_NAMES},
        "mean_rank": {g: float(np.mean([h["rank"][g] for h in adapt])) for g in GROUP_NAMES},
        "history": history,
        "meta": dict(getattr(stream, "meta", {}) or {}),
    }
