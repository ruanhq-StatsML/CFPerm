"""Wu et al. CVPR 2018 nonparametric instance discrimination, on FSDS blocks.

Frozen window embeddings play the role of the encoder. A linear projector is
trained with a memory bank: the current mini-batch is compared to stored
instance vectors, then those slots are momentum-updated. Extra keys from the
previous batch are concatenated as MoCo-style negatives, so batch t's bank
is the negative dictionary for batch t+1.

This module does not change TSS. It reads P(X) instance geometry, not P(Y|X).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from msrvtt_continuous_trainer import _softmax, _split_modalities
from msrvtt_multimodal_attribution import GROUP_NAMES, SEED

_HEAD_SEED = {"video": 1, "audio": 2, "text": 3}


def l2_normalize(X, axis=1):
    X = np.asarray(X, dtype=float)
    nrm = np.linalg.norm(X, axis=axis, keepdims=True)
    return X / (nrm + 1e-12), nrm


class MemoryBank:
    """One unit vector per instance. No gradient through the bank (Wu 2018)."""

    def __init__(self, n, dim, seed=SEED, momentum=0.5):
        rng = np.random.default_rng(seed)
        V = rng.normal(size=(int(n), int(dim)))
        self.V, _ = l2_normalize(V)
        self.momentum = float(momentum)

    def update(self, ids, f):
        ids = np.asarray(ids, dtype=int)
        f, _ = l2_normalize(f)
        v = self.momentum * self.V[ids] + (1.0 - self.momentum) * f
        self.V[ids], _ = l2_normalize(v)


class LinearProjector:
    """g: R^{d} → S^{k}, the only trained map on frozen X."""

    def __init__(self, d_in, d_out, seed=SEED):
        rng = np.random.default_rng(seed)
        scale = 1.0 / np.sqrt(max(int(d_in), 1))
        self.W = rng.normal(scale=scale, size=(int(d_in), int(d_out)))

    def encode(self, X):
        z = np.asarray(X, dtype=float) @ self.W
        f, nrm = l2_normalize(z)
        return f, z, nrm


def nonparametric_softmax_step(
    projector,
    bank,
    X,
    ids,
    tau=0.07,
    eta=0.05,
    extra_keys=None,
):
    """Wu nonparametric softmax. extra_keys are previous-batch embeddings."""
    X = np.asarray(X, dtype=float)
    ids = np.asarray(ids, dtype=int)
    f, _, nrm = projector.encode(X)
    keys = bank.V
    pos = ids
    if extra_keys is not None:
        extra = np.asarray(extra_keys, dtype=float)
        if extra.size:
            extra, _ = l2_normalize(extra)
            keys = np.vstack([bank.V, extra])
    logits = np.clip(f @ keys.T / float(tau), -50.0, 50.0)
    P = _softmax(logits)
    n = max(len(ids), 1)
    Glog = P / n
    Glog[np.arange(len(ids)), pos] -= 1.0 / n
    gf = (Glog @ keys) / float(tau)
    dots = (gf * f).sum(axis=1, keepdims=True)
    gz = (gf - f * dots) / (nrm + 1e-12)
    projector.W = projector.W - float(eta) * (X.T @ gz)
    bank.update(ids, f)
    p_pos = np.clip(P[np.arange(len(ids)), pos], 1e-8, 1.0)
    loss = float(-np.mean(np.log(p_pos)))
    acc = float(np.mean(logits.argmax(axis=1) == pos))
    return loss, acc


def instance_accuracy(projector, bank, X, ids):
    f, _, _ = projector.encode(X)
    pred = (f @ bank.V.T).argmax(axis=1)
    return float(np.mean(pred == np.asarray(ids, dtype=int)))


def max_sim_to_keys(f, keys):
    if keys is None or np.asarray(keys).size == 0:
        return float("nan")
    f, _ = l2_normalize(f)
    keys, _ = l2_normalize(keys)
    return float((f @ keys.T).max(axis=1).mean())


def mean_offdiag_cosine(X):
    Xn, _ = l2_normalize(np.asarray(X, dtype=float))
    S = Xn @ Xn.T
    n = S.shape[0]
    if n < 2:
        return float("nan")
    return float((S.sum() - np.trace(S)) / (n * (n - 1)))


@dataclass
class InstDiscState:
    projector: LinearProjector
    bank: MemoryBank
    keys: np.ndarray | None = None


def train_batch_instdisc(
    X,
    projector,
    bank,
    extra_keys=None,
    epochs=4,
    tau=0.07,
    eta=0.05,
    seed=SEED,
    chunk=None,
):
    X = np.asarray(X, dtype=float)
    n = X.shape[0]
    ids = np.arange(n, dtype=int)
    rng = np.random.default_rng(seed)
    chunk = int(chunk) if chunk is not None else max(8, n // 2)
    losses, accs = [], []
    for ep in range(int(epochs)):
        order = rng.permutation(n)
        el, ea = [], []
        for start in range(0, n, chunk):
            sl = order[start : start + chunk]
            loss, acc = nonparametric_softmax_step(
                projector,
                bank,
                X[sl],
                ids[sl],
                tau=tau,
                eta=eta,
                extra_keys=extra_keys,
            )
            el.append(loss)
            ea.append(acc)
        losses.append(float(np.mean(el)))
        accs.append(float(np.mean(ea)))
    f, _, _ = projector.encode(X)
    return {
        "loss": float(losses[-1]) if losses else float("nan"),
        "acc": float(accs[-1]) if accs else float("nan"),
        "acc_path": accs,
        "loss_path": losses,
        "features": f,
        "collision": mean_offdiag_cosine(X),
        "proj_collision": mean_offdiag_cosine(f),
    }


def run_instdisc_stream(
    stream,
    dim=32,
    epochs=4,
    tau=0.07,
    eta=0.08,
    seed=SEED,
):
    """Per-modality Wu bank. Previous-batch keys are next-batch negatives."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=int)
    batch = np.asarray(stream.batch, dtype=int)
    n_batches = int(batch.max()) + 1
    Xs_all = _split_modalities(X)
    states = {}
    history = []
    for t in range(n_batches):
        idx = np.flatnonzero(batch == t)
        row = {"round": int(t), "n": int(idx.size)}
        for g in GROUP_NAMES:
            Xm = Xs_all[g][idx]
            if t == 0:
                hs = _HEAD_SEED[g]
                proj = LinearProjector(Xm.shape[1], dim, seed=seed + hs)
                bank = MemoryBank(len(idx), dim, seed=seed + 17 + hs)
                extra = None
            else:
                st = states[g]
                if st.bank.V.shape[0] != len(idx):
                    bank = MemoryBank(len(idx), dim, seed=seed + t + _HEAD_SEED[g])
                else:
                    bank = st.bank
                proj = st.projector
                extra = st.keys
            rec = train_batch_instdisc(
                Xm,
                proj,
                bank,
                extra_keys=extra,
                epochs=epochs,
                tau=tau,
                eta=eta,
                seed=seed + 3 * t + _HEAD_SEED[g],
            )
            stale = max_sim_to_keys(rec["features"], extra)
            row[g] = {
                "loss": rec["loss"],
                "acc": rec["acc"],
                "stale_max_sim": stale,
                "collision": rec["collision"],
                "proj_collision": rec["proj_collision"],
            }
            states[g] = InstDiscState(projector=proj, bank=bank, keys=rec["features"].copy())
        history.append(row)
    adapt = history[1:] if len(history) > 1 else history
    summary = {
        "n_batches": n_batches,
        "dim": int(dim),
        "mean_acc": {g: float(np.mean([h[g]["acc"] for h in adapt])) for g in GROUP_NAMES},
        "mean_stale": {g: float(np.nanmean([h[g]["stale_max_sim"] for h in adapt])) for g in GROUP_NAMES},
        "mean_collision": {g: float(np.mean([h[g]["collision"] for h in history])) for g in GROUP_NAMES},
        "history": history,
        "meta": dict(getattr(stream, "meta", {}) or {}),
    }
    return summary, states
