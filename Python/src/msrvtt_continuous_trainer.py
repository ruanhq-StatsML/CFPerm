"""Prototype continuous trainer: one head per modality, then rotate.

  1. Add a linear head for video / audio / text.
  2. Warmup: train each head once (the other heads stay frozen).
  3. Rotate: cycle video → audio → text; only the active head gets SGD.

Signed LR (locked):
  covariate shift c_m large  →  η_m down   (do not chase P(X|W))
  concept drift   δ_m large  →  η_m up     (relearn P(Y|X))
  both quiet (clip-level text) → η_m ≈ 0   (do not invert π)

The formal stepsize is ``typed_shift_stepsize`` (TSS): covariate intensity
lowers η_m, concept intensity raises it, both-quiet freezes the head.
This file keeps the warmup-then-rotate loop and the η∝π budget path as
a comparator, not as the signed map.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from msrvtt_multimodal_attribution import (
    GROUP_NAMES,
    GROUPS,
    SEED,
    assign_temporal_batches,
    modality_mass,
    benchmark_feature_selection,
    standardize_columns,
)


def _softmax(z):
    z = np.asarray(z, dtype=float)
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(np.clip(z, -50, 50))
    return e / (e.sum(axis=1, keepdims=True) + 1e-12)


def _split_modalities(X):
    return {name: np.asarray(X[:, sl], dtype=float) for name, sl in GROUPS.items()}


@dataclass
class SeparateHeadProbe:
    """Late-sum of three linear heads. SGD updates one ``active`` head at a time."""

    n_classes: int
    seed: int = SEED
    W: dict = field(default_factory=dict)

    def __post_init__(self):
        rng = np.random.default_rng(self.seed)
        if not self.W:
            self.W = {
                name: rng.normal(scale=0.01, size=(sl.stop - sl.start, self.n_classes))
                for name, sl in GROUPS.items()
            }

    def logits(self, Xs):
        z = None
        for name in GROUP_NAMES:
            term = Xs[name] @ self.W[name]
            z = term if z is None else z + term
        return z

    def step(self, Xs, y, lrs, active=None):
        """Update only ``active``; pass a name during warmup and rotate."""
        y = np.asarray(y, dtype=int)
        P = _softmax(self.logits(Xs))
        n = len(y)
        G = P.copy()
        G[np.arange(n), y] -= 1.0
        G /= max(n, 1)
        names = (active,) if active else GROUP_NAMES
        for name in names:
            self.W[name] = self.W[name] - float(lrs.get(name, 0.0)) * (Xs[name].T @ G)

    def acc(self, Xs, y):
        y = np.asarray(y, dtype=int)
        pred = self.logits(Xs).argmax(axis=1)
        return float(np.mean(pred == y)) if y.size else float("nan")


def consecutive_rf_share(X0, X1, seed=SEED, n_estimators=40):
    if len(X0) < 8 or len(X1) < 8:
        return {g: 1.0 / 3.0 for g in GROUP_NAMES}, float("nan")
    vimp, auc = benchmark_feature_selection(X0, X1, seed=seed, n_estimators=n_estimators)
    _, share = modality_mass(vimp)
    return share, float(auc)


def _minibatch(rng, n, chunk):
    return rng.choice(n, size=min(chunk, n), replace=False)


def run_continuous_trainer(
    bundle,
    n_batches=10,
    eta0=0.05,
    steps_per_batch=8,
    warmup_steps=4,
    n_estimators=40,
    seed=SEED,
    equal_lr=False,
):
    """Warmup each head once, then rotate.

    Code path: ``lrs[m] = eta0 * π_m`` is a which-head budget from RF-Domain
    VIMP, not the signed step. Locked signs: large covariate intensity lowers
    ``η_m``; large concept drift raises it; both quiet freezes the head.
    """
    X = standardize_columns(bundle.X)
    y_raw = np.asarray(bundle.video_id)
    classes, y = np.unique(y_raw, return_inverse=True)
    batch = assign_temporal_batches(bundle.window_idx, n_batches=n_batches)
    n_batches = int(batch.max()) + 1 if batch.size else int(n_batches)
    probe = SeparateHeadProbe(n_classes=int(classes.size), seed=seed)
    rng = np.random.default_rng(seed)
    pi_prev = {g: 1.0 / 3.0 for g in GROUP_NAMES}
    history = []
    updates = {g: 0 for g in GROUP_NAMES}

    i_warm = np.flatnonzero(batch == 0)
    Xs_w, yw = _split_modalities(X[i_warm]), y[i_warm]
    chunk_w = max(4, i_warm.size // 2)
    warm_lr = {g: eta0 / 3.0 for g in GROUP_NAMES}
    for head in GROUP_NAMES:
        for _ in range(int(max(1, warmup_steps))):
            sl = _minibatch(rng, i_warm.size, chunk_w)
            Xs = {g: Xs_w[g][sl] for g in GROUP_NAMES}
            probe.step(Xs, yw[sl], warm_lr, active=head)
            updates[head] += 1
    history.append(
        {
            "round": 0,
            "phase": "warmup",
            "batch": 0,
            "active": list(GROUP_NAMES),
            "pi": dict(pi_prev),
            "delta_pi": {g: 0.0 for g in GROUP_NAMES},
            "lr": {g: float(warm_lr[g]) for g in GROUP_NAMES},
            "auc_shift": float("nan"),
            "acc_video_id": probe.acc(Xs_w, yw),
            "n": int(i_warm.size),
            "updates": dict(updates),
        }
    )

    for t in range(1, n_batches):
        i0, i1 = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        share, auc = consecutive_rf_share(X[i0], X[i1], seed=seed + t, n_estimators=n_estimators)
        if equal_lr:
            lrs = {g: eta0 / 3.0 for g in GROUP_NAMES}
        else:
            lrs = {g: eta0 * float(share[g]) for g in GROUP_NAMES}
        Xs_all = _split_modalities(X[i1])
        yb = y[i1]
        chunk = max(4, i1.size // 2)
        active_seq = []
        n_steps = max(1, int(steps_per_batch) // 3)
        for head in GROUP_NAMES:
            for _ in range(n_steps):
                sl = _minibatch(rng, i1.size, chunk)
                Xs = {g: Xs_all[g][sl] for g in GROUP_NAMES}
                probe.step(Xs, yb[sl], lrs, active=head)
                updates[head] += 1
                active_seq.append(head)
        acc = probe.acc(Xs_all, yb)
        delta = {g: float(share[g] - pi_prev[g]) for g in GROUP_NAMES}
        history.append(
            {
                "round": t,
                "phase": "rotate",
                "batch_prev": t - 1,
                "batch": t,
                "active": active_seq,
                "pi": {g: float(share[g]) for g in GROUP_NAMES},
                "delta_pi": delta,
                "lr": {g: float(lrs[g]) for g in GROUP_NAMES},
                "auc_shift": auc,
                "acc_video_id": acc,
                "n": int(i1.size),
                "updates": dict(updates),
            }
        )
        pi_prev = {g: float(share[g]) for g in GROUP_NAMES}

    hold = np.flatnonzero(batch == n_batches - 1)
    first = np.flatnonzero(batch == 0)
    rotate = [h for h in history if h.get("phase") == "rotate"]
    pi_src = rotate or history
    summary = {
        "n_batches": n_batches,
        "n_classes": int(classes.size),
        "equal_lr": bool(equal_lr),
        "eta0": float(eta0),
        "schedule": "warmup_each_head_then_rotate",
        "mean_pi": {g: float(np.mean([h["pi"][g] for h in pi_src])) for g in GROUP_NAMES},
        "last_batch_acc": probe.acc(_split_modalities(X[hold]), y[hold]) if hold.size else float("nan"),
        "first_batch_acc": probe.acc(_split_modalities(X[first]), y[first]) if first.size else float("nan"),
        "holdout_acc": probe.acc(_split_modalities(X[first]), y[first]) if first.size else float("nan"),
        "history": history,
    }
    return summary, probe


def plot_continuous_trainer(summary, path):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import COLORS, GRID, INK, MUTED, _save, _style

    _style()
    hist = summary["history"]
    rounds = [h["round"] for h in hist]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
    fig.suptitle(
        "Continuous trainer  ·  warmup, rotate  ·  π is which-head budget, not signed η",
        fontsize=12.5,
        fontweight="bold",
        color=INK,
    )
    ax = axes[0]
    for g in GROUP_NAMES:
        ax.plot(rounds, [h["pi"][g] for h in hist], color=COLORS[g], lw=2.1, marker="o", ms=4.5, label=g.capitalize())
    ax.axhline(1.0 / 3.0, color=MUTED, ls="--", lw=0.9)
    ax.set_xlabel("round  (B_{t-1} vs B_t)")
    ax.set_ylabel("RF share  π_m")
    ax.set_title("Covariate share each round vs last", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, color=GRID)
    ax = axes[1]
    for g in GROUP_NAMES:
        ax.plot(rounds, [h["lr"][g] for h in hist], color=COLORS[g], lw=2.1, marker="o", ms=4.5, label=g.capitalize())
    ax.set_xlabel("round")
    ax.set_ylabel("prototype budget  η0 π_m")
    ax.set_title("budget only  ·  signed η: cov↓  concept↑", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(True, color=GRID)
    fig.text(
        0.06,
        -0.02,
        "Warmup: train video, then audio, then text once.  Rotate: same order each batch, others frozen.  "
        "π from benchmark_feature_selection (RF-Domain VIMP).  "
        "mean π = video %.3f, audio %.3f, text %.3f.  first-batch acc = %.3f."
        % (
            summary["mean_pi"]["video"],
            summary["mean_pi"]["audio"],
            summary["mean_pi"]["text"],
            summary.get("first_batch_acc", summary["holdout_acc"]),
        ),
        fontsize=8.2,
        color=MUTED,
    )
    return _save(fig, path)
