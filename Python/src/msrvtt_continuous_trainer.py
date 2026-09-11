"""Prototype continuous trainer: per-head LR from FSDS shares.

Streaming loop (order of heads inside a step does not matter; all three
are updated together):

  for t = 1..K-1:
      π_t ← RF-Domain VIMP on consecutive batches (B_{t-1} vs B_t)
      SGD on B_t, video-id probe, η_m = η0 * π_{t,m}
      log π_t against π_{t-1}

This is a prototype of using attribution scores as training controls,
not a production online learner. Cosine decay is a separate board.
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
    rf_domain,
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
    """Late-sum of three linear heads. One SGD step updates all heads jointly."""

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

    def step(self, Xs, y, lrs):
        y = np.asarray(y, dtype=int)
        P = _softmax(self.logits(Xs))
        n = len(y)
        G = P.copy()
        G[np.arange(n), y] -= 1.0
        G /= max(n, 1)
        for name in GROUP_NAMES:
            self.W[name] = self.W[name] - float(lrs[name]) * (Xs[name].T @ G)

    def acc(self, Xs, y):
        y = np.asarray(y, dtype=int)
        pred = self.logits(Xs).argmax(axis=1)
        return float(np.mean(pred == y)) if y.size else float("nan")


def consecutive_rf_share(X0, X1, seed=SEED, n_estimators=40):
    if len(X0) < 8 or len(X1) < 8:
        return {g: 1.0 / 3.0 for g in GROUP_NAMES}, float("nan")
    vimp, auc = rf_domain(X0, X1, seed=seed, n_estimators=n_estimators)
    _, share = modality_mass(vimp)
    return share, float(auc)


def run_continuous_trainer(
    bundle,
    n_batches=10,
    eta0=0.05,
    steps_per_batch=8,
    n_estimators=40,
    seed=SEED,
    equal_lr=False,
):
    """Stream temporal bins; set per-head LR from consecutive-batch RF shares."""
    X = standardize_columns(bundle.X)
    y_raw = np.asarray(bundle.video_id)
    classes, y = np.unique(y_raw, return_inverse=True)
    batch = assign_temporal_batches(bundle.window_idx, n_batches=n_batches)
    n_batches = int(batch.max()) + 1 if batch.size else int(n_batches)
    probe = SeparateHeadProbe(n_classes=int(classes.size), seed=seed)
    rng = np.random.default_rng(seed)
    pi_prev = {g: 1.0 / 3.0 for g in GROUP_NAMES}
    history = []

    for t in range(1, n_batches):
        i0, i1 = np.flatnonzero(batch == t - 1), np.flatnonzero(batch == t)
        share, auc = consecutive_rf_share(X[i0], X[i1], seed=seed + t, n_estimators=n_estimators)
        if equal_lr:
            lrs = {g: eta0 / 3.0 for g in GROUP_NAMES}
        else:
            lrs = {g: eta0 * float(share[g]) for g in GROUP_NAMES}
        idx = i1
        Xs_all = _split_modalities(X[idx])
        yb = y[idx]
        chunk = max(4, idx.size // 2)
        for _ in range(steps_per_batch):
            sl = rng.choice(idx.size, size=min(chunk, idx.size), replace=False)
            Xs = {g: Xs_all[g][sl] for g in GROUP_NAMES}
            probe.step(Xs, yb[sl], lrs)
        acc = probe.acc(Xs_all, yb)
        delta = {g: float(share[g] - pi_prev[g]) for g in GROUP_NAMES}
        history.append(
            {
                "round": t,
                "batch_prev": t - 1,
                "batch": t,
                "pi": {g: float(share[g]) for g in GROUP_NAMES},
                "delta_pi": delta,
                "lr": {g: float(lrs[g]) for g in GROUP_NAMES},
                "auc_shift": auc,
                "acc_video_id": acc,
                "n": int(idx.size),
            }
        )
        pi_prev = {g: float(share[g]) for g in GROUP_NAMES}

    hold = np.flatnonzero(batch == n_batches - 1)
    first = np.flatnonzero(batch == 0)
    summary = {
        "n_batches": n_batches,
        "n_classes": int(classes.size),
        "equal_lr": bool(equal_lr),
        "eta0": float(eta0),
        "mean_pi": {g: float(np.mean([h["pi"][g] for h in history])) for g in GROUP_NAMES},
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
    fig.suptitle("Continuous trainer prototype  ·  per-head LR from RF VIMP", fontsize=13, fontweight="bold", color=INK)
    ax = axes[0]
    for g in GROUP_NAMES:
        ax.plot(rounds, [h["pi"][g] for h in hist], color=COLORS[g], lw=2.1, marker="o", ms=4.5, label=g.capitalize())
    ax.axhline(1.0 / 3.0, color=MUTED, ls="--", lw=0.9)
    ax.set_xlabel("round  (B_{t-1} vs B_t)")
    ax.set_ylabel("RF share  π_m")
    ax.set_title("Attribution each round vs last", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, color=GRID)
    ax = axes[1]
    for g in GROUP_NAMES:
        ax.plot(rounds, [h["lr"][g] for h in hist], color=COLORS[g], lw=2.1, marker="o", ms=4.5, label=g.capitalize())
    ax.set_xlabel("round")
    ax.set_ylabel("head LR  η_m")
    ax.set_title("η_m = η0 · π_m   (heads updated together)", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(True, color=GRID)
    fig.text(
        0.06,
        -0.02,
        "Prototype only.  Training order inside a step is simultaneous, not sequential.  "
        "Round-to-round Δπ is the exploration hook; cosine decay is a separate board.  "
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
