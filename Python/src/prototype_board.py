"""Inspectable prototype board: tensor shapes plus one hop of class means.

The FSDS layout is fixed in ``msrvtt_multimodal_attribution.GROUPS``:

    s  : (n, 2049) = [768 video | 512 audio | 768 text | 1 label]
    X  : (n, 2048)
    μ  : (C, d_m)  class mean of one block, the next-batch NCM prototype

This module dumps those shapes and draws one covariate hop so the stale
prototype, the SDC-compensated prototype, and the true new mean can be
looked at in the same plane. Not a new method.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from gpm_fsds import representation_bases
from instance_discrimination import LinearProjector, MemoryBank, l2_normalize
from msrvtt_continuous_trainer import _split_modalities
from msrvtt_multimodal_attribution import (
    D_AUDIO,
    D_TEXT,
    D_VIDEO,
    GROUP_NAMES,
    GROUPS,
    P_X,
    SEED,
)
from prototype_drift import (
    class_prototypes,
    location_shift_field,
    prototype_cosine,
    sdc_compensate,
)

PROJ_DIM = 32
GPM_K = 8


def layout_spec():
    """Where the shapes are declared. Slices are half-open."""
    return {
        "declared_in": "Python/src/msrvtt_multimodal_attribution.py",
        "constants": {"D_VIDEO": D_VIDEO, "D_AUDIO": D_AUDIO, "D_TEXT": D_TEXT, "P_X": P_X, "P_s": P_X + 1},
        "groups": {g: {"start": int(sl.start), "stop": int(sl.stop), "dim": int(sl.stop - sl.start)} for g, sl in GROUPS.items()},
        "s": "[n, 2049] = concat(video 768, audio 512, text 768, label 1)",
        "msrvtt_extract": {"s": [320, 2049], "n_videos": 16, "n_windows": 20, "on_disk": "data/msrvtt/features_window/meta.json (npy not in this checkout)"},
    }


def _pca_fit(X):
    X = np.asarray(X, dtype=float)
    mu = X.mean(axis=0)
    Z = X - mu
    _, _, Vt = np.linalg.svd(Z, full_matrices=False)
    W = Vt[:2].T
    return mu, W


def _pca_apply(X, mu, W):
    return (np.asarray(X, dtype=float) - mu) @ W


def prototype_hop(stream, t=1):
    """One consecutive pair: μ_{t-1}, SDC(μ_{t-1}), μ_t. Arrays keep native d_m."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=int)
    batch = np.asarray(stream.batch, dtype=int)
    n_classes = int(y.max()) + 1
    ip = np.flatnonzero(batch == int(t) - 1)
    ic = np.flatnonzero(batch == int(t))
    out = {
        "t": int(t),
        "n_prev": int(ip.size),
        "n_curr": int(ic.size),
        "n_classes": n_classes,
        "y_prev": y[ip],
        "y_curr": y[ic],
        "blocks": {},
    }
    Xsp, Xsc = _split_modalities(X[ip]), _split_modalities(X[ic])
    for g in GROUP_NAMES:
        mu_old, _ = class_prototypes(Xsp[g], y[ip], n_classes=n_classes)
        mu_true, _ = class_prototypes(Xsc[g], y[ic], n_classes=n_classes)
        field = location_shift_field(Xsp[g], Xsc[g])
        mu_sdc = sdc_compensate(mu_old, Xsc[g] - field, Xsc[g])
        out["blocks"][g] = {
            "X_prev": Xsp[g],
            "X_curr": Xsc[g],
            "mu_old": mu_old,
            "mu_sdc": mu_sdc,
            "mu_true": mu_true,
            "field": field,
            "cos_stale": prototype_cosine(mu_old, mu_true),
            "cos_sdc": prototype_cosine(mu_sdc, mu_true),
            "per_class_cos_stale": _per_class_cosine(mu_old, mu_true),
            "per_class_cos_sdc": _per_class_cosine(mu_sdc, mu_true),
        }
    return out


def _per_class_cosine(a, b):
    an, _ = l2_normalize(a)
    bn, _ = l2_normalize(b)
    return (an * bn).sum(axis=1)


def instdisc_shapes(Xm, extra_n, dim=PROJ_DIM, seed=SEED):
    n, d = np.asarray(Xm).shape
    proj = LinearProjector(d, dim, seed=seed)
    bank = MemoryBank(n, dim, seed=seed + 1)
    f, _, _ = proj.encode(Xm)
    extra = np.zeros((int(extra_n), dim)) if extra_n else np.zeros((0, dim))
    keys = np.vstack([bank.V, extra]) if extra.size else bank.V
    logits = f @ keys.T
    return {
        "X": [n, d],
        "W": list(proj.W.shape),
        "f": list(f.shape),
        "V": list(bank.V.shape),
        "extra_keys": [int(extra_n), dim],
        "logits": list(logits.shape),
    }


def inspect_shapes(stream, t=1, dim=PROJ_DIM, seed=SEED):
    """Live shapes for the stream plus one InstDisc / SDC / GPM hop."""
    X = np.asarray(stream.X, dtype=float)
    y = np.asarray(stream.y, dtype=int)
    batch = np.asarray(stream.batch, dtype=int)
    hop = prototype_hop(stream, t=t)
    n_classes = hop["n_classes"]
    spec = layout_spec()
    blocks = {}
    for g in GROUP_NAMES:
        b = hop["blocks"][g]
        d = b["X_curr"].shape[1]
        M, explained, k = representation_bases(b["X_prev"], thresh=0.90, max_k=GPM_K)
        blocks[g] = {
            "X_prev": list(b["X_prev"].shape),
            "X_curr": list(b["X_curr"].shape),
            "mu_old": list(b["mu_old"].shape),
            "mu_sdc": list(b["mu_sdc"].shape),
            "mu_true": list(b["mu_true"].shape),
            "field": list(b["field"].shape),
            "cos_stale": float(b["cos_stale"]),
            "cos_sdc": float(b["cos_sdc"]),
            "instdisc": instdisc_shapes(b["X_curr"], hop["n_prev"], dim=dim, seed=seed),
            "gpm_M": list(M.shape),
            "gpm_k": int(k),
            "gpm_explained": float(explained),
            "dW": [d, n_classes],
        }
    return {
        "layout": spec,
        "stream": {
            "X": list(X.shape),
            "y": list(y.shape),
            "batch": list(batch.shape),
            "n_batches": int(batch.max()) + 1,
            "n_classes": int(n_classes),
            "n_per": int(hop["n_curr"]),
            "meta": dict(getattr(stream, "meta", {}) or {}),
        },
        "hop": {"t": hop["t"], "n_prev": hop["n_prev"], "n_curr": hop["n_curr"]},
        "blocks": blocks,
    }


def arrays_for_npz(hop, head="video"):
    b = hop["blocks"][head]
    return {
        "X_prev": b["X_prev"],
        "X_curr": b["X_curr"],
        "y_prev": hop["y_prev"],
        "y_curr": hop["y_curr"],
        "mu_old": b["mu_old"],
        "mu_sdc": b["mu_sdc"],
        "mu_true": b["mu_true"],
        "field": b["field"],
        "per_class_cos_stale": b["per_class_cos_stale"],
        "per_class_cos_sdc": b["per_class_cos_sdc"],
    }


def plot_prototype_board(hop, shapes, path, head="video"):
    from msrvtt_attribution_plots import COLORS, GRID, INK, MUTED, _save, _style
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    _style()
    b = hop["blocks"][head]
    pooled = np.vstack([b["X_prev"], b["X_curr"]])
    pca_mu, pca_w = _pca_fit(pooled)
    z_curr = _pca_apply(b["X_curr"], pca_mu, pca_w)
    z_old = _pca_apply(b["mu_old"], pca_mu, pca_w)
    z_sdc = _pca_apply(b["mu_sdc"], pca_mu, pca_w)
    z_true = _pca_apply(b["mu_true"], pca_mu, pca_w)
    y = hop["y_curr"]
    n_classes = hop["n_classes"]
    cmap = plt.get_cmap("tab10")

    fig = plt.figure(figsize=(12.2, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.35], hspace=0.42, wspace=0.28)
    ax0 = fig.add_subplot(gs[0, :])
    ax0.set_axis_off()
    rows = _shape_rows(shapes, hop, head)
    table = ax0.table(
        cellText=rows,
        colLabels=["tensor", "shape", "where / what"],
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.2)
    table.scale(1.0, 1.35)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#E6E9EE")
        if r == 0:
            cell.set_facecolor("#F4F6F8")
            cell.set_text_props(fontweight="bold", color=INK)
        else:
            cell.set_text_props(color=INK if c else INK)
    ax0.set_title("Data layout (declared in GROUPS) and one hop of tensors", loc="left", fontsize=12, fontweight="bold")

    ax1 = fig.add_subplot(gs[1, 0])
    ax1.scatter(z_curr[:, 0], z_curr[:, 1], c=[cmap(int(c) % 10) for c in y], s=14, alpha=0.35, linewidths=0)
    for c in range(n_classes):
        col = cmap(c % 10)
        ax1.scatter(*z_old[c], marker="s", s=70, facecolors="none", edgecolors=col, linewidths=1.4, zorder=3)
        ax1.scatter(*z_sdc[c], marker="^", s=78, c=[col], zorder=4)
        ax1.scatter(*z_true[c], marker="o", s=78, c=[col], zorder=4)
        ax1.annotate("", xy=z_sdc[c], xytext=z_old[c], arrowprops=dict(arrowstyle="->", color=col, lw=1.1))
    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")
    ax1.grid(True, color=GRID)
    ax1.set_title("%s prototypes, hop t=%d→%d" % (head, hop["t"] - 1, hop["t"]), loc="left", fontsize=12, fontweight="bold")
    ax1.legend(
        handles=[
            Line2D([0], [0], marker="s", color="none", markeredgecolor=INK, markerfacecolor="none", label="stale μ_{t−1}"),
            Line2D([0], [0], marker="^", color="none", markerfacecolor=INK, label="SDC μ"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=INK, label="true μ_t"),
        ],
        frameon=False,
        fontsize=8,
        loc="best",
    )

    ax2 = fig.add_subplot(gs[1, 1])
    xs = np.arange(n_classes)
    w = 0.36
    ax2.bar(xs - w / 2, b["per_class_cos_stale"], w, color="#9AA3AE", label="stale")
    ax2.bar(xs + w / 2, b["per_class_cos_sdc"], w, color=COLORS[head], label="SDC")
    ax2.set_xticks(xs)
    ax2.set_xticklabels(["c%d" % c for c in xs])
    ax2.set_ylim(0.0, 1.05)
    ax2.set_ylabel("cos(μ, true μ_t)")
    ax2.grid(True, color=GRID, axis="y")
    ax2.set_title("Per-class prototype cosine", loc="left", fontsize=12, fontweight="bold")
    ax2.legend(frameon=False, fontsize=8)
    fig.suptitle("Prototype board", fontsize=13.2, fontweight="bold", color=INK, x=0.04, ha="left")
    fig.text(
        0.04,
        0.01,
        "Squares are last-batch class means. Arrows are SDC. Circles are the arriving batch's true means. "
        "s is (n, 2049); prototypes are (C, d_m).",
        fontsize=8.2,
        color=MUTED,
    )
    return _save(fig, path)


def _shape_rows(shapes, hop, head):
    st = shapes["stream"]
    bl = shapes["blocks"][head]
    idc = bl["instdisc"]
    return [
        ["s / X", "%s / %s" % (shapes["layout"]["s"].split(" = ")[0], st["X"]), "GROUPS in msrvtt_multimodal_attribution.py"],
        ["X^(%s) this hop" % head, str(bl["X_curr"]), "one FSDS block, n=%d windows" % hop["n_curr"]],
        ["μ_old, μ_sdc, μ_true", str(bl["mu_old"]), "C=%d class means, next-batch NCM" % hop["n_classes"]],
        ["field Δ", str(bl["field"]), "μ_t − μ_{t−1} in R^{d_m}"],
        ["projector W, f, bank V", "%s, %s, %s" % (idc["W"], idc["f"], idc["V"]), "Wu InstDisc, k=%d" % idc["f"][1]],
        ["logits vs V ⊕ extra", str(idc["logits"]), "current slots + previous-batch keys"],
        ["GPM M, ΔW", "%s, %s" % (bl["gpm_M"], bl["dW"]), "CGS of X_{t−1}, linear-head grad"],
    ]


def write_prototype_npz(hop, path, head="video"):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays_for_npz(hop, head=head))
    return path
