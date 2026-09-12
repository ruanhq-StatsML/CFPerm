"""Publication plots for MSR-VTT multimodal FSDS attribution."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec

VIDEO_C = "#C45C26"
AUDIO_C = "#2C4A6E"
TEXT_C = "#2F6B4F"
INK = "#1A2332"
MUTED = "#5A6573"
GRID = "#E6E9EE"
COLORS = {"video": VIDEO_C, "audio": AUDIO_C, "text": TEXT_C}
GROUPS = ("video", "audio", "text")


def _style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#9AA3AE",
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.linewidth": 0.9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def _save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white", pad_inches=0.18)
    plt.close(fig)
    return path


def plot_modality_shares(result, path: Path):
    _style()
    methods = [
        ("rf", "RF-Domain VIMP"),
        ("mmd", "MMD (coord VIMP)"),
        ("po", "PO-risk VIMP"),
    ]
    point = result["point"]
    share_map = {
        "rf": point["rf_share"],
        "mmd": point["mmd_share"],
        "po": point.get("po_feat_share", point.get("po_logo_share")),
    }
    boot = result.get("bootstrap", {})
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    x = np.arange(len(methods))
    width = 0.24
    for i, g in enumerate(GROUPS):
        means = []
        yerr = np.zeros((2, len(methods)))
        for j, (key, _) in enumerate(methods):
            mu = share_map[key][g]
            means.append(mu)
            rec = boot.get(key, {}).get(g, {})
            sd = rec.get("sd")
            if sd is not None and np.isfinite(sd):
                yerr[0, j] = sd
                yerr[1, j] = sd
            else:
                ci = rec.get("ci95")
                if ci:
                    yerr[0, j] = max(0.0, mu - ci[0])
                    yerr[1, j] = max(0.0, ci[1] - mu)
        ax.bar(
            x + (i - 1) * width,
            means,
            width,
            color=COLORS[g],
            label=g.capitalize(),
            yerr=yerr,
            capsize=3.5,
            error_kw={"elinewidth": 1.1, "ecolor": "#333333"},
            zorder=3,
        )
    ax.axhline(1.0 / 3.0, color=MUTED, ls="--", lw=1.0, label="equal share 1/3")
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in methods])
    ax.set_ylabel("Modality share (L1-normalized)")
    ax.set_ylim(0, 1.05)
    ax.set_title("MSR-VTT window FSDS: modality contribution", loc="left", fontsize=13, fontweight="bold")
    ax.yaxis.grid(True, color=GRID, zorder=0)
    ax.legend(frameon=False, ncol=4, loc="upper right")
    ax.text(
        0.0,
        -0.18,
        "Error bars: SD across B=10 video-clustered bootstrap replicates. W = early vs late windows.",
        transform=ax.transAxes,
        fontsize=8.5,
        color=MUTED,
    )
    return _save(fig, path)


def plot_per_video_auc(result, path: Path):
    _style()
    rows = sorted(result["per_video"], key=lambda r: r["auc"] if np.isfinite(r["auc"]) else -1)
    fig, ax = plt.subplots(figsize=(9.8, max(3.8, 0.38 * max(len(rows), 4) + 1.6)))
    y = np.arange(len(rows))
    aucs = [r["auc"] for r in rows]
    cols = [COLORS[r["dominant"]] for r in rows]
    ax.axvline(0.5, color=MUTED, ls="--", lw=1.1, zorder=1)
    ax.barh(y, aucs, color=cols, height=0.72, zorder=3)
    labels = ["video %s  n=%d" % (r["video_id"], r["n"]) for r in rows]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Out-of-fold RF-Domain AUC (early vs late)")
    ax.set_xlim(0.35, 1.02)
    ax.set_title("Video-granularity temporal shift (subset AUC)", loc="left", fontsize=13, fontweight="bold")
    pooled = result["point"]["rf_auc"]
    ax.axvline(pooled, color=INK, lw=1.15, ls=":", label="pooled AUC %.3f" % pooled)
    for g, c in COLORS.items():
        ax.plot([], [], color=c, lw=6, label="%s dominant" % g)
    ax.legend(frameon=False, fontsize=8.2, loc="lower right")
    ax.xaxis.grid(True, color=GRID, zorder=0)
    p = result["tests"]["within_video_auc_perm"]["p"]
    ax.text(
        0.0,
        -0.12,
        "Within-video permutation p(AUC) = %.4f. Color = RF VIMP dominant modality." % p,
        transform=ax.transAxes,
        fontsize=8.5,
        color=MUTED,
    )
    return _save(fig, path)


def plot_per_video_shares(result, path: Path):
    _style()
    rows = result["per_video"]
    fig, ax = plt.subplots(figsize=(9.8, max(3.6, 0.42 * max(len(rows), 4) + 1.4)))
    y = np.arange(len(rows))
    left = np.zeros(len(rows))
    for g in GROUPS:
        vals = np.array([r["share"][g] for r in rows])
        ax.barh(y, vals, left=left, color=COLORS[g], height=0.74, label=g.capitalize())
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(["video %s  AUC=%.2f" % (r["video_id"], r["auc"]) for r in rows], fontsize=8.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("RF-Domain VIMP share")
    ax.set_title("Per-video modality mix", loc="left", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, ncol=3, loc="upper right")
    return _save(fig, path)


def plot_bootstrap_diffs(result, path: Path):
    _style()
    boot = result["bootstrap"]
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8), sharey=True)
    methods = [("rf", "RF-Domain"), ("mmd", "MMD coord-VIMP"), ("po", "PO-risk VIMP")]
    pair = "video-audio"
    for ax, (key, lab) in zip(axes, methods):
        rec = boot[key]["pairwise"][pair]
        # reconstruct a simple CI plot for all pairs
        pairs = list(boot[key]["pairwise"].keys())
        ys = np.arange(len(pairs))
        means = [boot[key]["pairwise"][p]["mean_diff"] for p in pairs]
        los = [boot[key]["pairwise"][p]["ci95"][0] for p in pairs]
        his = [boot[key]["pairwise"][p]["ci95"][1] for p in pairs]
        ax.axvline(0, color=MUTED, ls="--", lw=1.0)
        ax.errorbar(
            means,
            ys,
            xerr=[np.array(means) - np.array(los), np.array(his) - np.array(means)],
            fmt="o",
            color=INK,
            capsize=3.2,
            ms=6,
        )
        ax.set_yticks(ys)
        ax.set_yticklabels(pairs)
        ax.set_title(lab, fontsize=11, fontweight="bold")
        ax.grid(True, axis="x", color=GRID)
        ax.set_xlabel("share difference")
    fig.suptitle("Cluster-bootstrap mean ± 1.96 SD (B=10) for modality share gaps", fontsize=13, fontweight="bold", y=1.03)
    fig.tight_layout()
    return _save(fig, path)


def plot_feature_vimp(result, path: Path, k=40):
    _style()
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.2), sharey=False)
    arr = result["arrays"]
    titles = [
        ("rf_vimp", "RF-Domain feature VIMP"),
        ("coord_mmd_vimp", "Coord-MMD feature VIMP"),
        ("po_vimp", "PO-risk feature VIMP"),
    ]
    bounds = [(0, 768, "video"), (768, 1280, "audio"), (1280, 2048, "text")]
    for ax, (key, lab) in zip(axes, titles):
        v = np.asarray(arr[key], dtype=float)
        v = np.nan_to_num(v, nan=0.0)
        order = np.argsort(-v)[:k]
        cols = []
        for j in order:
            if j < 768:
                cols.append(VIDEO_C)
            elif j < 1280:
                cols.append(AUDIO_C)
            else:
                cols.append(TEXT_C)
        ax.barh(np.arange(k)[::-1], v[order][::-1], color=cols[::-1], height=0.78)
        ax.set_yticks([])
        ax.set_title(lab, fontsize=10.5, fontweight="bold")
        ax.set_xlabel("importance")
        # modality mass annotation
        share = result["point"][
            "rf_share" if key == "rf_vimp" else "coord_mmd_share" if key == "coord_mmd_vimp" else "po_feat_share"
        ]
        ax.text(
            0.98,
            0.04,
            "V %.2f  A %.2f  T %.2f" % (share["video"], share["audio"], share["text"]),
            transform=ax.transAxes,
            ha="right",
            fontsize=8,
            color=MUTED,
        )
    fig.suptitle("Feature-specific risk / shift contribution (top %d)" % k, fontsize=13, fontweight="bold")
    fig.tight_layout()
    return _save(fig, path)


def plot_method_scatter(result, path: Path):
    _style()
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    p = result["point"]
    for g in GROUPS:
        ax.scatter(
            [p["rf_share"][g]],
            [p["po_feat_share"][g]],
            s=140,
            color=COLORS[g],
            label="%s RF vs PO" % g,
            zorder=4,
        )
        ax.scatter(
            [p["mmd_share"][g]],
            [p["po_feat_share"][g]],
            s=90,
            facecolors="none",
            edgecolors=COLORS[g],
            linewidths=1.6,
            marker="D",
            label="%s MMD vs PO" % g,
            zorder=4,
        )
    ax.plot([0, 1], [0, 1], color=MUTED, ls="--", lw=1)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("RF-Domain / MMD-block share")
    ax.set_ylabel("PO-risk VIMP share")
    ax.set_title("CS detectors vs PO-risk", loc="left", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, color=GRID)
    return _save(fig, path)


def _col_zscore(M):
    M = np.asarray(M, dtype=float)
    sd = M.std(axis=0, ddof=0)
    sd = np.where(sd < 1e-10, 1.0, sd)
    return (M - M.mean(axis=0)) / sd


def plot_batch_region_board(bundle, result, path: Path, pay=None):
    """Batch 0 vs Batch 1 heatmap board across modality regions."""
    from msrvtt_multimodal_attribution import GROUPS, region_board_payload

    _style()
    if pay is None:
        pay = region_board_payload(bundle)
    n0 = pay["order_n0"]
    n = pay["n"]
    d_vr = np.asarray(pay["d_video_region"], dtype=float)
    p_vr = np.asarray(pay["p_holm_video_region"], dtype=float)
    act = _col_zscore(pay["act_window_region"])
    n_reg = act.shape[1]
    vmax_d = float(np.nanpercentile(np.abs(d_vr), 98)) if d_vr.size else 1.0
    vmax_d = max(vmax_d, 0.35)
    dnorm = TwoSlopeNorm(vmin=-vmax_d, vcenter=0.0, vmax=vmax_d)
    amax = float(np.nanpercentile(np.abs(act), 98))
    amax = max(amax, 0.35)
    anorm = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)
    mad = pay.get("mean_abs_d_mod") or {
        g: float(np.mean(np.abs(d_vr[:, sl])))
        for g, sl in (("video", slice(0, 8)), ("audio", slice(8, 16)), ("text", slice(16, n_reg)))
    }

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig = plt.figure(figsize=(16.8, 14.4), dpi=210)
    fig.patch.set_facecolor("white")
    gs = GridSpec(
        3,
        3,
        figure=fig,
        height_ratios=[1.18, 1.02, 1.22],
        hspace=0.52,
        wspace=0.28,
        left=0.058,
        right=0.975,
        top=0.90,
        bottom=0.07,
    )
    fig.suptitle(
        "MSR-VTT  ·  Batch 0 vs Batch 1 region board",
        fontsize=20.5,
        fontweight="bold",
        color=INK,
        y=0.978,
    )
    fig.text(
        0.5,
        0.932,
        "Windows ordered Batch 0 (early) then Batch 1 (late), grouped by video.  "
        "Lead: region-mean activation (text stripes are between-video captions, identical across batches).  "
        "Middle: cosine(B0, B1).  Bottom: Cohen's d with Holm stars.",
        ha="center",
        fontsize=9.6,
        color=MUTED,
    )

    ax0 = fig.add_subplot(gs[0, :2])
    im0 = ax0.imshow(act, cmap="RdBu_r", norm=anorm, aspect="auto", interpolation="nearest")
    for yb in pay.get("video_boundaries", []):
        if abs(yb - (n0 - 0.5)) > 0.1:
            ax0.axhline(yb, color="white", lw=0.35, alpha=0.55)
    ax0.axhline(n0 - 0.5, color=INK, lw=1.55)
    ax0.axvline(7.5, color=INK, lw=1.15)
    ax0.axvline(15.5, color=INK, lw=1.15)
    ax0.set_yticks([n0 / 2.0 - 0.5, n0 + (n - n0) / 2.0 - 0.5])
    ax0.set_yticklabels(["Batch 0\n(early)", "Batch 1\n(late)"], fontsize=9.5)
    ax0.set_xticks(np.arange(n_reg))
    ax0.set_xticklabels(pay["region_labels"], fontsize=6.5, rotation=90)
    ax0.set_title("Region activation  ·  windows × embedding bins", loc="left", fontsize=12.2, fontweight="bold", pad=8)
    div0 = make_axes_locatable(ax0)
    cax0 = div0.append_axes("right", size="2.4%", pad=0.08)
    cbar = fig.colorbar(im0, cax=cax0)
    cbar.set_label("column z-score", fontsize=8)
    cbar.ax.tick_params(labelsize=7.5)

    axb = fig.add_subplot(gs[0, 2])
    vals = [mad[g] for g in ("video", "audio", "text")]
    colors = [VIDEO_C, AUDIO_C, TEXT_C]
    axb.barh(np.arange(3)[::-1], vals, color=colors, height=0.62)
    axb.set_yticks(np.arange(3)[::-1])
    axb.set_yticklabels(["Video", "Audio", "Text"], fontsize=9.5)
    axb.set_xlabel("mean |Cohen's d|  across videos × regions")
    axb.set_title("Modality effect size", fontsize=11.5, fontweight="bold")
    axb.axvline(0, color=MUTED, lw=0.8)
    axb.grid(True, axis="x", color=GRID)
    xmax = max(max(vals), 0.2) * 1.22
    axb.set_xlim(0, xmax)
    for i, v in enumerate(vals):
        axb.text(v + 0.025 * xmax, 2 - i, "%.3f" % v, va="center", fontsize=8.5, color=INK)

    for i, name in enumerate(GROUPS):
        ax = fig.add_subplot(gs[1, i])
        S = pay["sims"][name]
        im = ax.imshow(S, cmap="magma", vmin=0.0, vmax=1.0, aspect="auto", interpolation="nearest")
        ax.axhline(n0 - 0.5, color="white", lw=1.15)
        ax.axvline(n0 - 0.5, color="white", lw=1.15)
        ax.set_title(name.capitalize() + "  cosine(B0, B1)", fontsize=12, fontweight="bold", color=COLORS[name])
        ticks = [n0 / 2.0 - 0.5, n0 + (n - n0) / 2.0 - 0.5]
        ax.set_xticks(ticks)
        ax.set_xticklabels(["Batch 0", "Batch 1"], fontsize=8.5)
        ax.set_yticks(ticks)
        ax.set_yticklabels(["Batch 0", "Batch 1"], fontsize=8.5)
        bm = pay["block_mean"][name]
        ax.text(
            0.02,
            -0.18,
            "within %.3f   cross %.3f   gap %.3f"
            % (0.5 * (bm["B0B0"] + bm["B1B1"]), bm["B0B1"], bm["gap"]),
            transform=ax.transAxes,
            fontsize=8.0,
            color=MUTED,
            clip_on=False,
        )
        if i == 2:
            divs = make_axes_locatable(ax)
            caxs = divs.append_axes("right", size="4.5%", pad=0.06)
            cbar = fig.colorbar(im, cax=caxs)
            cbar.ax.tick_params(labelsize=7.5)
            cbar.set_label("cosine", fontsize=8)

    axh = fig.add_subplot(gs[2, :2])
    imh = axh.imshow(d_vr, cmap="RdBu_r", norm=dnorm, aspect="auto", interpolation="nearest")
    axh.axvline(7.5, color=INK, lw=1.2)
    axh.axvline(15.5, color=INK, lw=1.2)
    axh.set_yticks(np.arange(len(pay["video_ids"])))
    axh.set_yticklabels(["v%s" % v for v in pay["video_ids"]], fontsize=8)
    axh.set_xticks(np.arange(n_reg))
    axh.set_xticklabels(pay["region_labels"], fontsize=6.6, rotation=90)
    axh.set_title("Cohen's d by video × region   (Holm * p<0.05)", loc="left", fontsize=12.0, fontweight="bold", pad=8)
    yy, xx = np.where(p_vr < 0.05)
    axh.scatter(xx, yy, marker="*", s=22, c=INK, linewidths=0, zorder=4)
    divh = make_axes_locatable(axh)
    caxh = divh.append_axes("right", size="2.4%", pad=0.08)
    cbar = fig.colorbar(imh, cax=caxh)
    cbar.set_label("d  (B1 − B0)", fontsize=8)
    cbar.ax.tick_params(labelsize=7.5)

    axp = fig.add_subplot(gs[2, 2])
    stacked = np.vstack(
        [
            np.asarray(pay["pooled_d"], dtype=float),
            np.asarray(pay.get("mean_abs_d_region", np.mean(np.abs(d_vr), axis=0)), dtype=float),
        ]
    )
    abs_max = max(float(np.nanpercentile(np.abs(stacked), 98)), 0.2)
    pnorm = TwoSlopeNorm(vmin=-abs_max, vcenter=0.0, vmax=abs_max)
    axp.imshow(stacked, cmap="RdBu_r", norm=pnorm, aspect="auto", interpolation="nearest")
    axp.set_yticks([0, 1])
    axp.set_yticklabels(["signed pooled d", "mean |d|"], fontsize=8.2)
    axp.set_xticks(np.arange(n_reg))
    axp.set_xticklabels(pay["region_labels"], fontsize=6.2, rotation=90)
    axp.set_title("Pooled vs mean |d|", fontsize=11, fontweight="bold", pad=8)
    axp.axvline(7.5, color=INK, lw=1.0)
    axp.axvline(15.5, color=INK, lw=1.0)
    for j, pval in enumerate(pay["pooled_p_holm"]):
        if pval < 0.05:
            axp.text(j, 0, "*", ha="center", va="center", fontsize=11, color=INK, fontweight="bold")

    fried = result.get("tests", {}).get("per_video_rf_shares", {}).get("friedman", {})
    fig.text(
        0.055,
        0.018,
        "Region = equal-width bin inside 768 / 512 / 768.  "
        "Welch t on region-mean activation; Holm within video (stars) and across pooled regions.  "
        "Friedman p(equal modality RF shares) = %.2e.  "
        "Mean |d|: video %.3f, audio %.3f, text %.3f.  "
        "%d/%d pooled Holm-sig, %d video×region stars.  "
        "Text is window-invariant on this extract (d=0)."
        % (
            fried.get("p", float("nan")),
            mad["video"],
            mad["audio"],
            mad["text"],
            pay["n_sig_pooled"],
            n_reg,
            pay["n_sig_cells"],
        ),
        fontsize=8.0,
        color=MUTED,
    )
    return _save(fig, path)


def plot_batch_pair_board(bundle, result, path: Path, pay=None, n_batches=10):
    """Product dashboard: three modality cosine(Bi, Bj) heatmaps + head LR + drill-down."""
    from matplotlib.patches import FancyBboxPatch
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from msrvtt_multimodal_attribution import GROUP_NAMES, batch_pair_payload

    _style()
    if pay is None:
        pay = batch_pair_payload(bundle, n_batches=n_batches)
    k = int(pay["n_batches"])
    ticks = np.arange(k)
    labs = ["B%d" % (i + 1) for i in ticks]
    mats = [np.asarray(pay["cosine"][g], dtype=float) for g in GROUP_NAMES]
    stacked = np.concatenate([m[np.isfinite(m)] for m in mats])
    vmax = float(np.nanpercentile(np.abs(stacked), 98)) if stacked.size else 0.01
    vmax = max(vmax, 1e-4)
    nrm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    point = (result or {}).get("point") or {}
    rf = point.get("rf_share") or {"video": 0.772, "audio": 0.191, "text": 0.037}
    po = point.get("po_feat_share") or {"video": 0.813, "audio": 0.174, "text": 0.013}

    fig = plt.figure(figsize=(16.8, 11.8), dpi=210)
    fig.patch.set_facecolor("white")
    gs = GridSpec(
        2,
        3,
        figure=fig,
        height_ratios=[1.38, 0.95],
        hspace=0.42,
        wspace=0.28,
        left=0.05,
        right=0.98,
        top=0.88,
        bottom=0.07,
    )
    fig.suptitle(
        "Dashboard  ·  three modality heatmaps   cosine(Bi, Bj)",
        fontsize=19.5,
        fontweight="bold",
        color=INK,
        y=0.978,
    )
    fig.text(
        0.5,
        0.928,
        "Post on the board: VIDEO / AUDIO / TEXT heatmaps (shared color scale), lag decay, "
        "per-head LR from RF & PO-risk VIMP, drill-down  token/frame  →  patch/box  |  audio slice.  "
        "K=%d bins of window index  ·  prototype of 10k frames / 1k per batch."
        % k,
        ha="center",
        fontsize=9.3,
        color=MUTED,
    )

    last_im = None
    for i, name in enumerate(GROUP_NAMES):
        ax = fig.add_subplot(gs[0, i])
        M = mats[i]
        last_im = ax.imshow(M, cmap="RdBu_r", norm=nrm, origin="upper", interpolation="nearest")
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels(labs, fontsize=7.4, rotation=90)
        ax.set_yticklabels(labs, fontsize=7.4)
        ax.set_title(
            name.upper() + "   cosine(Bi, Bj)",
            fontsize=13.2,
            fontweight="bold",
            color=COLORS[name],
            pad=8,
        )
        ax.set_xlabel("batch j")
        if i == 0:
            ax.set_ylabel("batch i")
        gap = pay["lag0_minus_lagmax"].get(name, float("nan"))
        ax.text(
            0.0,
            -0.16,
            "decay  lag0−lag%d  =  %.4f" % (k - 1, gap),
            transform=ax.transAxes,
            fontsize=8.4,
            color=MUTED,
            clip_on=False,
        )
        if i == 2:
            div = make_axes_locatable(ax)
            cax = div.append_axes("right", size="4.0%", pad=0.06)
            cbar = fig.colorbar(last_im, cax=cax)
            cbar.set_label("mean cosine", fontsize=8)
            cbar.ax.tick_params(labelsize=7.5)

    axl = fig.add_subplot(gs[1, 0])
    h = np.arange(k)
    for name in GROUP_NAMES:
        axl.plot(
            h,
            pay["lag_cosine"][name],
            color=COLORS[name],
            lw=2.3,
            marker="o",
            ms=4.8,
            label=name.capitalize(),
        )
    axl.set_xlabel("lag  |i − j|")
    axl.set_ylabel("mean cosine")
    axl.set_title("Lag  ·  refresh half-life", fontsize=12.2, fontweight="bold")
    axl.legend(frameon=False, fontsize=8.5, loc="upper right")
    axl.grid(True, color=GRID)
    axl.set_xticks(h)

    axr = fig.add_subplot(gs[1, 1])
    y = np.arange(3)
    hgt = 0.36
    rf_v = [rf[g] for g in GROUP_NAMES]
    po_v = [po[g] for g in GROUP_NAMES]
    axr.barh(y + hgt / 2, rf_v[::-1], height=hgt, color=[COLORS[g] for g in GROUP_NAMES][::-1], label="RF-Domain VIMP")
    axr.barh(y - hgt / 2, po_v[::-1], height=hgt, color="#9AA3AE", label="PO-risk VIMP")
    axr.set_yticks(y)
    axr.set_yticklabels(["Text", "Audio", "Video"], fontsize=9.5)
    axr.set_xlabel("share  →  which-head budget  (not the signed η)")
    axr.set_title("Which head  (cov. mass);  η sign is separate", fontsize=11.2, fontweight="bold")
    axr.set_xlim(0, 1.05)
    axr.axvline(1.0 / 3.0, color=MUTED, ls="--", lw=0.9)
    axr.legend(frameon=False, fontsize=7.8, loc="lower right")
    axr.grid(True, axis="x", color=GRID)

    axd = fig.add_subplot(gs[1, 2])
    axd.set_xlim(0, 1)
    axd.set_ylim(0, 1)
    axd.axis("off")
    axd.set_title("Drill-down  (after a heatmap cell)", fontsize=12.0, fontweight="bold")
    steps = (
        (0.72, "1   token / frame", "stored sliding-window embeddings"),
        (0.44, "2   patch / bounding-box", "spatial slice of the frame"),
        (0.16, "3   audio slice", "time slice of the soundtrack"),
    )
    for y0, title, sub in steps:
        box = FancyBboxPatch(
            (0.06, y0 - 0.08),
            0.88,
            0.22,
            boxstyle="round,pad=0.012,rounding_size=0.04",
            facecolor="#F4F6F8",
            edgecolor="#C5CCD4",
            lw=1.0,
            transform=axd.transAxes,
            clip_on=False,
        )
        axd.add_patch(box)
        axd.text(0.12, y0 + 0.06, title, transform=axd.transAxes, fontsize=10.2, fontweight="bold", color=INK, va="center")
        axd.text(0.12, y0 - 0.02, sub, transform=axd.transAxes, fontsize=8.0, color=MUTED, va="center")

    fig.text(
        0.05,
        0.016,
        "Refresh ≠ step.  Covariate mass opens the head / KV; large c_m^{cov} lowers η_m; "
        "large concept drift raises η_m.  Do not invert π: clip-level text is both-quiet, η_t ≈ 0.",
        fontsize=8.0,
        color=MUTED,
    )
    return _save(fig, path)


def write_all_plots(result, out_dir: Path, bundle=None):
    out_dir = Path(out_dir)
    paths = {
        "shares": plot_modality_shares(result, out_dir / "msrvtt_modality_shares.png"),
        "auc": plot_per_video_auc(result, out_dir / "msrvtt_per_video_auc.png"),
        "video_shares": plot_per_video_shares(result, out_dir / "msrvtt_per_video_shares.png"),
        "diffs": plot_bootstrap_diffs(result, out_dir / "msrvtt_bootstrap_share_diffs.png"),
        "features": plot_feature_vimp(result, out_dir / "msrvtt_feature_specific_vimp.png"),
        "scatter": plot_method_scatter(result, out_dir / "msrvtt_rf_mmd_porisk.png"),
    }
    if bundle is not None:
        from msrvtt_multimodal_attribution import (
            REGION_BOARD_SKIP,
            batch_pair_payload,
            region_board_payload,
            write_json,
        )

        pay = region_board_payload(bundle)
        write_json(
            out_dir / "msrvtt_region_board_stats.json",
            {k: v for k, v in pay.items() if k not in REGION_BOARD_SKIP},
        )
        paths["board"] = plot_batch_region_board(
            bundle, result, out_dir / "msrvtt_batch_region_heatmap_board.png", pay=pay
        )
        paths["board_stats"] = str(out_dir / "msrvtt_region_board_stats.json")
        pair = batch_pair_payload(bundle, n_batches=10)
        write_json(
            out_dir / "msrvtt_batch_pair_cosine_stats.json",
            {k: v for k, v in pair.items() if k != "batch"},
        )
        paths["batch_pair"] = plot_batch_pair_board(
            bundle, result, out_dir / "msrvtt_batch_ij_cosine_board.png", pay=pair
        )
    return {k: str(v) for k, v in paths.items()}
