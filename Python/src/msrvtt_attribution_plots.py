"""Publication plots for MSR-VTT multimodal FSDS attribution."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

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
        ("mmd", "MMD-LOCO"),
        ("po", "PO-risk LOGO"),
    ]
    point = result["point"]
    share_map = {
        "rf": point["rf_share"],
        "mmd": point["mmd_share"],
        "po": point["po_logo_share"],
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
            ci = boot.get(key, {}).get(g, {}).get("ci95")
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
        "Error bars: video-clustered bootstrap 95% CI. W = early vs late windows.",
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
    methods = [("rf", "RF-Domain"), ("mmd", "MMD-LOCO"), ("po", "PO-risk LOGO")]
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
    fig.suptitle("Cluster-bootstrap 95% CI for modality share gaps", fontsize=13, fontweight="bold", y=1.03)
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
            [p["po_logo_share"][g]],
            s=140,
            color=COLORS[g],
            label="%s RF vs PO" % g,
            zorder=4,
        )
        ax.scatter(
            [p["mmd_share"][g]],
            [p["po_logo_share"][g]],
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
    ax.set_xlabel("RF-Domain / MMD-LOCO share")
    ax.set_ylabel("PO-risk LOGO share")
    ax.set_title("CS detectors vs PO-risk", loc="left", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, color=GRID)
    return _save(fig, path)


def write_all_plots(result, out_dir: Path):
    out_dir = Path(out_dir)
    paths = {
        "shares": plot_modality_shares(result, out_dir / "msrvtt_modality_shares.png"),
        "auc": plot_per_video_auc(result, out_dir / "msrvtt_per_video_auc.png"),
        "video_shares": plot_per_video_shares(result, out_dir / "msrvtt_per_video_shares.png"),
        "diffs": plot_bootstrap_diffs(result, out_dir / "msrvtt_bootstrap_share_diffs.png"),
        "features": plot_feature_vimp(result, out_dir / "msrvtt_feature_specific_vimp.png"),
        "scatter": plot_method_scatter(result, out_dir / "msrvtt_rf_mmd_porisk.png"),
    }
    return {k: str(v) for k, v in paths.items()}
