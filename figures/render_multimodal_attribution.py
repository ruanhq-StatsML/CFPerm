#!/usr/bin/env python3
"""Render a polished 3-layer multimodal attribution procedure figure."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ORANGE = "#E07A3D"
ORANGE_DK = "#C45C26"
PEACH = "#F4B183"
NAVY = "#2C4A6E"
NAVY_DK = "#1B334D"
WHITE = "#FFFFFF"
INK = "#1A2332"
MUTED = "#5A6573"

L1_BG = "#FFF7F0"
L2_BG = "#F2F6FB"
L3_BG = "#F1F8F4"
HEADER_L1 = "#C45C26"
HEADER_L2 = "#2C4A6E"
HEADER_L3 = "#2F6B4F"


def rbox(ax, x, y, w, h, fc, ec, lw=1.35, rs=0.11, z=3, ls="-"):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.01,rounding_size={rs}",
        facecolor="none" if fc == "none" else fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=z,
        linestyle=ls,
        mutation_aspect=0.85,
    )
    ax.add_patch(patch)
    return patch


def label(ax, x, y, text, **kw):
    opts = dict(ha="center", va="center", zorder=6, color=WHITE, fontweight="bold")
    opts.update(kw)
    ax.text(x, y, text, **opts)


def box_label(ax, x, y, w, h, text, **kw):
    label(ax, x + w / 2, y + h / 2, text, **kw)


def arrow(ax, p0, p1, color=NAVY, lw=1.65, rad=0.0, ms=13):
    ax.add_patch(
        FancyArrowPatch(
            p0,
            p1,
            arrowstyle="-|>,head_width=0.26,head_length=0.36",
            mutation_scale=ms,
            linewidth=lw,
            color=color,
            connectionstyle=f"arc3,rad={rad}",
            zorder=4,
            shrinkA=0.4,
            shrinkB=0.4,
        )
    )


def header_bar(ax, x, y, w, h, color, badge, title):
    rbox(ax, x, y, w, h, color, color, lw=0, rs=0.08, z=3)
    rbox(ax, x + 0.16, y + 0.09, 1.22, h - 0.18, WHITE, WHITE, lw=0, rs=0.08, z=4)
    label(ax, x + 0.77, y + h / 2, badge, color=color, fontsize=8.7)
    label(ax, x + w / 2 + 0.32, y + h / 2, title, color=WHITE, fontsize=13.0)


def vimp_card(ax, x, y, w, h, title, lines, title_c):
    rbox(ax, x, y, w, h, WHITE, ORANGE_DK, lw=1.15, rs=0.10)
    rbox(ax, x + 0.08, y + h - 0.34, w - 0.16, 0.26, title_c, title_c, lw=0, rs=0.07, z=5)
    label(ax, x + w / 2, y + h - 0.21, title, fontsize=8.6)
    ax.text(
        x + 0.28,
        y + h - 0.48,
        lines,
        ha="left",
        va="top",
        fontsize=8.7,
        color=NAVY_DK,
        family="DejaVu Sans Mono",
        zorder=6,
        linespacing=1.38,
    )


def main():
    out = Path(__file__).resolve().parent / "multimodal_attribution_procedure.png"
    W, H = 17.5, 13.15

    fig, ax = plt.subplots(figsize=(W, H), dpi=230)
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor(WHITE)
    ax.set_facecolor(WHITE)

    ax.text(
        W / 2,
        12.78,
        "Hierarchical Multimodal Attribution Procedure",
        ha="center",
        va="center",
        fontsize=21.5,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        W / 2,
        12.36,
        "Adapting Feature Selection for Distribution Shift (FSDS) to user-preference / mixture-shift attribution",
        ha="center",
        va="center",
        fontsize=11.3,
        color=MUTED,
        fontstyle="italic",
    )

    # Input
    rbox(ax, 0.32, 11.48, 16.86, 0.62, "#EEF2F7", "#C5D0DC", lw=1.15, rs=0.12, z=1)
    label(ax, 0.95, 11.79, "INPUT", color=NAVY, fontsize=8.6)
    rbox(ax, 1.55, 11.59, 3.55, 0.40, ORANGE, ORANGE_DK, lw=0.85, rs=0.09)
    box_label(ax, 1.55, 11.59, 3.55, 0.40, "Sentence / Image Batch 0", fontsize=9.6)
    rbox(ax, 5.25, 11.59, 3.55, 0.40, ORANGE, ORANGE_DK, lw=0.85, rs=0.09)
    box_label(ax, 5.25, 11.59, 3.55, 0.40, "Sentence / Image Batch 1", fontsize=9.6)
    label(ax, 9.05, 11.79, "+", color=NAVY, fontsize=14)
    rbox(ax, 9.30, 11.59, 4.55, 0.40, PEACH, ORANGE_DK, lw=0.85, rs=0.09)
    box_label(ax, 9.30, 11.59, 4.55, 0.40, "Sentiment  or  Image Category", fontsize=9.6, color=NAVY_DK)
    label(
        ax,
        15.60,
        11.79,
        "preference / mixture shift",
        color=MUTED,
        fontsize=8.3,
        fontstyle="italic",
        fontweight="normal",
    )

    # ── LAYER 1 ──────────────────────────────────────────────────
    x, y, w, h = 0.28, 7.28, 16.94, 3.98
    rbox(ax, x, y, w, h, L1_BG, "#E8D5C6", lw=1.2, rs=0.16, z=1)
    header_bar(
        ax,
        x + 0.08,
        y + h - 0.50,
        w - 0.16,
        0.42,
        HEADER_L1,
        "LAYER 1",
        "Modality Attribution   —   which modality drives the shift?",
    )
    ax.text(
        x + 0.28,
        y + h - 0.72,
        "Score every modality with FSDS-style detectors and rank features. Localization instead of a unique CS / CD decomposition.",
        fontsize=9.5,
        color=MUTED,
        va="center",
        zorder=5,
    )

    ax.text(0.72, 10.12, "DETECTORS", fontsize=8.15, fontweight="bold", color=HEADER_L1, zorder=5)
    rbox(ax, 0.58, 8.95, 4.55, 0.95, ORANGE, ORANGE_DK, lw=1.2, rs=0.12)
    box_label(ax, 0.58, 8.95, 4.55, 0.95, "Covariate Shift    P(X)\nMMD  ·  RF Domain Classifier", fontsize=10.6)
    rbox(ax, 0.58, 7.62, 4.55, 0.95, ORANGE, ORANGE_DK, lw=1.2, rs=0.12)
    box_label(
        ax,
        0.58,
        7.62,
        4.55,
        0.95,
        "Concept Drift    P(Y | X)\nMeta-Learner as Distance Estimator",
        fontsize=10.6,
    )

    # detectors merge into one trunk
    ax.plot([5.18, 5.55], [9.42, 8.70], color=NAVY, lw=1.45, zorder=4)
    ax.plot([5.18, 5.55], [8.10, 8.70], color=NAVY, lw=1.45, zorder=4)
    arrow(ax, (5.55, 8.70), (6.22, 8.70), lw=1.7, ms=12)

    ax.text(6.40, 10.12, "MODALITIES  (score all three)", fontsize=8.15, fontweight="bold", color=HEADER_L1, zorder=5)
    mw, mh, my = 2.08, 1.88, 7.55
    for i, lab in enumerate(["Image", "Audio", "Text\nEmbedding"]):
        rbox(ax, 6.32 + i * 2.28, my, mw, mh, ORANGE, ORANGE_DK, lw=1.2, rs=0.14)
        box_label(ax, 6.32 + i * 2.28, my, mw, mh, lab, fontsize=13.2)

    ax.text(14.95, 10.12, "RANKED VIMP", fontsize=8.15, fontweight="bold", color=HEADER_L1, zorder=5, ha="center")
    vimp_card(
        ax,
        13.20,
        8.88,
        3.70,
        1.10,
        "Covariate Shift",
        "Image  f1–f10\nAudio  f1–f5\nText   f1–f3",
        HEADER_L1,
    )
    vimp_card(
        ax,
        13.20,
        7.52,
        3.70,
        1.10,
        "Concept Drift",
        "Image  f1–f10\nAudio  f1–f5\nText   f1–f3",
        HEADER_L1,
    )
    arrow(ax, (13.00, 8.70), (13.18, 8.70), lw=1.55, ms=11)

    # ── LAYER 2 ──────────────────────────────────────────────────
    x, y, w, h = 0.28, 3.62, 16.94, 3.32
    rbox(ax, x, y, w, h, L2_BG, "#C9D4E4", lw=1.2, rs=0.16, z=1)
    header_bar(
        ax,
        x + 0.08,
        y + h - 0.50,
        w - 0.16,
        0.42,
        HEADER_L2,
        "LAYER 2",
        "Instance Localization   —   which samples carry the shift?",
    )
    ax.text(
        x + 0.28,
        y + h - 0.72,
        "Condition on the selected modality. Score samples in Batch 0 vs. Batch 1 and keep those that drive the discrepancy.",
        fontsize=9.5,
        color=MUTED,
        va="center",
        zorder=5,
    )

    ax.text(0.72, 5.78, "SHIFT SCORERS", fontsize=8.15, fontweight="bold", color=HEADER_L2, zorder=5)
    rbox(ax, 0.50, 3.82, 4.72, 1.82, "none", "#7A8BA0", lw=1.15, rs=0.14, z=2, ls=(0, (3.5, 2.4)))
    for i, m in enumerate(["RF-Domain Classifier", "Pseudo-Outcome Risk", "MMD-LOCO"]):
        rbox(ax, 0.66, 5.00 - i * 0.54, 4.40, 0.46, ORANGE, ORANGE_DK, lw=1.1, rs=0.09)
        box_label(ax, 0.66, 5.00 - i * 0.54, 4.40, 0.46, m, fontsize=10.4)

    arrow(ax, (5.26, 4.73), (6.18, 4.73), lw=1.7, ms=12)

    ax.text(8.58, 5.78, "SELECTED INSTANCES", fontsize=8.15, fontweight="bold", color=HEADER_L2, zorder=5, ha="center")
    rbox(ax, 6.28, 4.82, 4.55, 0.72, ORANGE, ORANGE_DK, lw=1.2, rs=0.11)
    box_label(ax, 6.28, 4.82, 4.55, 0.72, "Sentence Selected", fontsize=12.2)
    rbox(ax, 6.28, 3.88, 4.55, 0.72, ORANGE, ORANGE_DK, lw=1.2, rs=0.11)
    box_label(ax, 6.28, 3.88, 4.55, 0.72, "Image Selected", fontsize=12.2)

    arrow(ax, (10.86, 5.18), (11.62, 4.95), lw=1.5, ms=11)
    arrow(ax, (10.86, 4.24), (11.62, 4.55), lw=1.5, ms=11)

    rbox(ax, 11.68, 3.88, 5.22, 1.72, WHITE, NAVY, lw=1.25, rs=0.12)
    ax.text(
        14.29,
        5.18,
        "Salient instances\nin the selected modality",
        ha="center",
        va="center",
        fontsize=11.0,
        fontweight="bold",
        color=NAVY_DK,
        zorder=6,
        linespacing=1.25,
    )
    ax.text(
        14.29,
        4.28,
        "Keep samples that drive the\nbatch discrepancy; drop the rest.",
        ha="center",
        va="center",
        fontsize=10.0,
        color=NAVY_DK,
        zorder=6,
        fontweight="normal",
        linespacing=1.35,
    )

    # ── LAYER 3 ──────────────────────────────────────────────────
    x, y, w, h = 0.28, 0.22, 16.94, 3.06
    rbox(ax, x, y, w, h, L3_BG, "#C5DCCE", lw=1.2, rs=0.16, z=1)
    header_bar(
        ax,
        x + 0.08,
        y + h - 0.50,
        w - 0.16,
        0.42,
        HEADER_L3,
        "LAYER 3",
        "Fine-grained Post-Hoc Localization   —   where inside the sample?",
    )
    ax.text(
        x + 0.28,
        y + h - 0.72,
        "Condition on the target concept. Localize tokens for text, spatial regions for images.",
        fontsize=9.5,
        color=MUTED,
        va="center",
        zorder=5,
    )

    ax.text(0.72, 2.18, "TARGET CONCEPT", fontsize=8.15, fontweight="bold", color=HEADER_L3, zorder=5)
    rbox(ax, 0.58, 0.98, 4.55, 1.02, PEACH, ORANGE_DK, lw=1.15, rs=0.12)
    box_label(ax, 0.58, 0.98, 4.55, 1.02, "Sentiment  ·  Image Category", fontsize=11.2, color=NAVY_DK)

    arrow(ax, (5.16, 1.49), (6.05, 1.49), lw=1.7, ms=12)

    ax.text(
        11.20,
        2.18,
        "FINE-GRAINED ATTRIBUTION",
        fontsize=8.15,
        fontweight="bold",
        color=HEADER_L3,
        zorder=5,
        ha="center",
    )
    rbox(ax, 6.12, 0.88, 10.20, 1.18, "none", "#7A8BA0", lw=1.05, rs=0.14, z=2, ls=(0, (3.5, 2.4)))
    rbox(ax, 6.28, 0.98, 4.72, 1.00, ORANGE, ORANGE_DK, lw=1.2, rs=0.12)
    box_label(ax, 6.28, 0.98, 4.72, 1.00, "Token-Level Localization\n(text / transcript)", fontsize=11.6)
    rbox(ax, 11.28, 0.98, 4.85, 1.00, ORANGE, ORANGE_DK, lw=1.2, rs=0.12)
    box_label(ax, 11.28, 0.98, 4.85, 1.00, "Bounding-Box Attribution\n(image / spatial regions)", fontsize=11.6)

    ax.text(
        W / 2,
        0.52,
        "Output:  token spans  ·  spatial boxes  ·  ranked drivers of the observed preference / mixture shift",
        ha="center",
        va="center",
        fontsize=10.1,
        color=HEADER_L3,
        fontweight="bold",
        zorder=6,
    )

    def vflow(y0, y1, note):
        ax.annotate(
            "",
            xy=(8.75, y1),
            xytext=(8.75, y0),
            arrowprops=dict(arrowstyle="-|>", color=NAVY, lw=1.85, mutation_scale=13),
            zorder=7,
        )
        if note:
            ax.text(8.95, (y0 + y1) / 2, note, fontsize=8.4, color=NAVY, va="center", fontstyle="italic", zorder=7)

    vflow(11.48, 11.28, "")
    vflow(7.28, 6.96, "selected modality + top features")
    vflow(3.62, 3.30, "selected instances + target concept")

    fig.savefig(out, dpi=230, bbox_inches="tight", facecolor=WHITE, pad_inches=0.16)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
