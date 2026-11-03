#!/usr/bin/env python3
"""Plot streaming C-GCD performance curves: with weights vs without.

  python3 scripts/plot_streaming_performance.py
  python3 scripts/plot_streaming_performance.py --quick
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from category_shift_boards import cold_start_efficiency, run_cgcd  # noqa: E402

OUT = ROOT / "results" / "repr_sep_cold_start"
FIGS = OUT / "figures"
ARTIFACTS = Path("/opt/cursor/artifacts")
DOCS = ROOT / "docs" / "method"

# with-weights vs no-weights contrast
MODES = ("uniform", "weighted", "repr_sep")
COLORS = {
    "uniform": "#6B7280",      # no weights
    "weighted": "#2563EB",     # sample-level weights
    "repr_sep": "#DC2626",     # opportunity-gated sep + weights at expand
}
LABELS = {
    "uniform": "no weights (uniform)",
    "weighted": "with sample weights",
    "repr_sep": "repr-sep + weights",
}
MARKERS = {"uniform": "o", "weighted": "s", "repr_sep": "D"}


def _collect(seeds, *, quick: bool, expand_min: int) -> dict:
    """Per-stage curves + efficiency curves over seeds."""
    stage_metrics = {m: {"novel": [], "known": [], "forget": [], "O_local": []} for m in MODES}
    eff_curves = {m: [] for m in MODES}

    for seed in seeds:
        for mode in MODES:
            r = run_cgcd(
                seed=int(seed), mode=mode, quick=bool(quick), expand_min=int(expand_min),
            )
            stages = r["stages"]
            stage_metrics[mode]["novel"].append([s.get("novel_acc", np.nan) for s in stages])
            stage_metrics[mode]["known"].append([s.get("known_acc", np.nan) for s in stages])
            stage_metrics[mode]["forget"].append([s.get("forgetting", np.nan) for s in stages])
            stage_metrics[mode]["O_local"].append([s.get("O_local", np.nan) for s in stages])
        eff = cold_start_efficiency(seed=int(seed), quick=bool(quick), target=0.60)
        for mode in MODES:
            pts = eff["curves"][mode]["points"]
            eff_curves[mode].append([(p["budget"], p["novel_acc"]) for p in pts])

    def stack(metric_lists):
        arr = np.asarray(metric_lists, float)  # (n_seeds, n_stages)
        if arr.size == 0:
            return {"mean": [], "std": [], "seeds": []}
        # all-NaN columns (e.g. O_local on uniform) → nan mean without warn
        with np.errstate(all="ignore"):
            mu = np.nanmean(arr, axis=0)
            sd = np.nanstd(arr, axis=0)
        return {
            "mean": np.where(np.isnan(mu), 0.0, mu).tolist() if False else mu.tolist(),
            "std": sd.tolist(),
            "seeds": arr.tolist(),
        }

    out_stages = {}
    for mode in MODES:
        out_stages[mode] = {
            "novel": stack(stage_metrics[mode]["novel"]),
            "known": stack(stage_metrics[mode]["known"]),
            "forget": stack(stage_metrics[mode]["forget"]),
            "O_local": stack(stage_metrics[mode]["O_local"]),
        }

    # efficiency: align budgets
    out_eff = {}
    for mode in MODES:
        # assume same budgets across seeds
        budgets = [b for b, _ in eff_curves[mode][0]]
        mat = np.asarray([[acc for _, acc in run] for run in eff_curves[mode]], float)
        out_eff[mode] = {
            "budgets": budgets,
            "mean": np.nanmean(mat, axis=0).tolist(),
            "std": np.nanstd(mat, axis=0).tolist(),
        }

    return {
        "seeds": list(seeds),
        "quick": bool(quick),
        "stages": out_stages,
        "efficiency": out_eff,
        "n_stages": len(out_stages["uniform"]["novel"]["mean"]),
    }


def _style_axes(ax, title, ylabel, xlabel):
    ax.set_title(title, fontsize=12, pad=8)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.grid(True, alpha=0.25, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_stage_curves(data: dict, path: Path) -> None:
    n = data["n_stages"]
    xs = np.arange(n)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6), constrained_layout=True)

    specs = [
        ("novel", "Novel accuracy (streaming)", "novel_acc"),
        ("known", "Known retention", "known_acc"),
        ("forget", "Forgetting on pinned REF", "forgetting"),
    ]
    for ax, (key, title, ylab) in zip(axes, specs):
        for mode in MODES:
            mu = np.asarray(data["stages"][mode][key]["mean"], float)
            sd = np.asarray(data["stages"][mode][key]["std"], float)
            ax.plot(
                xs, mu, marker=MARKERS[mode], color=COLORS[mode],
                label=LABELS[mode], linewidth=2.0, markersize=7,
            )
            ax.fill_between(xs, mu - sd, mu + sd, color=COLORS[mode], alpha=0.15)
        ax.set_xticks(xs)
        ax.set_xticklabels([f"s{i}" for i in xs])
        ax.set_ylim(0.0, 1.05)
        _style_axes(ax, title, ylab, "streaming stage")
        if key == "novel":
            ax.legend(frameon=False, fontsize=8, loc="best")

    fig.suptitle(
        "Streaming C-GCD proxy — with weights vs no weights",
        fontsize=13, fontweight="bold", y=1.02,
    )
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_efficiency(data: dict, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    for mode in MODES:
        budgets = np.asarray(data["efficiency"][mode]["budgets"], float)
        mu = np.asarray(data["efficiency"][mode]["mean"], float)
        sd = np.asarray(data["efficiency"][mode]["std"], float)
        ax.plot(
            budgets, mu, marker=MARKERS[mode], color=COLORS[mode],
            label=LABELS[mode], linewidth=2.0, markersize=7,
        )
        ax.fill_between(budgets, mu - sd, mu + sd, color=COLORS[mode], alpha=0.15)
    ax.axhline(0.60, color="#9CA3AF", linestyle=":", linewidth=1.2, label="target 0.60")
    ax.set_ylim(0.0, 1.05)
    _style_axes(ax, "Cold-start efficiency: novel budget → novel_acc", "novel_acc", "n novel samples")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_delta_vs_uniform(data: dict, path: Path) -> None:
    """Δnovel vs no-weights at each stage."""
    n = data["n_stages"]
    xs = np.arange(n)
    base = np.asarray(data["stages"]["uniform"]["novel"]["mean"], float)
    fig, ax = plt.subplots(figsize=(6.2, 3.8), constrained_layout=True)
    width = 0.35
    for i, mode in enumerate(("weighted", "repr_sep")):
        mu = np.asarray(data["stages"][mode]["novel"]["mean"], float) - base
        ax.bar(
            xs + (i - 0.5) * width, mu, width=width,
            color=COLORS[mode], label=LABELS[mode], alpha=0.9,
        )
    ax.axhline(0.0, color="#374151", linewidth=1.0)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"s{i}" for i in xs])
    _style_axes(ax, "Δ novel_acc vs no weights (uniform)", "Δ novel_acc", "streaming stage")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_combo(data: dict, path: Path) -> None:
    """One-pager: stage novel + efficiency + delta."""
    fig = plt.figure(figsize=(11.5, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)

    # (0,0)-(0,1): novel streaming curve wide
    ax0 = fig.add_subplot(gs[0, :])
    n = data["n_stages"]
    xs = np.arange(n)
    for mode in MODES:
        mu = np.asarray(data["stages"][mode]["novel"]["mean"], float)
        sd = np.asarray(data["stages"][mode]["novel"]["std"], float)
        ax0.plot(
            xs, mu, marker=MARKERS[mode], color=COLORS[mode],
            label=LABELS[mode], linewidth=2.2, markersize=8,
        )
        ax0.fill_between(xs, mu - sd, mu + sd, color=COLORS[mode], alpha=0.14)
    ax0.set_xticks(xs)
    ax0.set_xticklabels([f"stage {i}" for i in xs])
    ax0.set_ylim(0.0, 1.05)
    _style_axes(ax0, "Streaming novel_acc — weights vs no weights", "novel_acc", "")
    ax0.legend(frameon=False, fontsize=9, loc="best")

    ax1 = fig.add_subplot(gs[1, 0])
    for mode in MODES:
        budgets = np.asarray(data["efficiency"][mode]["budgets"], float)
        mu = np.asarray(data["efficiency"][mode]["mean"], float)
        sd = np.asarray(data["efficiency"][mode]["std"], float)
        ax1.plot(
            budgets, mu, marker=MARKERS[mode], color=COLORS[mode],
            label=LABELS[mode], linewidth=2.0, markersize=6,
        )
        ax1.fill_between(budgets, mu - sd, mu + sd, color=COLORS[mode], alpha=0.14)
    ax1.axhline(0.60, color="#9CA3AF", linestyle=":", linewidth=1.2)
    ax1.set_ylim(0.0, 1.05)
    _style_axes(ax1, "Efficiency curve (budget → novel_acc)", "novel_acc", "n novel samples")

    ax2 = fig.add_subplot(gs[1, 1])
    base = np.asarray(data["stages"]["uniform"]["novel"]["mean"], float)
    width = 0.35
    for i, mode in enumerate(("weighted", "repr_sep")):
        mu = np.asarray(data["stages"][mode]["novel"]["mean"], float) - base
        ax2.bar(
            xs + (i - 0.5) * width, mu, width=width,
            color=COLORS[mode], label=LABELS[mode], alpha=0.9,
        )
    ax2.axhline(0.0, color="#374151", linewidth=1.0)
    ax2.set_xticks(xs)
    ax2.set_xticklabels([f"s{i}" for i in xs])
    _style_axes(ax2, "Δ novel vs no weights", "Δ novel_acc", "streaming stage")
    ax2.legend(frameon=False, fontsize=8)

    fig.suptitle(
        "Streaming category-shift cold-start — performance curves",
        fontsize=14, fontweight="bold",
    )
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--expand-min", type=int, default=12)
    ap.add_argument("--n-seeds", type=int, default=3)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    seeds = [args.seed + i for i in range(args.n_seeds)]
    data = _collect(seeds, quick=bool(args.quick), expand_min=args.expand_min)
    (OUT / "streaming_curves.json").write_text(json.dumps(data, indent=2))

    p_stage = FIGS / "streaming_stage_curves.png"
    p_eff = FIGS / "streaming_efficiency_curves.png"
    p_delta = FIGS / "streaming_delta_vs_noweight.png"
    p_combo = FIGS / "streaming_performance_curves.png"
    plot_stage_curves(data, p_stage)
    plot_efficiency(data, p_eff)
    plot_delta_vs_uniform(data, p_delta)
    plot_combo(data, p_combo)

    # copy hero figures to artifacts + docs
    DOCS.mkdir(parents=True, exist_ok=True)
    try:
        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        for p in (p_stage, p_eff, p_delta, p_combo):
            shutil.copy2(p, ARTIFACTS / p.name)
            shutil.copy2(p, DOCS / p.name)
        shutil.copy2(OUT / "streaming_curves.json", ARTIFACTS / "streaming_curves.json")
    except Exception as e:
        print("artifact copy warn:", e)

    md = "\n".join([
        "# Streaming performance curves — weights vs no weights",
        "",
        f"Seeds {seeds}, quick={args.quick}.",
        "",
        "## Novel / known / forgetting by stage",
        "",
        f'<img alt="stage curves" src="{p_stage}" width="900" />',
        "",
        "## One-pager",
        "",
        f'<img alt="combo" src="{p_combo}" width="900" />',
        "",
        "## Efficiency (budget → novel_acc)",
        "",
        f'<img alt="efficiency" src="{p_eff}" width="560" />',
        "",
        "## Δ novel vs no weights",
        "",
        f'<img alt="delta" src="{p_delta}" width="560" />',
        "",
        "```bash",
        "python3 scripts/plot_streaming_performance.py",
        "```",
        "",
    ])
    (OUT / "STREAMING_CURVES.md").write_text(md)
    (DOCS / "Streaming_performance_curves.md").write_text(md)
    print(md)
    print("wrote", p_combo)


if __name__ == "__main__":
    main()
