#!/usr/bin/env python3
"""Next-batch group LR from π. Two frozen batches — not online learning."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from clever_covariate_gap import (  # noqa: E402
    decompose_modality_gap,
    inject_block_shift,
    make_synthetic_shift,
)
from pi_next_batch import (  # noqa: E402
    block_learning_rates,
    gd_domain_path,
    next_batch_lr_eval,
)
from sklearn.model_selection import train_test_split  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402


def plot_paths(out_dir: Path, seed: int = 8) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    X, W, Y, spec = make_synthetic_shift(
        n=720, d_text=36, d_vad=4, gt="valence", mean_shift=1.15, seed=seed,
    )
    idx = np.arange(len(W))
    A, rest = train_test_split(idx, test_size=0.6, random_state=seed, stratify=W)
    B, C = train_test_split(rest, test_size=0.45, random_state=seed + 1, stratify=W[rest])
    gap = decompose_modality_gap(
        X[A], W[A], spec, seed=seed, n_splits=3, n_estimators=40, light=True,
    )
    pi = gap.pi_consensus
    scaler = StandardScaler()
    XB = scaler.fit_transform(X[B])
    XC = scaler.transform(X[C])
    n_steps = 28
    eta0 = 0.45
    curves = {
        "uniform": block_learning_rates(spec, pi, eta0=eta0, mode="uniform"),
        "π-boost": block_learning_rates(spec, pi, eta0=eta0, mode="boost"),
        "π-damp": block_learning_rates(spec, pi, eta0=eta0, mode="damp"),
    }
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    colors = {"uniform": "#9aa0a6", "π-boost": "#ff7f0e", "π-damp": "#1f77b4"}
    for name, eta in curves.items():
        aucs = gd_domain_path(XB, W[B], spec, eta, n_steps=n_steps, test=(XC, W[C]))
        ax.plot(np.arange(1, n_steps + 1), aucs, "o-", ms=3, label=name, color=colors[name])
    ax.set_xlabel("GD steps on the next batch (π frozen from batch A)")
    ax.set_ylabel("holdout two-sample AUC")
    ax.set_title("Group learning rates · not a stream, not online")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p = out_dir / "next_batch_lr_path.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    return str(p)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--quick", action="store_true")
    args = p.parse_args()
    n = 400 if args.quick else 720
    n_est = 25 if args.quick else 40
    n_steps = 12 if args.quick else 18
    n_seeds = 1 if args.quick else 3
    out_dir = ROOT / "results" / "pi_next_batch"
    out_dir.mkdir(parents=True, exist_ok=True)

    def synth(seed):
        return make_synthetic_shift(
            n=n, d_text=36 if not args.quick else 20, d_vad=4,
            gt="valence", mean_shift=1.15, seed=seed,
        )

    def inject(seed):
        X, W, Y, spec = make_synthetic_shift(
            n=n, d_text=36 if not args.quick else 20, d_vad=4,
            gt="text", mean_shift=0.12, seed=seed,
        )
        return inject_block_shift(X, W, spec, "valence", alpha=1.2), W, Y, spec

    board = {}
    for name, maker, gt in (
        ("synthetic GT=valence", synth, "valence"),
        ("inject valence", inject, "valence"),
    ):
        rows = [
            next_batch_lr_eval(
                *maker(7 + 11 * s), gt=gt, seed=7 + s,
                n_steps=n_steps, n_estimators=n_est,
            )
            for s in range(n_seeds)
        ]
        keys = rows[0]["domain_auc"].keys()
        board[name] = {
            "domain_auc": {
                k: round(float(np.mean([r["domain_auc"][k] for r in rows])), 4) for k in keys
            },
            "mass_on_gt": {
                k: round(float(np.mean([r["mass_on_gt"][k] for r in rows])), 4) for k in keys
            },
            "y_mse": {
                k: round(float(np.mean([r["y_mse"][k] for r in rows])), 4)
                for k in rows[0]["y_mse"]
            },
        }
    path = plot_paths(out_dir)
    lines = [
        "# Next-batch group learning rates from π",
        "",
        "Not online learning. Freeze π on batch A; finite-step GD on batch B with",
        "η_m = η0 M π_m (mean-preserving). Same simplex as packing / π-RF / adaptive ridge.",
        "",
        "| setting | uniform | π-boost | VIMP-boost | π-damp | oracle |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, row in board.items():
        a = row["domain_auc"]
        lines.append(
            f"| {name} | {a['uniform']:.3f} | {a['pi_boost']:.3f} | "
            f"{a['vimp_boost']:.3f} | {a['pi_damp']:.3f} | {a.get('oracle', float('nan')):.3f} |"
        )
    lines += [
        "",
        "Coefficient mass on GT after 18 GD steps:",
        "",
        "| setting | uniform | π-boost | VIMP-boost | π-damp |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, row in board.items():
        a = row["mass_on_gt"]
        lines.append(
            f"| {name} | {a['uniform']:.3f} | {a['pi_boost']:.3f} | "
            f"{a['vimp_boost']:.3f} | {a['pi_damp']:.3f} |"
        )
    lines += ["", f"Path plot: `{Path(path).name}`", ""]
    (out_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")
    (out_dir / "summary.json").write_text(json.dumps({"board": board, "path": path}, indent=2), encoding="utf-8")
    print(json.dumps({"board": board, "path": path}, indent=2))


if __name__ == "__main__":
    main()
