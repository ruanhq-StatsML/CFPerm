#!/usr/bin/env python3
"""Budgeted-block stand-in for gap-guided LLM routing (no vendor API)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from clever_covariate_gap import inject_block_shift, make_synthetic_shift, shuffle_block  # noqa: E402
from gap_guided_inference import (  # noqa: E402
    budgeted_inference_eval,
    worked_example_packet,
)


def _mean_eval(maker, *, n_seeds: int, **kw) -> dict:
    rows = []
    for s in range(n_seeds):
        X, W, Y, spec = maker(seed=kw.get("base_seed", 7) + 17 * s)
        rows.append(
            budgeted_inference_eval(
                X, W, spec,
                gt=kw["gt"],
                seed=kw.get("base_seed", 7) + s,
                n_estimators=kw.get("n_estimators", 50),
                test_size=0.35,
            )
        )
    keys = rows[0]["auc"].keys()
    auc = {k: round(float(np.mean([r["auc"][k] for r in rows])), 4) for k in keys}
    hit_keys = rows[0]["hit_gt"].keys()
    hit = {k: round(float(np.mean([r["hit_gt"][k] for r in rows])), 4) for k in hit_keys}
    return {
        "n_seeds": n_seeds,
        "gt": kw["gt"],
        "auc": auc,
        "hit_gt": hit,
        "pi_last": rows[-1]["pi"],
        "pi_vimp_last": rows[-1]["pi_vimp"],
        "selected_last": rows[-1]["selected_block"],
    }


def plot_board(board: dict, out_dir: Path) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    settings = list(board["settings"].keys())
    # AUC: blend / pi / vimp / random / oracle
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    labs = settings
    x = np.arange(len(labs))
    w = 0.15
    series = [
        ("oracle_gt", "oracle GT block", "#2ca02c"),
        ("blend_top1", "blend r top-1", "#111111"),
        ("pi_top1", "π top-1", "#ff7f0e"),
        ("vimp_top1", "VIMP top-1", "#1f77b4"),
        ("random_row", "random block", "#9aa0a6"),
    ]
    for i, (key, lab, c) in enumerate(series):
        vals = [board["settings"][s]["auc"].get(key, np.nan) for s in settings]
        ax.bar(x + (i - 2) * w, vals, w, label=lab, color=c)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=12, ha="right")
    ax.set_ylabel("held-out domain AUC (one-block budget)")
    ax.set_ylim(0.45, 1.02)
    ax.axhline(0.5, ls="--", lw=0.8, color="0.5")
    ax.set_title("Budgeted inference · routed block expert vs VIMP / random")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p = out_dir / "budgeted_auc.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths["auc"] = p

    fig, ax = plt.subplots(figsize=(8.4, 4.3))
    hit_series = [
        ("blend_top1", "blend r", "#111111"),
        ("pi_top1", "π", "#ff7f0e"),
        ("instance_top1", "instance s", "#d62728"),
        ("vimp_top1", "VIMP", "#1f77b4"),
        ("random_row", "random", "#9aa0a6"),
    ]
    w = 0.15
    for i, (key, lab, c) in enumerate(hit_series):
        vals = [board["settings"][s]["hit_gt"][key] for s in settings]
        ax.bar(x + (i - 2) * w, vals, w, label=lab, color=c)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=12, ha="right")
    ax.set_ylabel("P(selected block = GT)")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Channel hit rate under a one-tool budget")
    ax.legend(frameon=False, fontsize=8, ncol=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p = out_dir / "gt_hit_rate.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths["hit"] = p
    return {k: str(v) for k, v in paths.items()}


def write_readme(board: dict, out_dir: Path) -> None:
    lines = [
        "# Gap-guided LLM inference (budgeted-block stand-in)",
        "",
        "Frozen π / s_m(x) route a one-block expert — the same decision as enabling one tool or packing one modality into the prompt.",
        "No vendor LLM is called. Wide text + concentrated valence shift is the regime where raw VIMP overweights text.",
        "",
        "## Held-out domain AUC (one-block budget, mean over seeds)",
        "",
        "| setting | oracle | blend | π | VIMP | random |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, row in board["settings"].items():
        a = row["auc"]
        lines.append(
            f"| {name} | {a.get('oracle_gt', float('nan')):.3f} | "
            f"{a['blend_top1']:.3f} | {a['pi_top1']:.3f} | "
            f"{a['vimp_top1']:.3f} | {a['random_row']:.3f} |"
        )
    lines += [
        "",
        "## P(selected block = GT)",
        "",
        "| setting | blend | π | instance | VIMP | random |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, row in board["settings"].items():
        h = row["hit_gt"]
        lines.append(
            f"| {name} | {h['blend_top1']:.3f} | {h['pi_top1']:.3f} | "
            f"{h['instance_top1']:.3f} | {h['vimp_top1']:.3f} | {h['random_row']:.3f} |"
        )
    lines += [
        "",
        "Copy-paste prompt: `example_system_prompt.txt`. Packet JSON: `example_packet.json`.",
        "",
        "```bash",
        "python3 scripts/run_gap_guided_inference.py",
        "```",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--quick", action="store_true")
    p.add_argument("--seed", type=int, default=2026)
    args = p.parse_args()
    n = 240 if args.quick else 520
    d_text = 16 if args.quick else 28
    n_est = 30 if args.quick else 55
    n_seeds = 1 if args.quick else 3
    out_dir = ROOT / "results" / "gap_guided_inference"
    out_dir.mkdir(parents=True, exist_ok=True)

    def synth_valence(seed: int):
        return make_synthetic_shift(
            n=n, d_text=d_text, d_vad=4, gt="valence", mean_shift=1.25, seed=seed,
        )

    def synth_inject(seed: int):
        X, W, Y, spec = make_synthetic_shift(
            n=n, d_text=d_text, d_vad=4, gt="text", mean_shift=0.12, seed=seed,
        )
        X = inject_block_shift(X, W, spec, "valence", alpha=1.15)
        return X, W, Y, spec

    def synth_shuffle_neg(seed: int):
        X, W, Y, spec = make_synthetic_shift(
            n=n, d_text=d_text, d_vad=4, gt="valence", mean_shift=1.25, seed=seed,
        )
        X = shuffle_block(X, spec, "valence", seed=seed + 3)
        return X, W, Y, spec

    board = {
        "settings": {
            "synthetic GT=valence": _mean_eval(
                synth_valence, n_seeds=n_seeds, gt="valence",
                base_seed=args.seed, n_estimators=n_est,
            ),
            "inject valence": _mean_eval(
                synth_inject, n_seeds=n_seeds, gt="valence",
                base_seed=args.seed + 3, n_estimators=n_est,
            ),
            "shuffle valence (neg.)": _mean_eval(
                synth_shuffle_neg, n_seeds=n_seeds, gt="valence",
                base_seed=args.seed + 9, n_estimators=n_est,
            ),
        }
    }
    X, W, Y, spec = synth_valence(args.seed)
    pkt = worked_example_packet(X, W, spec, seed=args.seed, n_estimators=n_est)
    (out_dir / "example_packet.json").write_text(
        json.dumps({k: pkt[k] for k in pkt if k != "system_prompt"}, indent=2),
        encoding="utf-8",
    )
    (out_dir / "example_system_prompt.txt").write_text(pkt["system_prompt"], encoding="utf-8")
    paths = plot_board(board, out_dir)
    write_readme(board, out_dir)
    summary = {
        "board": board,
        "example_row": pkt["row_index"],
        "example_tools": pkt["openai_tools"],
        "example_must_cite": pkt["decision"]["critic_must_cite"],
        "plots": paths,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({
        "plots": paths,
        "must_cite": pkt["decision"]["critic_must_cite"],
        "board_auc": {k: v["auc"] for k, v in board["settings"].items()},
        "board_hit": {k: v["hit_gt"] for k, v in board["settings"].items()},
    }, indent=2))


if __name__ == "__main__":
    main()
