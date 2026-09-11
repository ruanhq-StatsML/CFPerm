#!/usr/bin/env python3
"""Same-expert + context-packing board (two-sample, no TMLE, no vendor LLM)."""
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
    make_diffuse_equal_shift,
    worked_example_packet,
)


def _mean_eval(maker, *, n_seeds: int, **kw) -> dict:
    rows = []
    for s in range(n_seeds):
        X, W, Y, spec = maker(seed=kw.get("base_seed", 7) + 17 * s)
        rows.append(
            budgeted_inference_eval(
                X, W, spec,
                gt=kw.get("gt"),
                seed=kw.get("base_seed", 7) + s,
                n_estimators=kw.get("n_estimators", 50),
                test_size=0.35,
            )
        )
    keys = rows[0]["auc"].keys()
    auc = {k: round(float(np.mean([r["auc"][k] for r in rows])), 4) for k in keys}
    hit = {}
    if kw.get("gt") is not None and rows[0].get("hit_gt"):
        hit_keys = rows[0]["hit_gt"].keys()
        hit = {k: round(float(np.mean([r["hit_gt"][k] for r in rows])), 4) for k in hit_keys}
    pack_keys = rows[0]["pack_auc"].keys()
    pack_auc = {
        k: round(float(np.mean([r["pack_auc"][k] for r in rows if k in r["pack_auc"]])), 4)
        for k in pack_keys
    }
    modes = [r["population_pack"]["mode"] for r in rows]
    bkeys = [
        "adapt_prior_auc", "adapt_shrink_auc", "adapt_query_only_auc",
        "adapt_prior_hit", "adapt_shrink_hit", "adapt_query_only_hit",
        "delta_auc_shrink_minus_prior", "delta_hit_shrink_minus_prior",
        "switch_rate",
    ]
    stage_B = {}
    for k in bkeys:
        vals = [
            r["stage_B_adapt"][k]
            for r in rows
            if r["stage_B_adapt"].get(k) is not None
            and r["stage_B_adapt"][k] == r["stage_B_adapt"][k]
        ]
        stage_B[k] = round(float(np.mean(vals)), 4) if vals else None
    return {
        "n_seeds": n_seeds,
        "gt": kw.get("gt"),
        "auc": auc,
        "hit_gt": hit,
        "pack_auc": pack_auc,
        "pi_entropy_last": rows[-1]["pi_entropy"],
        "population_pack_modes": modes,
        "pi_last": rows[-1]["pi"],
        "pi_vimp_last": rows[-1]["pi_vimp"],
        "selected_last": rows[-1]["selected_block"],
        "stage_A_last": rows[-1]["stage_A_block"],
        "stage_B": stage_B,
    }


def plot_board(board: dict, out_dir: Path) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    settings = list(board["settings"].keys())
    labs = settings
    x = np.arange(len(labs))

    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    w = 0.15
    series = [
        ("oracle_gt", "oracle GT specialist", "#2ca02c"),
        ("pi_top1", "same expert (π)", "#ff7f0e"),
        ("vimp_top1", "same expert (VIMP)", "#1f77b4"),
        ("instance_top1", "instance expert s(x)", "#d62728"),
        ("random_fixed", "same random expert", "#9aa0a6"),
    ]
    for i, (key, lab, c) in enumerate(series):
        vals = [board["settings"][s]["auc"].get(key, np.nan) for s in settings]
        ax.bar(x + (i - 2) * w, vals, w, label=lab, color=c)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=12, ha="right")
    ax.set_ylabel("held-out two-sample AUC")
    ax.set_ylim(0.45, 1.02)
    ax.axhline(0.5, ls="--", lw=0.8, color="0.5")
    ax.set_title("Same expert: one specialist ê_m for every query (m from π, not from x)")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p = out_dir / "budgeted_auc.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths["auc"] = p

    hit_settings = [s for s in settings if board["settings"][s].get("hit_gt")]
    if hit_settings:
        fig, ax = plt.subplots(figsize=(8.4, 4.3))
        xh = np.arange(len(hit_settings))
        hit_series = [
            ("pi_top1", "π (same expert)", "#ff7f0e"),
            ("vimp_top1", "VIMP (same expert)", "#1f77b4"),
            ("instance_top1", "instance s(x)", "#d62728"),
            ("random_fixed", "same random", "#9aa0a6"),
        ]
        w = 0.18
        for i, (key, lab, c) in enumerate(hit_series):
            vals = [board["settings"][s]["hit_gt"].get(key, np.nan) for s in hit_settings]
            ax.bar(xh + (i - 1.5) * w, vals, w, label=lab, color=c)
        ax.set_xticks(xh)
        ax.set_xticklabels(hit_settings, rotation=12, ha="right")
        ax.set_ylabel("P(selected block = GT)")
        ax.set_ylim(0.0, 1.05)
        ax.set_title("Localization: same-expert hit vs per-row instance hit")
        ax.legend(frameon=False, fontsize=8, ncol=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        p = out_dir / "gt_hit_rate.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        paths["hit"] = p

    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    pack_series = [
        ("same_expert_oracle", "pack GT block", "#2ca02c"),
        ("same_expert_pi", "pack π block", "#ff7f0e"),
        ("same_expert_vimp", "pack VIMP block", "#1f77b4"),
        ("concat_top2_pi", "concat top-2(π)", "#9467bd"),
        ("pack_all", "pack all (abstain-from-drop)", "#111111"),
    ]
    w = 0.15
    for i, (key, lab, c) in enumerate(pack_series):
        vals = [board["settings"][s]["pack_auc"].get(key, np.nan) for s in settings]
        ax.bar(x + (i - 2) * w, vals, w, label=lab, color=c)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, rotation=12, ha="right")
    ax.set_ylabel("held-out AUC of one RF on concatenated columns")
    ax.set_ylim(0.45, 1.02)
    ax.axhline(0.5, ls="--", lw=0.8, color="0.5")
    ax.set_title("Context packing: one reader on concat(E), E = E(π) for all rows")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p = out_dir / "packed_auc.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    paths["pack"] = p

    if any(board["settings"][s].get("stage_B") for s in settings):
        fig, ax = plt.subplots(figsize=(8.8, 4.3))
        w = 0.22
        qseries = [
            ("adapt_prior_auc", "prior E(π)", "#ff7f0e"),
            ("adapt_shrink_auc", "query update λ=0.4", "#111111"),
            ("adapt_query_only_auc", "query-only λ=0", "#d62728"),
        ]
        for i, (key, lab, c) in enumerate(qseries):
            vals = [board["settings"][s].get("stage_B", {}).get(key, np.nan) for s in settings]
            ax.bar(x + (i - 1) * w, vals, w, label=lab, color=c)
        ax.set_xticks(x)
        ax.set_xticklabels(labs, rotation=12, ha="right")
        ax.set_ylabel("held-out AUC after choosing E")
        ax.set_ylim(0.45, 1.02)
        ax.axhline(0.5, ls="--", lw=0.8, color="0.5")
        ax.set_title("Query update of E (π and f stay frozen)")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        p = out_dir / "query_update_auc.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        paths["query"] = p
    return {k: str(v) for k, v in paths.items()}


def write_readme(board: dict, out_dir: Path) -> None:
    lines = [
        "# Gap-guided reading budget (not TMLE, not causal)",
        "",
        "Two-sample localization of W on X. **Same expert**: one m̂=argmax π for every query.",
        "**Context packing**: one reader on concatenated columns of E=E(π). Not an opinion pool.",
        "",
        "## Same-expert specialist AUC (ê_m, one block, all rows)",
        "",
        "| setting | oracle | π | VIMP | instance s(x) | same random |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, row in board["settings"].items():
        a = row["auc"]
        lines.append(
            f"| {name} | {a.get('oracle_gt', float('nan')):.3f} | "
            f"{a.get('pi_top1', float('nan')):.3f} | {a.get('vimp_top1', float('nan')):.3f} | "
            f"{a.get('instance_top1', float('nan')):.3f} | {a.get('random_fixed', float('nan')):.3f} |"
        )
    lines += [
        "",
        "## Context packing AUC (one RF on concat E)",
        "",
        "| setting | pack GT | pack π | pack VIMP | concat top-2(π) | pack all |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, row in board["settings"].items():
        p = row["pack_auc"]
        lines.append(
            f"| {name} | {p.get('same_expert_oracle', float('nan')):.3f} | "
            f"{p.get('same_expert_pi', float('nan')):.3f} | {p.get('same_expert_vimp', float('nan')):.3f} | "
            f"{p.get('concat_top2_pi', float('nan')):.3f} | {p.get('pack_all', float('nan')):.3f} |"
        )
    lines += [
        "",
        "## Stage B: query update of E (π, f frozen)",
        "",
        "| setting | prior AUC | shrink AUC | query-only AUC | switch rate | Δ AUC |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, row in board["settings"].items():
        b = row.get("stage_B") or {}
        lines.append(
            f"| {name} | {b.get('adapt_prior_auc', float('nan')):.3f} | "
            f"{b.get('adapt_shrink_auc', float('nan')):.3f} | "
            f"{b.get('adapt_query_only_auc', float('nan')):.3f} | "
            f"{b.get('switch_rate', float('nan')):.3f} | "
            f"{b.get('delta_auc_shrink_minus_prior', float('nan')):+.3f} |"
        )
    lines += [
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

    def synth_diffuse(seed: int):
        return make_diffuse_equal_shift(
            n=n, d_text=d_text, d_vad=4, mean_shift=0.95, seed=seed,
        )

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
            "diffuse equal shift": _mean_eval(
                synth_diffuse, n_seeds=n_seeds, gt=None,
                base_seed=args.seed + 6, n_estimators=n_est,
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
    (out_dir / "example_packed_user.txt").write_text(pkt["user_message"], encoding="utf-8")
    paths = plot_board(board, out_dir)
    write_readme(board, out_dir)
    summary = {
        "board": board,
        "example_row": pkt["row_index"],
        "example_pack_mode": pkt["decision"].get("pack_mode"),
        "example_packed_blocks": pkt["decision"].get("packed_blocks"),
        "plots": paths,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({
        "plots": paths,
        "pack_mode": pkt["decision"].get("pack_mode"),
        "packed_blocks": pkt["decision"].get("packed_blocks"),
        "board_pack": {k: v["pack_auc"] for k, v in board["settings"].items()},
        "board_auc": {k: v["auc"] for k, v in board["settings"].items()},
        "board_hit": {k: v["hit_gt"] for k, v in board["settings"].items()},
        "board_stage_B": {k: v.get("stage_B") for k, v in board["settings"].items()},
    }, indent=2))


if __name__ == "__main__":
    main()
