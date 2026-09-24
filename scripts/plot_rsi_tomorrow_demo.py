#!/usr/bin/env python3
"""Tomorrow-demo plot: PO rank_eff + √/∛ + freeze Pareto (one page).

  PYTHONPATH=. python3 scripts/plot_rsi_tomorrow_demo.py \\
    --out results/agod_po_eff/rsi_tomorrow_demo.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def _load(path: Path) -> dict:
    # Python may have written NaN; tolerate via a tiny replace
    text = path.read_text().replace("NaN", "null")
    return json.loads(text)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--eff",
        type=Path,
        default=ROOT / "results/agod_po_eff/po_eff_scorecard.json",
    )
    ap.add_argument(
        "--power",
        type=Path,
        default=ROOT / "results/agod_po_power_eff/po_power_eff.json",
    )
    ap.add_argument(
        "--freeze",
        type=Path,
        default=ROOT / "results/agod_freeze_eff/freeze_eff.json",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results/agod_po_eff/rsi_tomorrow_demo.png",
    )
    args = ap.parse_args()
    eff = _load(args.eff)
    power = _load(args.power)
    freeze = _load(args.freeze) if args.freeze.exists() else {}

    names, probe_e, refit_e = [], [], []
    for c in eff.get("cards") or []:
        if not c.get("ok"):
            continue
        by = {m["mode"]: m for m in c["modes"]}
        names.append(c["dataset"].replace("_", "\n"))
        probe_e.append(float(by.get("probe", {}).get("rank_eff", np.nan)))
        refit_e.append(float(by.get("refit", {}).get("rank_eff", np.nan)))

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))

    # A: rank_eff
    ax = axes[0]
    x = np.arange(len(names))
    w = 0.38
    ax.bar(x - w / 2, probe_e, w, label="probe", color="#4C78A8")
    ax.bar(x + w / 2, refit_e, w, label="refit", color="#F58518")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=7)
    ax.set_ylabel("rank_eff (Δρ / M-FLOPs)")
    ax.set_title("A. Hard-rank / adapt FLOP\nrefit wins at low duty")
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0, color="gray", lw=0.6)

    # B: IPTW softness
    ax = axes[1]
    best = power.get("best_sig_counts") or {}
    soft = float(power.get("soft_win_rate_cbrt_le_sqrt") or 0)
    labels = list(best.keys()) + ["∛≤√"]
    n_packs = sum(best.values()) or 1
    vals = [best[k] for k in best.keys()] + [soft * n_packs]
    colors = ["#54A24B", "#E45756", "#B279A2", "#72B7B2"][: len(labels)]
    ax.bar(range(len(labels)), vals, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("pack count")
    ax.set_title(f"B. IPTW √ vs ∛ (same FLOPs)\n∛≤√ {soft:.0%} · uniform often best")

    # C: freeze Pareto scatter
    ax = axes[2]
    markers = {"electricity": "o", "synthetic": "s"}
    for c in freeze.get("cards") or []:
        ds = c.get("dataset", "?")
        for r in c.get("policies") or []:
            if r.get("policy") in ("always_adapt", "no_adapt"):
                continue
            mr, fr = r.get("mse_rel"), r.get("flops_rel")
            if mr is None or fr is None:
                continue
            col = "#54A24B" if r.get("dominates_always") else "#E45756"
            ax.scatter(
                fr,
                mr,
                s=70,
                c=col,
                marker=markers.get(ds, "D"),
                zorder=3,
                label=f"{ds}:{r['policy'][:10]}",
            )
    ax.axhline(1.0, color="gray", ls="--", lw=0.8)
    ax.axvline(1.0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("flops_rel vs always_adapt")
    ax.set_ylabel("mse_rel vs always_adapt")
    ax.set_title("C. Freeze Pareto\ngreen = MSE↓ & FLOPs↓")
    ax.set_xlim(0.4, 1.15)
    ax.set_ylim(0.7, 1.4)
    # de-dupe legend
    handles, labs = ax.get_legend_handles_labels()
    by_label = dict(zip(labs, handles))
    ax.legend(
        by_label.values(),
        by_label.keys(),
        fontsize=6,
        frameon=False,
        loc="upper right",
    )
    # shade Pareto quadrant
    ax.fill_between([0.4, 1.0], 0.7, 1.0, color="#54A24B", alpha=0.08, zorder=0)

    fig.suptitle(
        f"RSI efficiency demo · duty≈{eff.get('mean_gate_duty', float('nan')):.2f} · "
        f"E[refit]/E[probe]≈{eff.get('mean_budget_ratio_refit_vs_probe', float('nan')):.2f} · "
        f"freeze Pareto wins={freeze.get('n_pareto_dominances', 0)}",
        fontsize=10,
    )
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
