#!/usr/bin/env python3
"""Tomorrow-demo plot: PO rank_eff + √/∛ soft wins (no retrain).

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
        "--out",
        type=Path,
        default=ROOT / "results/agod_po_eff/rsi_tomorrow_demo.png",
    )
    args = ap.parse_args()
    eff = json.loads(args.eff.read_text())
    power = json.loads(args.power.read_text())

    names, probe_e, refit_e, duties = [], [], [], []
    for c in eff.get("cards") or []:
        if not c.get("ok"):
            continue
        by = {m["mode"]: m for m in c["modes"]}
        names.append(c["dataset"].replace("_", "\n"))
        probe_e.append(float(by.get("probe", {}).get("rank_eff", np.nan)))
        refit_e.append(float(by.get("refit", {}).get("rank_eff", np.nan)))
        duties.append(float(c.get("gate_duty", np.nan)))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    # Panel A: rank_eff probe vs refit
    ax = axes[0]
    x = np.arange(len(names))
    w = 0.38
    ax.bar(x - w / 2, probe_e, w, label="probe", color="#4C78A8")
    ax.bar(x + w / 2, refit_e, w, label="refit", color="#F58518")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("rank_eff (Δρ / M-FLOPs)")
    ax.set_title("A. Hard-rank per adapt FLOP\n(refit wins when duty is low)")
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0, color="gray", lw=0.6)

    # Panel B: IPTW softness
    ax = axes[1]
    best = power.get("best_sig_counts") or {}
    soft = float(power.get("soft_win_rate_cbrt_le_sqrt") or 0)
    labels = list(best.keys()) + ["∛≤√ rate"]
    vals = [best[k] for k in best.keys()] + [soft * sum(best.values())]
    colors = ["#54A24B", "#E45756", "#B279A2", "#72B7B2"][: len(labels)]
    ax.bar(range(len(labels)), vals, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("pack count (rate×n for ∛≤√)")
    ax.set_title(
        f"B. IPTW √ vs ∛ (same FLOPs)\n"
        f"∛≤√ on {soft:.0%} · uniform still often best"
    )

    fig.suptitle(
        f"RSI PO efficiency demo · mean duty="
        f"{eff.get('mean_gate_duty', float('nan')):.2f} · "
        f"E[refit]/E[probe]≈{eff.get('mean_budget_ratio_refit_vs_probe', float('nan')):.2f}",
        fontsize=11,
    )
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
