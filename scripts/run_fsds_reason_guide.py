#!/usr/bin/env python3
"""Scorecard: FSDS → adjust next step → ROI / incremental value."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.fsds_reason_guide import justify_bullets, run_suite


def render_md(suite: dict) -> str:
    b, g, inc = suite["baseline"], suite["guided"], suite["incremental"]
    top = suite["top_features"]
    lines = [
        "# FSDS Reasoning Empowerment — Adjust Next Step",
        "",
        "> Skill: `fsds_reasoning_empowerment` · Y=1{node on oracle path} · Guide={prune,expand,reorder,backtrack}",
        "",
        "## Verdict (simulation ledger)",
        "",
    ]
    for L in justify_bullets(suite):
        lines.append(f"- {L}")
    lines += [
        "",
        "## Regime & attribution",
        "",
        f"- **overlap** = {suite['overlap']:.4f}",
        f"- **regime** = `{suite['regime']}`",
        f"- **components** = {suite['attr']['components']}",
        f"- **AUC P(Y=1|X)** = {suite['auc_y']:.4f}",
        "",
        "### Top fused features",
        "",
        "| feature | importance |",
        "|---|---:|",
    ]
    for name, w in top:
        lines.append(f"| `{name}` | {w:.4f} |")
    lines += [
        "",
        "## Search effect (held-out tasks)",
        "",
        f"n_train={suite['n_train']} · n_test={suite['n_test']}",
        "",
        "| policy | success | tokens/task | cost/task | net/task | waste_off_opt |",
        "|---|---:|---:|---:|---:|---:|",
        f"| score-greedy | {b['success_rate']:.1%} | {b['tokens']:.1f} | {b['cost']:.3f} | {b['net']:.3f} | {b['waste_off_opt']:.3f} |",
        f"| FSDS-guided | {g['success_rate']:.1%} | {g['tokens']:.1f} | {g['cost']:.3f} | {g['net']:.3f} | {g['waste_off_opt']:.3f} |",
        f"| **Δ** | {inc['delta_success']:+.1%} | {inc['delta_tokens']:+.1f} | {b['cost']-g['cost']:+.3f} | {inc['delta_net_per_task']:+.3f} | {b['waste_off_opt']-g['waste_off_opt']:+.3f} |",
        "",
        "### Guide action counts (test)",
        "",
        "```",
        json.dumps(suite["action_counts"], indent=2),
        "```",
        "",
        "## Economic ledger (skill §7)",
        "",
        "```",
        "NetValue = Σ Y_i · V_task · N_tasks − Σ (1−Y_i) · C_node",
        "Incremental ≈ ΔRevenue + ΔCost_savings + ΔAUC·V_auc − Cost_FSDS",
        "```",
        "",
        f"| term | value / task |",
        f"|---|---:|",
        f"| ΔRevenue (via success) | {inc['delta_revenue']:+.4f} |",
        f"| ΔCost savings | {inc['delta_cost_savings']:+.4f} |",
        f"| ΔAUC · V_auc | {inc['delta_auc_value']:+.4f} |",
        f"| Cost_FSDS | −{inc['cost_fsds']:.4f} |",
        f"| **Net incremental** | **{inc['net_incremental']:+.4f}** |",
        f"| **ROI** (Δ / invest) | **{inc['roi']:.2f}** |",
        "",
        f"Params: `{suite['econ']}`",
        "",
        "## Justification",
        "",
        "1. **Why adjust next step with FSDS**: search nodes carry X (struct/semantic/score/…) "
        "and a labelable Y (on-oracle-path). FSDS predicts Y under shift; importance∝P(Y=1) "
        "is exactly the signal Guide needs for expand vs prune.",
        "2. **Why ROI is computable here**: success→V_task, off-opt nodes→C_node, tokens/calls "
        "metered; Δ vs score-greedy is an A/B-style ledger on the same trees.",
        "3. **What this does *not* claim**: dollar ROI in production RAP/LATS, or that fused "
        "importance is a causal CATE without further identification. Numbers are **simulation "
        "scorecard** under the skill's NetValue formula.",
        "4. **Incremental value**: positive when guided raises success and/or cuts wasted "
        "tokens enough to beat Cost_FSDS; see table above.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-train", type=int, default=40)
    ap.add_argument("--n-test", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-depth", type=int, default=4)
    ap.add_argument("--out", type=Path, default=Path("results/fsds_reason_guide"))
    args = ap.parse_args()
    suite = run_suite(
        n_train=args.n_train,
        n_test=args.n_test,
        seed=args.seed,
        max_depth=args.max_depth,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    # strip non-json bits already done in run_suite
    payload = {k: v for k, v in suite.items()}
    (args.out / "suite.json").write_text(json.dumps(payload, indent=2, default=str))
    md = render_md(suite)
    (args.out / "FSDS_REASON_NEXTSTEP_SCORECARD.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
