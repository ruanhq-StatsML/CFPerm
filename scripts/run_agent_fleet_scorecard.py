#!/usr/bin/env python3
"""180-agent fleet scorecard: even force vs pile-on + CUPED/fraud CL.

  PYTHONPATH=. python3 scripts/run_agent_fleet_scorecard.py \\
    --out results/agod_agent_fleet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

import numpy as np

from agod.agent_fleet import (
    DEFAULT_COHORTS,
    compare_even_vs_pile,
    fleet_observation_spec,
    multi_tick_cl_demo,
)

ROOT = Path(__file__).resolve().parents[1]


def _clean(o):
    if isinstance(o, float) and (o != o or o in (float("inf"), float("-inf"))):
        return None
    if isinstance(o, (np.floating,)):
        v = float(o)
        return None if v != v else v
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    return o


def _fmt(x: Any, nd: int = 3) -> str:
    try:
        v = float(x)
        if v != v:
            return "—"
        return f"{v:.{nd}g}"
    except Exception:
        return "—"


def render_md(rep: Dict[str, Any], cl: Dict[str, Any]) -> str:
    even = rep["even"]
    pile = rep["pile_on"]
    fe, fp = even["fairness"], pile["fairness"]
    spec = fleet_observation_spec()
    lines = [
        "# Agent fleet CL scorecard (N≈180)",
        "",
        f"**Verdict:** {rep.get('verdict')}",
        "",
        f"Observation: `{spec['X']}` → Y=`Score`; intermediate=`{spec['intermediate']}`",
        "",
        "## Even force vs pile-on",
        "",
        "| mode | n | cohorts | Gini WIP | HHI WIP | HHI headcount | mean Score | fraud lift | CUPED ratio |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]

    def row(name: str, block: dict, fair: dict) -> str:
        t = block["tick"]
        fr = (t.get("fraud") or {}).get("lift_vs_baseline")
        cr = (t.get("cuped") or {}).get("ratio")
        return (
            f"| `{name}` | {fair['n_agents']} | {fair['n_cohorts']} | "
            f"{_fmt(fair['gini_wip'])} | {_fmt(fair['hhi_wip'])} | "
            f"{_fmt(fair['hhi_headcount'])} | {_fmt(block['mean_score'])} | "
            f"{_fmt(fr)} | {_fmt(cr)} |"
        )

    lines.append(row("even_force", even, fe))
    lines.append(row("pile_on", pile, fp))
    lines += ["", "## Default cohort mix", "", "| cohort | n |", "|---|---:|"]
    for c, n in DEFAULT_COHORTS.items():
        lines.append(f"| `{c}` | {n} |")
    lines += [
        "",
        "## Reflection (even tick)",
        "",
        f"- action: `{even['tick'].get('reflection', {}).get('action')}`",
        f"- reason: {even['tick'].get('reflection', {}).get('reason')}",
        "",
        "## Multi-tick CL reflection",
        "",
        f"- actions: `{cl.get('actions')}`",
        f"- fraud lifts: `{[round(float(x), 3) for x in cl.get('fraud_lifts') or []]}`",
        f"- CUPED ratios: `{[round(float(x), 3) for x in cl.get('cuped_ratios') or []]}`",
        f"- reading: {cl.get('reading')}",
        "",
        "## Reading",
        "",
        "- Soft-cap specialty WIP; Score only counts wired/used/roi.",
        "- Antifraud×10: optimize human useful-rate, not AUC.",
        "- CUPED×10: optimize Var ratio < 1, not another ranker.",
        "- Rehome only from `flex_reserve` → underserved cohorts.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results/agod_agent_fleet")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-agents", type=int, default=180)
    args = ap.parse_args()
    rep = compare_even_vs_pile(n_agents=args.n_agents, seed=args.seed)
    cl = multi_tick_cl_demo(n_ticks=5, seed=args.seed)
    rep["multi_tick_cl"] = {
        k: cl[k]
        for k in (
            "n_ticks",
            "actions",
            "mean_scores",
            "cuped_ratios",
            "fraud_lifts",
            "reading",
            "final_cohort_counts",
        )
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "fleet_scorecard.json").write_text(
        json.dumps(_clean(rep), indent=2) + "\n"
    )
    md = render_md(rep, cl)
    (args.out / "FLEET_SCORECARD.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
