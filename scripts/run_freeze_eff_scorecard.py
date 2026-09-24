#!/usr/bin/env python3
"""Freeze closed-loop MSE–FLOPs scorecard (Grad-OnlineRFPerm).

Uses canonical numbers from Grad_OnlineRFPerm_extras by default, or a
freeze bundle JSON if provided.

  PYTHONPATH=. python3 scripts/run_freeze_eff_scorecard.py \\
    --out results/agod_freeze_eff
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.freeze_eff import scorecard_canonical, scorecard_from_freeze_bundle

ROOT = Path(__file__).resolve().parents[1]


def render_md(card: dict) -> str:
    lines = [
        "# Grad-RFPerm freeze: MSE–FLOPs efficiency",
        "",
        card.get("headline", ""),
        "",
        f"- source=`{card.get('source')}` · datasets={card.get('n_datasets')} · "
        f"Pareto wins={card.get('n_pareto_dominances')}",
        "",
        "| dataset | policy | mse_rel | flops_rel | freeze_eff | Pareto? |",
        "|---|---|---:|---:|---:|:---:|",
    ]
    for c in card.get("cards") or []:
        for r in c.get("policies") or []:
            if r["policy"] == "no_adapt":
                continue
            pe = r.get("freeze_eff")
            pe_s = f"{pe:.3f}" if pe == pe else "—"
            lines.append(
                f"| `{c['dataset']}` | `{r['policy']}` | "
                f"{r['mse_rel']:.2f} | {r['flops_rel']:.2f} | {pe_s} | "
                f"{'Y' if r['dominates_always'] else 'N'} |"
            )
        lines.append(f"| | | | | | *{c.get('reading')}* |")
    lines += [
        "",
        "## Justify",
        "",
        "- After Grad reject, freezing low-share / early layers cuts adapt FLOPs.",
        "- **Pareto** = MSE better *and* FLOPs lower than always_adapt.",
        "- electricity: freeze policies ≈0.86× MSE @ ~0.7× FLOPs (Pareto).",
        "- synthetic: FLOPs save but MSE↑ — freeze_eff can be negative; still "
        "prefer freeze_low_share over freeze_early.",
        "- `no_adapt` collapses — do not equate freeze with stop.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", type=Path, default=None, help="optional freeze JSON")
    ap.add_argument("--out", type=Path, default=ROOT / "results/agod_freeze_eff")
    args = ap.parse_args()
    if args.bundle and args.bundle.exists():
        raw = json.loads(args.bundle.read_text())
        freeze = raw.get("freeze", raw)
        card = scorecard_from_freeze_bundle(freeze)
    else:
        card = scorecard_canonical()
    args.out.mkdir(parents=True, exist_ok=True)

    def _clean(o):
        if isinstance(o, float) and (o != o or o in (float("inf"), float("-inf"))):
            return None
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_clean(v) for v in o]
        return o

    (args.out / "freeze_eff.json").write_text(
        json.dumps(_clean(card), indent=2) + "\n"
    )
    md = render_md(card)
    (args.out / "FREEZE_EFF.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
