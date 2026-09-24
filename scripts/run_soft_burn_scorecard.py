#!/usr/bin/env python3
"""Soft IPTW burn decisions + FLOPs ledger (√/∛).

  PYTHONPATH=. python3 scripts/run_soft_burn_scorecard.py \\
    --summary results/agod_po_cbrt/summary.json \\
    --out results/agod_soft_burn
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.po_power_eff import power_scorecard_from_summary
from agod.soft_burn import summarize_burn_decisions

ROOT = Path(__file__).resolve().parents[1]


def _clean(o):
    if isinstance(o, float) and (o != o or o in (float("inf"), float("-inf"))):
        return None
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    return o


def render_md(rep: dict) -> str:
    lines = [
        "# Soft IPTW burn decisions",
        "",
        rep.get("headline", ""),
        "",
        rep.get("policy", ""),
        "",
        "| dataset | duty | √ rel | ∛ rel | decision | burn? | α | reason |",
        "|---|---:|---:|---:|---|:---:|---:|---|",
    ]
    for c in rep.get("cards") or []:
        by = {m["mode"]: m for m in c.get("modes") or []}
        gs, gc = by.get("gated_sqrt", {}), by.get("gated_cbrt", {})
        b = c.get("burn") or {}

        def f(x, nd=3):
            try:
                v = float(x)
                if v != v:
                    return "—"
                return f"{v:.{nd}g}" if abs(v) < 1e-2 or abs(v) > 1e3 else f"{v:.3f}"
            except Exception:
                return "—"

        alpha = b.get("alpha")
        alpha_s = f"{alpha:.3f}" if isinstance(alpha, float) else "—"
        lines.append(
            f"| `{c.get('dataset')}` | {f(c.get('gate_duty'))} | "
            f"{f(gs.get('rel_vs_uniform'))} | {f(gc.get('rel_vs_uniform'))} | "
            f"`{b.get('decision')}` | {'Y' if b.get('burn') else 'N'} | {alpha_s} | "
            f"{b.get('reason')} |"
        )
    lines += [
        "",
        "## 算力账",
        "",
        "- Adaptation FLOPs = duty × refit（真账单）",
        "- Weighting FLOPs ≈ O(n)（α 不改总账）",
        "- 因此：该不该烧看 MSE 风险，不看「省不算力」",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--summary",
        type=Path,
        default=ROOT / "results/agod_po_cbrt/summary.json",
    )
    ap.add_argument("--out", type=Path, default=ROOT / "results/agod_soft_burn")
    args = ap.parse_args()
    power = power_scorecard_from_summary(json.loads(args.summary.read_text()))
    rep = summarize_burn_decisions(power.get("cards") or [])
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "soft_burn.json").write_text(
        json.dumps(_clean(rep), indent=2) + "\n"
    )
    md = render_md(rep)
    (args.out / "SOFT_BURN.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
