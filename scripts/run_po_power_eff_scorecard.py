#!/usr/bin/env python3
"""Sqrt vs cbrt IPTW efficiency (same FLOPs, softer weights).

  PYTHONPATH=. python3 scripts/run_po_power_eff_scorecard.py \\
    --summary results/agod_po_cbrt/summary.json \\
    --out results/agod_po_power_eff
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.po_power_eff import power_scorecard_from_summary

ROOT = Path(__file__).resolve().parents[1]


def render_md(card: dict) -> str:
    lines = [
        "# IPTW power softness: √ vs ∛ (same FLOPs)",
        "",
        card.get("headline", ""),
        "",
        card.get("note", ""),
        "",
        f"- batches={card['n_batches']} · batch_size={card['batch_size']} · "
        f"datasets={card['n_datasets']} · "
        f"∛≤√ rate={card.get('soft_win_rate_cbrt_le_sqrt')}",
        "",
        "| dataset | duty | unif | gated√ rel | gated∛ rel | "
        "gated√ mse_eff | gated∛ mse_eff | ∛≤√? | best | burn | reading |",
        "|---|---:|---:|---:|---:|---:|---:|:---:|---|---|---|",
    ]
    from agod.soft_burn import attach_burn_to_power_card

    for c in card.get("cards") or []:
        c = attach_burn_to_power_card(c)
        by = {m["mode"]: m for m in c.get("modes") or []}
        u, gs, gc = by.get("uniform", {}), by.get("gated_sqrt", {}), by.get("gated_cbrt", {})
        b = c.get("burn") or {}

        def f(x, nd=4):
            try:
                v = float(x)
                if not (v == v):
                    return "—"
                return f"{v:.{nd}g}" if abs(v) < 1e-2 or abs(v) > 1e3 else f"{v:.4f}"
            except Exception:
                return "—"

        sw = c.get("soft_win_cbrt_le_sqrt")
        sw_s = "Y" if sw is True else ("N" if sw is False else "—")
        lines.append(
            f"| `{c['dataset']}` | {f(c.get('gate_duty'), 3)} | "
            f"{f(u.get('mse_mean_sig'))} | {f(gs.get('rel_vs_uniform'), 3)} | "
            f"{f(gc.get('rel_vs_uniform'), 3)} | {f(gs.get('mse_eff'))} | "
            f"{f(gc.get('mse_eff'))} | {sw_s} | `{c.get('best_sig_mode')}` | "
            f"`{b.get('decision')}` | {c.get('reading')} |"
        )
    lines += [
        "",
        "## Burn policy (α ≠ FLOPs)",
        "",
        "- See [`Soft_Weight_Burn_Logic.md`](../../docs/summaries/Soft_Weight_Burn_Logic.md).",
        "- Default: do not burn; `SOFTEN_ONLY` → always ∛ if gate is mandatory.",
        "- Burn only when gated_α beats uniform on sig-MSE.",
        "",
        "## Tomorrow-demo takeaway",
        "",
        "1. Softness ≠ free lunch: uniform still wins most packs.",
        "2. When you *do* gate IPTW, ∛ ≤ √ on a majority of packs here.",
        "3. FLOPs identical → choose α by MSE risk, not by compute.",
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
    ap.add_argument(
        "--out", type=Path, default=ROOT / "results/agod_po_power_eff"
    )
    args = ap.parse_args()
    summary = json.loads(args.summary.read_text())
    card = power_scorecard_from_summary(summary)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "po_power_eff.json").write_text(json.dumps(card, indent=2) + "\n")
    md = render_md(card)
    (args.out / "PO_POWER_EFF.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
