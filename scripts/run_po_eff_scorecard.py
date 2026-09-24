#!/usr/bin/env python3
"""PO-risk efficiency scorecard: ranking/MSE per adaptation FLOP.

Reads an existing ``run_agod_po_ref_vs_refit`` summary (no re-train) and
reports which of ref / probe / refit buys skill per relative FLOPs.

  PYTHONPATH=. python3 scripts/run_po_eff_scorecard.py \\
    --summary results/agod_po_ref_vs_refit/summary.json \\
    --out results/agod_po_eff
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.po_eff import scorecard_from_summary

ROOT = Path(__file__).resolve().parents[1]


def render_md(card: dict) -> str:
    lines = [
        "# PO-risk adaptation efficiency",
        "",
        card.get("headline", ""),
        "",
        f"- batches={card['n_batches']} · batch_size={card['batch_size']} · "
        f"n_control={card['n_control']} · datasets={card['n_datasets']}",
        "",
        "## Per dataset",
        "",
        "| dataset | n_rej | ref ρ | probe Δρ | refit Δρ | "
        "probe rank_eff | refit rank_eff | probe mse_eff | refit mse_eff | reading |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for c in card.get("cards") or []:
        if not c.get("ok"):
            lines.append(f"| {c.get('dataset')} | — | | | | | | | | {c.get('reason')} |")
            continue
        by = {m["mode"]: m for m in c["modes"]}
        p, r = by.get("probe", {}), by.get("refit", {})

        def f(x, nd=4):
            try:
                v = float(x)
                return f"{v:.{nd}g}" if abs(v) < 1e-2 or abs(v) > 1e3 else f"{v:.4f}"
            except Exception:
                return "—"

        lines.append(
            f"| `{c['dataset']}` | {c['n_reject']} | {f(c['ref_spearman'], 3)} | "
            f"{f(p.get('delta_spearman_vs_ref'), 3)} | {f(r.get('delta_spearman_vs_ref'), 3)} | "
            f"{f(p.get('rank_eff'))} | {f(r.get('rank_eff'))} | "
            f"{f(p.get('mse_eff'))} | {f(r.get('mse_eff'))} | {c.get('reading')} |"
        )
    lines += [
        "",
        "## How to read",
        "",
        "- **rank_eff**: ΔSpearman vs frozen ref per million adaptation FLOPs.",
        "- **mse_eff**: (uniform − mode) sig-only MSE per million FLOPs "
        "(positive = cheaper error drop).",
        "- Probe pays *always-on* fits; refit pays only on reject — so a tiny "
        "Δρ for refit can still beat probe on rank_eff.",
        "- Negative mse_eff with positive rank_eff ⇒ keep for triage, not for IPTW.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--summary",
        type=Path,
        default=ROOT / "results/agod_po_ref_vs_refit/summary.json",
    )
    ap.add_argument(
        "--out", type=Path, default=ROOT / "results/agod_po_eff"
    )
    args = ap.parse_args()
    summary = json.loads(args.summary.read_text())
    card = scorecard_from_summary(summary)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "po_eff_scorecard.json").write_text(
        json.dumps(card, indent=2) + "\n"
    )
    md = render_md(card)
    (args.out / "PO_EFF_SCORECARD.md").write_text(md)
    print(md)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
