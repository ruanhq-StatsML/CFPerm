#!/usr/bin/env python3
"""LOCO PO-risk tip timing monitor — synthetic smoke.

After attribution delivers tip features J*, stream batches and detect
*when* those tips start driving PO-risk via OnlineRFPerm on tip-group LOCO.

  PYTHONPATH=. python3 scripts/run_loco_po_timing_monitor.py \\
    --shift-at 20 --burn-in 8 --n-batches 40
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from agod.loco_po_monitor import make_tip_shift_stream, run_loco_po_stream


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--d", type=int, default=8)
    ap.add_argument("--tips", type=str, default="0,1", help="comma tip indices J*")
    ap.add_argument("--shift-at", type=int, default=20)
    ap.add_argument("--burn-in", type=int, default=8)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-estimators", type=int, default=25)
    ap.add_argument("--out", type=Path, default=Path("results/loco_po_timing"))
    args = ap.parse_args()

    tips = [int(x) for x in args.tips.split(",") if x.strip() != ""]
    stream, tip_idx = make_tip_shift_stream(
        n_batches=args.n_batches,
        batch_size=args.batch_size,
        d=args.d,
        tip_idx=tips,
        shift_at=args.shift_at,
        seed=args.seed,
    )
    summary = run_loco_po_stream(
        stream,
        tip_idx,
        burn_in=args.burn_in,
        alpha=args.alpha,
        seed=args.seed,
        n_estimators=args.n_estimators,
        per_tip=True,
        known_shift=args.shift_at,
    )

    args.out.mkdir(parents=True, exist_ok=True)
    # drop heavy row blobs for json top-level compactness
    dump = {k: v for k, v in summary.items() if k != "rows"}
    dump["n_rows"] = len(summary["rows"])
    (args.out / "summary.json").write_text(json.dumps(dump, indent=2))

    # light traj csv
    lines = [
        "t,tip_group_loco,tip_change,score,po_full,tip_share,hard_reject,tip_reject,po_reject"
    ]
    for i, row in enumerate(summary["rows"]):
        lines.append(
            f"{row['t']},{row['tip_group_loco']:.6f},{row['tip_change']:.6f},"
            f"{row['score']:.6f},{row['po_full']:.6f},{row['tip_share']:.6f},"
            f"{int(summary['hard_reject_hist'][i])},"
            f"{int(summary['tip_reject_hist'][i])},"
            f"{int(summary['po_reject_hist'][i])}"
        )
    (args.out / "traj.csv").write_text("\n".join(lines) + "\n")

    md = [
        "# LOCO PO-risk tip timing monitor",
        "",
        f"- tips J* = `{tip_idx}`",
        f"- known shift @ batch **{args.shift_at}**",
        f"- burn-in = {args.burn_in}, α = {args.alpha}, "
        f"L̄ = {summary.get('loco_center')}, thr = {summary.get('score_threshold')}",
        f"- **hard-gate tip change-point t** = **{summary['tip_first_reject_t']}** "
        f"(delay={summary['detection_delay_tip']})",
        f"- tip OnlineRFPerm first reject t = {summary.get('tip_rfperm_first_reject_t')}",
        f"- full-PO first reject t = **{summary['po_first_reject_t']}** "
        f"(delay={summary['detection_delay_po']})",
        f"- lead tip−PO = `{summary['lead_tip_vs_po']}` (negative ⇒ tip earlier)",
        "",
        "Score: `S_t = |L_tip(t) − L̄| + PO_t`; timing = first "
        "`S > mean(S_burn)+k·std(S_burn)`.",
        "",
        f"Artifacts: `{args.out}/summary.json`, `{args.out}/traj.csv`.",
    ]
    (args.out / "LOCO_PO_TIMING_REPORT.md").write_text("\n".join(md) + "\n")

    print(json.dumps(dump, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
