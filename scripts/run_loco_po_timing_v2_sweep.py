#!/usr/bin/env python3
"""LOCO PO-risk tip timing — v2 multi-mode / multi-seed sweep.

  PYTHONPATH=. python3 scripts/run_loco_po_timing_v2_sweep.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.loco_po_monitor import evaluate_score_modes, make_tip_shift_stream, run_loco_po_stream


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("results/loco_po_timing_v2"))
    ap.add_argument("--shift-at", type=int, default=20)
    ap.add_argument("--seeds", type=str, default="0,1,2,3,4")
    args = ap.parse_args()
    seeds = [int(x) for x in args.seeds.split(",") if x.strip() != ""]
    args.out.mkdir(parents=True, exist_ok=True)

    sweep = evaluate_score_modes(
        seeds=seeds,
        modes=("sum", "abs_dev", "po", "collapse", "cusum", "confirm"),
        shift_at=args.shift_at,
        true_tips=(0, 1),
        wrong_tips=(6, 7),
        n_estimators=20,
    )
    (args.out / "sweep.json").write_text(json.dumps(sweep, indent=2))

    # one detailed traj with default collapse
    stream, tips = make_tip_shift_stream(
        n_batches=40, batch_size=128, d=8, tip_idx=(0, 1), shift_at=args.shift_at, seed=0
    )
    detail = run_loco_po_stream(
        stream,
        tips,
        burn_in=8,
        seed=0,
        n_estimators=25,
        per_tip=True,
        score_mode="confirm",
        known_shift=args.shift_at,
    )
    dump = {k: v for k, v in detail.items() if k != "rows"}
    (args.out / "confirm_seed0.json").write_text(json.dumps(dump, indent=2))

    lines = [
        "# LOCO PO-risk tip timing v2 sweep",
        "",
        f"shift_at={args.shift_at}, seeds={seeds}",
        "",
        "| mode | true delay mean±std | hit±1 | miss | early/run | wrong hit±1 | wrong miss |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in sweep["rows"]:
        dm = r["true_delay_mean"]
        ds = r["true_delay_std"]
        dms = "—" if dm is None else f"{dm:.2f}±{(ds or 0):.2f}"
        lines.append(
            f"| `{r['mode']}` | {dms} | {r['true_hit_pm1']:.0%} | {r['true_miss']:.0%} | "
            f"{r['true_early_alarms_mean']:.2f} | {r['wrong_hit_pm1']:.0%} | {r['wrong_miss']:.0%} |"
        )
    best = max(
        sweep["rows"],
        key=lambda r: (
            r["true_hit_pm1"],
            r["wrong_miss"],
            -(r["true_early_alarms_mean"] or 0),
            -abs(r["true_delay_mean"] or 99),
        ),
    )
    lines += [
        "",
        f"**Preferred default:** `{best['mode']}` "
        f"(true hit±1={best['true_hit_pm1']:.0%}, wrong miss={best['wrong_miss']:.0%}, "
        f"early/run={best['true_early_alarms_mean']:.2f}).",
        "",
        f"confirm seed0 tip t*={detail['tip_first_reject_t']} delay={detail['detection_delay_tip']}.",
    ]
    (args.out / "LOCO_PO_TIMING_V2_REPORT.md").write_text("\n".join(lines) + "\n")

    # latex fragment
    tex = [
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{LOCO PO-risk tip timing v2: score-mode sweep "
        f"(shift@{args.shift_at}, {len(seeds)} seeds; true tips vs wrong tips)."+"}",
        r"\label{tab:loco-v2}",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        r"mode & true delay & hit$\pm1$ & early/run & wrong hit$\pm1$ & wrong miss \\",
        r"\midrule",
    ]
    for r in sweep["rows"]:
        dm = r["true_delay_mean"]
        ds = r["true_delay_std"]
        dms = "---" if dm is None else f"{dm:.2f}$\\pm${(ds or 0):.2f}"
        tex.append(
            f"{r['mode']} & {dms} & {100*r['true_hit_pm1']:.0f}\\% & "
            f"{r['true_early_alarms_mean']:.2f} & {100*r['wrong_hit_pm1']:.0f}\\% & "
            f"{100*r['wrong_miss']:.0f}\\% \\\\"
        )
    tex += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (args.out / "table_loco_v2.tex").write_text("\n".join(tex))
    print(json.dumps(sweep, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
