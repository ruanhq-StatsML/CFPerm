#!/usr/bin/env python3
"""Confirm + tip-cmean joint formulation prototype smoke.

  PYTHONPATH=. python3 scripts/run_cmean_confirm_prototype.py --shift-at 20
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agod.cmean_confirm import run_cmean_confirm_prototype
from agod.loco_po_monitor import make_tip_shift_stream


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shift-at", type=int, default=20)
    ap.add_argument("--burn-in", type=int, default=8)
    ap.add_argument("--hard-k", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/cmean_confirm_proto"))
    args = ap.parse_args()

    stream, tips = make_tip_shift_stream(
        n_batches=40,
        batch_size=128,
        d=8,
        tip_idx=(0, 1),
        shift_at=args.shift_at,
        seed=args.seed,
        shift_mean=2.5,
    )
    true = run_cmean_confirm_prototype(
        stream,
        tips,
        burn_in=args.burn_in,
        hard_k=args.hard_k,
        seed=args.seed,
        n_estimators=20,
        known_shift=args.shift_at,
    )
    wrong = run_cmean_confirm_prototype(
        stream,
        (6, 7),
        burn_in=args.burn_in,
        hard_k=args.hard_k,
        seed=args.seed,
        n_estimators=20,
        known_shift=args.shift_at,
    )

    args.out.mkdir(parents=True, exist_ok=True)
    dump = {
        "true_tips": {k: v for k, v in true.items() if k != "rows"},
        "wrong_tips": {k: v for k, v in wrong.items() if k != "rows"},
    }
    (args.out / "summary.json").write_text(json.dumps(dump, indent=2))

    thr = true["thresholds"]
    md = [
        "# Confirm + tip-cmean formulation (prototype)",
        "",
        "## Joint fire",
        "`Fire = 1{Δ_tip ≥ ε_Δ} · 1{confirm}` with",
        "`confirm = 1{collapse ≥ thr_col} · 1{PO ≥ thr_PO}`.",
        "",
        "## How magnitude thresholds are set",
        f"- recipe: **`{thr['recipe']}`** with `k={args.hard_k}`",
        f"- `ε_Δ = mean(Δ_tip_burn) + k·sd = {thr['delta_burn_mean']:.4f} "
        f"+ {args.hard_k}·{thr['delta_burn_std']:.4f} = **{thr['eps_delta']:.4f}**`",
        f"- `thr_collapse = {thr['collapse_burn_mean']:.4f} + k·"
        f"{thr['collapse_burn_std']:.4f} = **{thr['thr_collapse']:.4f}**`",
        f"- `thr_PO = {thr['po_burn_mean']:.4f} + k·{thr['po_burn_std']:.4f} "
        f"= **{thr['thr_po']:.4f}**`",
        "",
        "Burn-only calibration; engineering gates (not Type-I α).",
        "",
        "## Result",
        f"- true tips {tips}: t*={true['t_star']} delay={true['detection_delay']}",
        f"- wrong tips (6,7): t*={wrong['t_star']} delay={wrong['detection_delay']}",
        f"- signed tip shifts at fire: `{true['signed_tips_at_fire']}`",
    ]
    (args.out / "CMEAN_CONFIRM_PROTO.md").write_text("\n".join(md) + "\n")
    print(json.dumps(dump, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
