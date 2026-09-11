#!/usr/bin/env python3
"""Run the Chronoberg modality-gap clever-covariate prototype."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from chronoberg_clever_cov import run_prototype  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n-per-batch", type=int, default=1400)
    p.add_argument("--d-text", type=int, default=32)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--n-estimators", type=int, default=70)
    p.add_argument("--quick", action="store_true", help="smaller n for a smoke run")
    args = p.parse_args()
    n = 500 if args.quick else args.n_per_batch
    ns = (150, 300) if args.quick else (200, 400, 800)
    summary = run_prototype(
        n_per_batch=n,
        d_text=24 if args.quick else args.d_text,
        seed=args.seed,
        n_estimators=40 if args.quick else args.n_estimators,
        efficiency_ns=ns,
    )
    print(json.dumps({
        "n0": summary["n0"],
        "n1": summary["n1"],
        "gap": summary["gap"]["pi_consensus"],
        "comparisons": [
            {
                "name": r["name"],
                "auc_raw": r["domain_auc_raw"],
                "auc_clever": r["domain_auc_clever"],
                "delta": r["domain_auc_delta"],
                "gt": r.get("gt"),
                "mass_on_gt_raw": r.get("mass_on_gt_raw"),
                "z_share_on_gt": r.get("z_share_on_gt"),
            }
            for r in summary["comparisons"]
        ],
        "sample_efficiency": summary["sample_efficiency"],
    }, indent=2))


if __name__ == "__main__":
    main()
