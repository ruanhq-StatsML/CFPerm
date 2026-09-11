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
    p.add_argument("--n-per-batch", type=int, default=1000)
    p.add_argument("--d-text", type=int, default=48)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--n-estimators", type=int, default=100)
    p.add_argument("--quick", action="store_true", help="smaller n for a smoke run")
    args = p.parse_args()
    n = 400 if args.quick else args.n_per_batch
    ns = (160, 320) if args.quick else (200, 400, 800, 1200)
    summary = run_prototype(
        n_per_batch=n,
        d_text=24 if args.quick else args.d_text,
        seed=args.seed,
        n_estimators=40 if args.quick else args.n_estimators,
        efficiency_ns=ns,
        n_repeats=2 if args.quick else 4,
        n_gt_seeds=2 if args.quick else 3,
    )
    print(json.dumps({
        "n0": summary["n0"],
        "n1": summary["n1"],
        "gap": summary["gap"]["pi_consensus"],
        "comparisons": [
            {
                "name": r["name"],
                "auc_raw": r["domain_auc_raw"],
                "auc_pool": r.get("domain_auc_pool"),
                "auc_subspace": r.get("domain_auc_subspace"),
                "auc_bawf": r.get("domain_auc_bawf"),
                "auc_adapt": r.get("domain_auc_adapt"),
                "auc_logit_unif": r.get("domain_auc_logit_unif"),
                "auc_stack": r.get("domain_auc_stack_xz"),
                "delta_bawf_vs_sub": round(
                    r.get("domain_auc_bawf", r["domain_auc_raw"])
                    - r.get("domain_auc_subspace", r["domain_auc_raw"]),
                    4,
                ),
                "gt": r.get("gt"),
                "mass_on_gt_raw": r.get("mass_on_gt_raw"),
                "pi_on_gt": r.get("pi_consensus_on_gt"),
                "bawf_on_gt": r.get("bawf_on_gt"),
            }
            for r in summary["comparisons"]
        ],
        "sample_efficiency": summary["sample_efficiency"],
    }, indent=2))


if __name__ == "__main__":
    main()
