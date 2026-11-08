#!/usr/bin/env python3
"""FICO weak-5 + nonlinear concept-drift CFPerm vs median-split baseline.

  python3 scripts/run_fico_weak_cd_cfperm.py
  python3 scripts/run_fico_weak_cd_cfperm.py --quick
  python3 scripts/run_fico_weak_cd_cfperm.py --n-draws 5 --n-perm 200 --n-estimators 150
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fico_heloc_weak_cd import run_weak_cd_experiment, write_weak_cd_outputs  # noqa: E402

OUT = ROOT / "results" / "fico_heloc_cfperm" / "weak_cd"
CACHE = ROOT / "data" / "heloc_dataset_v1.csv"
BASELINE = ROOT / "results" / "fico_heloc_cfperm" / "power_summary.csv"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--n-draws", type=int, default=5)
    ap.add_argument("--k-drivers", type=int, default=5)
    ap.add_argument("--intensity", type=float, default=2.5)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-estimators", type=int, default=150)
    ap.add_argument("--n-perm-cate", type=int, default=50)
    ap.add_argument("--max-n-per-arm", type=int, default=1500)
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--vimps", type=str, default="loco,permucate,grf")
    ap.add_argument("--no-posthoc", action="store_true")
    args = ap.parse_args()

    if args.quick:
        n_draws, n_perm, n_estimators, n_perm_cate, max_n = 2, 8, 40, 8, 250
    else:
        n_draws = int(args.n_draws)
        n_perm = int(args.n_perm)
        n_estimators = int(args.n_estimators)
        n_perm_cate = int(args.n_perm_cate)
        max_n = int(args.max_n_per_arm)

    vimps = tuple(v.strip().lower() for v in args.vimps.split(",") if v.strip())
    print(
        f"weak-cd draws={n_draws} k={args.k_drivers} intensity={args.intensity} "
        f"vimps={vimps} B={n_estimators} n_perm={n_perm}",
        flush=True,
    )
    payload = run_weak_cd_experiment(
        cache_path=CACHE,
        n_draws=n_draws,
        k_drivers=int(args.k_drivers),
        intensity=float(args.intensity),
        vimps=vimps,
        n_perm=n_perm,
        n_estimators=n_estimators,
        n_perm_cate=n_perm_cate,
        max_n_per_arm=max_n,
        seed=int(args.seed),
        n_jobs=int(args.n_jobs),
        run_posthoc=not args.no_posthoc,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    write_weak_cd_outputs(payload, OUT, baseline_power_csv=BASELINE)
    print("power:", payload["power"], flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
