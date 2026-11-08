#!/usr/bin/env python3
"""Run FICO HELOC CFPerm (LOCO / permucate / grf) + subset-localization power.

  python3 scripts/run_fico_heloc_cfperm.py
  python3 scripts/run_fico_heloc_cfperm.py --quick
  python3 scripts/run_fico_heloc_cfperm.py --n-perm 200 --n-estimators 150
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fico_heloc_cfperm import run_fico_experiment, write_outputs  # noqa: E402

OUT = ROOT / "results" / "fico_heloc_cfperm"
CACHE = ROOT / "data" / "heloc_dataset_v1.csv"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="smoke settings")
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-estimators", type=int, default=150, help="B / num.trees")
    ap.add_argument("--n-perm-cate", type=int, default=50)
    ap.add_argument("--max-n-per-arm", type=int, default=1500)
    ap.add_argument("--n-jobs", type=int, default=2)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument(
        "--vimps",
        type=str,
        default="loco,permucate,grf",
        help="comma-separated VIMP notions",
    )
    args = ap.parse_args()

    n_perm = 8 if args.quick else int(args.n_perm)
    n_estimators = 40 if args.quick else int(args.n_estimators)
    n_perm_cate = 8 if args.quick else int(args.n_perm_cate)
    max_n = 250 if args.quick else int(args.max_n_per_arm)
    vimps = tuple(v.strip().lower() for v in args.vimps.split(",") if v.strip())

    print(
        f"fico-cfperm vimps={vimps} B={n_estimators} n_perm={n_perm} "
        f"n_perm_cate={n_perm_cate} max_n_per_arm={max_n}",
        flush=True,
    )
    payload = run_fico_experiment(
        cache_path=CACHE,
        vimps=vimps,
        n_perm=n_perm,
        n_estimators=n_estimators,
        n_perm_cate=n_perm_cate,
        max_n_per_arm=max_n,
        seed=int(args.seed),
        n_jobs=int(args.n_jobs),
    )
    OUT.mkdir(parents=True, exist_ok=True)
    write_outputs(payload, OUT)
    power = payload["power"]
    print("power overall:", power["overall"], flush=True)
    for v, d in power["by_vimp"].items():
        print(f"  {v}:", d, flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
