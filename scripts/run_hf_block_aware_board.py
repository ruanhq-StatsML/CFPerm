#!/usr/bin/env python3
"""Run π-guided Stage-2 on labeled Hugging Face datasets and write the LaTeX dashboard."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from block_aware_hf import run_all  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n-per-batch", type=int, default=800)
    p.add_argument("--d-text", type=int, default=32)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--n-estimators", type=int, default=80)
    p.add_argument("--n-gt-seeds", type=int, default=2)
    p.add_argument("--quick", action="store_true")
    p.add_argument(
        "--datasets",
        default="adult,imdb_rt,mnli",
        help="comma-separated: adult,imdb_rt,mnli",
    )
    args = p.parse_args()
    ds = tuple(x.strip() for x in args.datasets.split(",") if x.strip())
    board = run_all(
        n_per_batch=400 if args.quick else args.n_per_batch,
        d_text=16 if args.quick else args.d_text,
        seed=args.seed,
        n_estimators=40 if args.quick else args.n_estimators,
        n_gt_seeds=1 if args.quick else args.n_gt_seeds,
        datasets=ds,
    )
    print(json.dumps(board, indent=2))


if __name__ == "__main__":
    main()
