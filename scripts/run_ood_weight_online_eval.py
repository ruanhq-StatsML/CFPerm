#!/usr/bin/env python3
"""Run OOD-risk-as-weights online eval (MSE / regret / BWT) across boards.

  python3 scripts/run_ood_weight_online_eval.py
  python3 scripts/run_ood_weight_online_eval.py --quick
  python3 scripts/run_ood_weight_online_eval.py --boards amazon,diabetes,interstate,stock
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from ood_weight_online_eval import (  # noqa: E402
    OUT,
    WEIGHT_METHODS,
    plot_ood_eval,
    plot_paths,
    run_multi_board,
    write_json,
    write_readme,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument(
        "--boards",
        type=str,
        default="amazon,synthetic,diabetes,interstate,stock",
        help="comma-separated boards",
    )
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-per", type=int, default=240)
    ap.add_argument("--n-batches", type=int, default=9)
    ap.add_argument(
        "--methods",
        type=str,
        default=",".join(WEIGHT_METHODS),
        help="comma-separated weight methods",
    )
    args = ap.parse_args()
    boards = [b.strip() for b in args.boards.split(",") if b.strip()]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    seeds = list(range(2026, 2026 + (2 if args.quick else int(args.seeds))))
    print(
        "ood-weight online eval",
        boards,
        methods,
        "seeds",
        seeds,
        "n_per",
        args.n_per,
        flush=True,
    )
    multi = run_multi_board(
        boards=boards,
        seeds=seeds,
        methods=methods,
        n_per=int(args.n_per),
        n_batches=int(args.n_batches),
        quick=bool(args.quick),
    )
    OUT.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "ood_weight_online_eval.json", multi)
    plot_ood_eval(multi, OUT / "ood_weight_online_eval.png")
    for board in boards:
        plot_paths(multi, OUT / ("path_%s.png" % board), board=board)
    write_readme(multi, OUT / "README.md")
    print("wrote", OUT, flush=True)
    for board in boards:
        print("==", board, "==", flush=True)
        for m in methods:
            t = multi["suites"][board]["table"][m]
            print(
                "  %-8s  mse=%.4f  regret=%+.4f  bwt=%.4f"
                % (
                    m,
                    t["online_mse"]["mean"],
                    t["regret"]["mean"],
                    t["bwt"]["mean"],
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
