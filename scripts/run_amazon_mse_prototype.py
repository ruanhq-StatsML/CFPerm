#!/usr/bin/env python3
"""Amazon MSE prototype: SDC star-prototypes vs GPM-Ridge vs plateau.

  python3 scripts/run_amazon_mse_prototype.py
  python3 scripts/run_amazon_mse_prototype.py --synthetic
  python3 scripts/run_amazon_mse_prototype.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from amazon_mse_prototype import METHODS, plot_mse_suite, run_mse_suite
from msrvtt_multimodal_attribution import write_json

OUT = ROOT / "results" / "amazon_mse_prototype"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-per", type=int, default=240)
    args = ap.parse_args()
    n_per = 80 if args.quick else int(args.n_per)
    n_seeds = 2 if args.quick else int(args.seeds)
    source = "synthetic" if args.synthetic else "amazon"
    seeds = list(range(2026, 2026 + n_seeds))
    print("amazon mse prototype", source, n_seeds, "seeds", n_per, "per", flush=True)
    suite = run_mse_suite(
        seeds=seeds,
        n_per=n_per,
        n_batches=9,
        source=source,
        steps_per_batch=4 if args.quick else 8,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    write_json(
        OUT / "amazon_mse_prototype.json",
        {
            "table": suite["table"],
            "shapes": suite["shapes"],
            "seeds": suite["seeds"],
            "source": suite["source"],
            "rows": suite["rows"],
        },
    )
    plot_mse_suite(suite, OUT / "amazon_mse_prototype.png")
    print("shapes", json.dumps(suite["shapes"]["amazon_live"], indent=2), flush=True)
    for method in METHODS:
        cell = suite["table"][method]
        print(
            "  %s  online_mse=%.3f (%.3f)"
            % (method, cell["online_mse"]["mean"], cell["online_mse"]["sd"]),
            flush=True,
        )


if __name__ == "__main__":
    main()
