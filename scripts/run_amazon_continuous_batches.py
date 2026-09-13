#!/usr/bin/env python3
"""Amazon review continuous batches: TSS vs LR clocks on rating MSE.

  python3 scripts/run_amazon_continuous_batches.py
  python3 scripts/run_amazon_continuous_batches.py --synthetic
  python3 scripts/run_amazon_continuous_batches.py --quick
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from amazon_continuous_batches import (  # noqa: E402
    CATEGORIES,
    METHODS,
    plot_amazon_suite,
    plot_relationship_heatmap,
    run_amazon_suite,
    strip_traces,
    write_tex_table,
)
from msrvtt_multimodal_attribution import write_json  # noqa: E402

OUT = ROOT / "results" / "amazon_continuous_batches"
DOCS = ROOT / "docs" / "method"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--n-per", type=int, default=240)
    ap.add_argument("--n-batches", type=int, default=9)
    args = ap.parse_args()

    n_batches = min(int(args.n_batches), len(CATEGORIES))
    n_per = 80 if args.quick else int(args.n_per)
    n_seeds = 3 if args.quick else int(args.seeds)
    source = "synthetic" if args.synthetic else "amazon"
    seeds = list(range(2026, 2026 + n_seeds))
    methods = list(METHODS)

    print("amazon suite", source, n_seeds, "seeds", n_batches, "batches", n_per, "per", flush=True)
    suite = run_amazon_suite(
        seeds=seeds,
        n_batches=n_batches,
        n_per=n_per,
        methods=methods,
        source=source,
        steps_per_batch=4 if args.quick else 8,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "amazon_vs_schedulers.json", strip_traces(suite))
    plot_amazon_suite(suite, OUT / "amazon_vs_schedulers.png")
    if suite.get("relationship"):
        plot_relationship_heatmap(suite["relationship"], OUT / "amazon_batch_relationship_heatmap.png")
        write_json(
            OUT / "amazon_batch_relationship.json",
            {
                k: v
                for k, v in suite["relationship"].items()
            },
        )
    write_tex_table(suite, OUT / "Amazon_continuous_batches.tex")
    write_tex_table(suite, DOCS / "Amazon_continuous_batches.tex")

    print("hindsight eta", [r["hindsight_eta"] for r in suite["hindsight"]], flush=True)
    print("identification", suite["identification"][0], flush=True)
    for method, cell in suite["table"].items():
            print(
                "  %s  mse=%.3f  regret=%.3f  bwt=%.3f  eta=%.4f  c=%.3f  d=%.3f"
                % (
                    method,
                    cell["online_mse"]["mean"],
                    cell["regret"]["mean"],
                    cell["bwt"]["mean"],
                    cell["mean_eta"]["mean"],
                    cell["mean_c"]["mean"],
                    cell["mean_delta"]["mean"],
                ),
                flush=True,
            )
    ident = suite["identification"]
    if ident:
        mean_c = sum(r["mean_c"] for r in ident) / len(ident)
        mean_d = sum(r["mean_delta"] for r in ident) / len(ident)
        print("mean_c %.3f  mean_delta %.3f" % (mean_c, mean_d), flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
