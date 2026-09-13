#!/usr/bin/env python3
"""Attribution → next-batch training adapter on Amazon rating MSE.

  python3 scripts/run_attribution_adapter.py
  python3 scripts/run_attribution_adapter.py --quick
  python3 scripts/run_attribution_adapter.py --synthetic
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from attribution_adapter import (  # noqa: E402
    METHODS,
    plot_adapter_suite,
    run_adapter_suite,
    write_adapter_tex,
)
from msrvtt_multimodal_attribution import write_json  # noqa: E402

OUT = ROOT / "results" / "attribution_adapter"
DOCS = ROOT / "docs" / "method"


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
    print("attribution adapter", source, n_seeds, "seeds", n_per, "per", flush=True)
    suite = run_adapter_suite(
        seeds=seeds,
        n_per=n_per,
        n_batches=9,
        source=source,
        steps_per_batch=4 if args.quick else 8,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {
        "table": suite["table"],
        "lifts": suite["lifts"],
        "shapes": suite["shapes"],
        "seeds": suite["seeds"],
        "source": suite["source"],
        "methods": suite["methods"],
        "rows": suite["rows"],
        "hop_logs": suite["hop_logs"],
    }
    write_json(OUT / "attribution_adapter.json", payload)
    plot_adapter_suite(suite, OUT / "attribution_adapter.png")
    write_adapter_tex(suite, OUT / "Attribution_adapter.tex")
    write_adapter_tex(suite, DOCS / "Attribution_adapter.tex")
    live = suite["shapes"].get("amazon_live") or suite["shapes"].get("amazon") or {}
    print("shapes", json.dumps(live, indent=2), flush=True)
    for method in METHODS:
        cell = suite["table"][method]["online_mse"]
        print(
            "  %s  online_mse=%.3f (%.3f)" % (method, cell["mean"], cell["sd"]),
            flush=True,
        )
    if "bank" in suite["lifts"]:
        print("lift vs bank:", flush=True)
        for m, cell in suite["lifts"]["bank"].items():
            print("  %s  Δ=%.3f" % (m, cell["mean_mse_drop"]), flush=True)


if __name__ == "__main__":
    main()
