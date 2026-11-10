#!/usr/bin/env python3
"""FSDS × Interaction Value Online Adaptation board.

  python3 scripts/run_fsds_interaction_online_adapt.py
  python3 scripts/run_fsds_interaction_online_adapt.py --quick
  python3 scripts/run_fsds_interaction_online_adapt.py --scenarios compound,interaction_break
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fsds_interaction_online_adapt import (  # noqa: E402
    METHODS,
    SCENARIOS,
    plot_suite,
    run_msrvtt_hop_probe,
    run_suite,
    write_readme,
)
from msrvtt_multimodal_attribution import write_json  # noqa: E402

OUT = ROOT / "results" / "fsds_interaction_online_adapt"
DOCS = ROOT / "docs" / "method"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-per", type=int, default=120)
    ap.add_argument("--n-batches", type=int, default=8)
    ap.add_argument(
        "--scenarios",
        type=str,
        default="",
        help="comma list; default all SCENARIOS",
    )
    ap.add_argument("--skip-msrvtt", action="store_true")
    args = ap.parse_args()

    scenarios = (
        [s.strip() for s in args.scenarios.split(",") if s.strip()]
        if args.scenarios
        else list(SCENARIOS)
    )
    n_seeds = 2 if args.quick else int(args.seeds)
    seeds = list(range(2026, 2026 + n_seeds))
    print(
        "fsds-interaction-online",
        "quick" if args.quick else "full",
        n_seeds,
        "seeds",
        args.n_per,
        "per",
        flush=True,
    )
    suite = run_suite(
        scenarios=scenarios,
        methods=list(METHODS),
        seeds=seeds,
        n_batches=int(args.n_batches),
        n_per=int(args.n_per),
        quick=bool(args.quick),
    )
    probe = None
    if not args.skip_msrvtt:
        print("msrvtt hop probe…", flush=True)
        # load from repo root so feature zip resolves
        import os

        os.chdir(ROOT)
        probe = run_msrvtt_hop_probe(n_videos=16 if args.quick else 32, n_windows=16)

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {
        "summary": suite["summary"],
        "lifts": suite["lifts"],
        "scenarios": suite["scenarios"],
        "methods": suite["methods"],
        "seeds": suite["seeds"],
        "n_batches": suite["n_batches"],
        "n_per": suite["n_per"],
        "msrvtt_probe": probe,
        # histories are large; keep aggregate rows without nested history for json
        "rows": [
            {k: v for k, v in r.items() if k != "history"} for r in suite["rows"]
        ],
    }
    write_json(OUT / "fsds_interaction_online_adapt.json", payload)
    plot_suite(suite, OUT / "fsds_interaction_online_adapt.png")
    write_readme(suite, probe, OUT / "README.md")
    write_readme(suite, probe, DOCS / "FSDS_Interaction_Online_Adapt.md")

    for sc in suite["scenarios"]:
        print("scenario", sc, flush=True)
        for method in suite["methods"]:
            cell = suite["summary"][sc][method]["online_mse"]
            lift = suite["lifts"][sc][method]
            print(
                "  %-16s mse=%.4f (%.4f)  vs_static=%+.4f  vs_loc=%+.4f"
                % (
                    method,
                    cell["mean"],
                    cell["sd"],
                    lift["vs_static"],
                    lift["vs_localize"],
                ),
                flush=True,
            )
    if probe:
        print(
            "msrvtt",
            probe["source"],
            "phi_TI=%.4f" % probe["mean_phi_TI"],
            "tri=%.3f" % probe["mean_triangle"],
            "pi=",
            json.dumps(probe["pi"]),
            flush=True,
        )


if __name__ == "__main__":
    main()
