#!/usr/bin/env python3
"""FSDS × Interaction Value board — primary eval is Φ(i,j) visualization.

  python3 scripts/run_fsds_interaction_online_adapt.py
  python3 scripts/run_fsds_interaction_online_adapt.py --quick
  python3 scripts/run_fsds_interaction_online_adapt.py --viz-only
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
    collect_interaction_viz,
    plot_interaction_diag_profiles,
    plot_interaction_viz,
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
    ap.add_argument("--viz-only", action="store_true", help="skip MSE suite; only Φ(i,j) board")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-per", type=int, default=120)
    ap.add_argument("--n-batches", type=int, default=8)
    ap.add_argument("--scenarios", type=str, default="")
    ap.add_argument("--skip-msrvtt", action="store_true")
    args = ap.parse_args()

    scenarios = (
        [s.strip() for s in args.scenarios.split(",") if s.strip()]
        if args.scenarios
        else list(SCENARIOS)
    )
    viz_scenarios = [s for s in scenarios if s != "concept_drift"] or scenarios
    if args.quick:
        viz_scenarios = [s for s in ("healthy", "interaction_break", "compound") if s in scenarios or not args.scenarios]
        if not viz_scenarios:
            viz_scenarios = ["healthy", "interaction_break", "compound"]

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    print("collect Φ(i,j) interaction viz…", flush=True)
    viz = collect_interaction_viz(
        scenarios=viz_scenarios,
        seed=2026,
        n_batches=6 if args.quick else 8,
        n_per=60 if args.quick else 100,
    )
    write_json(OUT / "interaction_viz.json", viz)
    plot_interaction_viz(viz, OUT / "interaction_TI_matrix.png", pair="TI")
    plot_interaction_viz(viz, OUT / "interaction_TA_matrix.png", pair="TA")
    plot_interaction_diag_profiles(viz, OUT / "interaction_diag_profiles.png")
    for sc in viz["scenarios"]:
        ti = viz["board"][sc]["pairs"]["TI"]
        print(
            "Φ[%s] TI diag_drop=%.4f  ‖Δ‖F=%.3f  topΔ=(%d,%d)%+.4f"
            % (
                sc,
                ti["diag_drop"],
                ti["frobenius_delta"],
                ti["top_delta_abs"][0]["i"],
                ti["top_delta_abs"][0]["j"],
                ti["top_delta_abs"][0]["phi"],
            ),
            flush=True,
        )

    suite = {
        "summary": {},
        "lifts": {},
        "interpret_board": {},
        "scenarios": viz_scenarios,
        "methods": [],
        "seeds": [],
        "n_batches": viz["n_batches"],
        "n_per": viz["n_per"],
    }
    if not args.viz_only:
        n_seeds = 2 if args.quick else int(args.seeds)
        seeds = list(range(2026, 2026 + n_seeds))
        print(
            "secondary MSE suite",
            "quick" if args.quick else "full",
            n_seeds,
            "seeds",
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
        plot_suite(suite, OUT / "fsds_interaction_online_adapt.png")

    probe = None
    if not args.skip_msrvtt:
        print("msrvtt hop probe…", flush=True)
        import os

        os.chdir(ROOT)
        probe = run_msrvtt_hop_probe(n_videos=16 if args.quick else 32, n_windows=16)

    payload = {
        "eval": "interaction(i,j) visualization is primary; MSE is footnote",
        "interaction_viz": {
            "scenarios": viz["scenarios"],
            "shared_dim": viz["shared_dim"],
            "summary": {
                sc: {
                    pair: {
                        "diag_drop": viz["board"][sc]["pairs"][pair]["diag_drop"],
                        "frobenius_delta": viz["board"][sc]["pairs"][pair]["frobenius_delta"],
                        "top_delta_abs": viz["board"][sc]["pairs"][pair]["top_delta_abs"][:3],
                    }
                    for pair in ("TI", "TA", "IA")
                }
                for sc in viz["scenarios"]
            },
        },
        "summary": suite.get("summary", {}),
        "lifts": suite.get("lifts", {}),
        "interpret_board": suite.get("interpret_board", {}),
        "scenarios": suite.get("scenarios", viz_scenarios),
        "methods": suite.get("methods", []),
        "seeds": suite.get("seeds", []),
        "msrvtt_probe": probe,
    }
    write_json(OUT / "fsds_interaction_online_adapt.json", payload)
    write_readme(suite, probe, OUT / "README.md", viz=viz)
    write_readme(suite, probe, DOCS / "FSDS_Interaction_Online_Adapt.md", viz=viz)

    if probe:
        print(
            "msrvtt",
            probe["source"],
            "phi_TI=%.4f" % probe["mean_phi_TI"],
            "r_po=%.4f" % probe.get("r_po", float("nan")),
            "read=",
            probe.get("interpret", {}).get("read", ""),
            flush=True,
        )


if __name__ == "__main__":
    main()
