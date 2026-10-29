#!/usr/bin/env python3
"""Synthetic instrumental-variable selection + sensitivity + subset localization.

  python3 scripts/run_iv_selection_synth.py
  python3 scripts/run_iv_selection_synth.py --quick
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from iv_selection_synth import (  # noqa: E402
    conf_sweep,
    plot_sensitivity_extra,
    plot_suite,
    run_sensitivity_suite,
    run_suite,
    write_json,
    write_readme,
)

OUT = ROOT / "results" / "iv_selection_synth"
DOCS = ROOT / "docs" / "method"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--conf", type=float, default=1.0)
    args = ap.parse_args()

    n_seeds = 3 if args.quick else int(args.seeds)
    n = 800 if args.quick else int(args.n)
    seeds = list(range(2026, 2026 + n_seeds))
    sweep_seeds = list(range(2026, 2026 + (2 if args.quick else 5)))

    print("iv-synth", n_seeds, "seeds", "n=", n, "conf=", args.conf, flush=True)
    suite = run_suite(seeds=seeds, n=n, conf_strength=float(args.conf))
    sweep = conf_sweep(
        conf_grid=[0.0, 0.5, 1.0, 1.5, 2.0],
        seeds=sweep_seeds,
        n=n,
    )
    print("sensitivity + subset localization…", flush=True)
    sens = run_sensitivity_suite(
        seeds=seeds, n=n, conf_strength=float(args.conf)
    )

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {
        "suite": {k: v for k, v in suite.items() if k != "score_book"},
        "sweep": sweep,
        "sensitivity": {
            k: v
            for k, v in sens["sensitivity"].items()
            if k not in ("rows",)
        },
        "subset": {
            k: v for k, v in sens["subset"].items() if k != "rows"
        },
    }
    payload["suite"]["score_example"] = suite["score_book"][0]["scores"]
    write_json(OUT / "iv_selection_synth.json", payload)
    plot_suite(suite, sweep, OUT / "iv_selection_synth.png", sens=sens)
    plot_sensitivity_extra(sens, OUT / "iv_sensitivity_extra.png")
    write_readme(suite, sweep, OUT / "README.md", sens=sens)
    write_readme(suite, sweep, DOCS / "IV_selection_synth.md", sens=sens)

    print("true_ate", suite["true_ate"], flush=True)
    print(
        "selection  P=%.3f  R=%.3f  F1=%.3f"
        % (
            suite["selection"]["precision"]["mean"],
            suite["selection"]["recall"]["mean"],
            suite["selection"]["f1"]["mean"],
        ),
        flush=True,
    )
    for m in suite["methods"]:
        t = suite["table"][m]
        print(
            "  %s  abs_err=%.3f (%.3f)  ATE=%.3f  nZ=%.1f  valid/invalid=%.1f/%.1f"
            % (
                m,
                t["abs_err"]["mean"],
                t["abs_err"]["sd"],
                t["ate_hat"]["mean"],
                t["mean_n_z"],
                t["mean_n_valid"],
                t["mean_n_invalid"],
            ),
            flush=True,
        )
    S = sens["sensitivity"]
    print(
        "E-value=%.3f  Γ*=%.2f  ρ*=%.3f  OVB_RV=%.3f (×obs=%.1f)  "
        "Conley_γ*/σ=%.2f  crossed_null_frac=%.2f"
        % (
            S["evalue_point"]["mean"],
            S["gamma_star"]["mean"],
            S.get("rho_star", {}).get("mean", float("nan")),
            S.get("ovb_rv", {}).get("mean", float("nan")),
            S.get("ovb_times_obs", {}).get("mean", float("nan")),
            S.get("conley_gamma_null_over_sigma", {}).get("mean", float("nan")),
            S.get("conley_crossed_frac", 0.0),
        ),
        flush=True,
    )
    print(
        "subset F het ratio=%.2fx" % sens["subset"]["f_heterogeneity_ratio"],
        flush=True,
    )
    for name in ("Q1_low", "Q2", "Q3", "Q4_high"):
        t = sens["subset"]["by_segment"][name]
        print(
            "  %s  F=%.1f  R2=%.3f  ATE=%.3f  scale=%.2f"
            % (
                name,
                t["stage1_f"]["mean"],
                t["stage1_r2"]["mean"],
                t["ate_hat"]["mean"],
                t["mean_seg_scale"]["mean"],
            ),
            flush=True,
        )
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
