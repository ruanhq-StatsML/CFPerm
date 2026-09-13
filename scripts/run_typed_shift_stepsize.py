#!/usr/bin/env python3
"""Run typed shift stepsize vs LR schedulers.

  python3 scripts/run_typed_shift_stepsize.py
  python3 scripts/run_typed_shift_stepsize.py --quick
  python3 scripts/run_typed_shift_stepsize.py --no-bundle
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_multimodal_attribution import find_feature_zip, load_window_bundle, write_json
from typed_shift_stepsize import (
    METHODS,
    graft_midclip_concept,
    plot_comparison,
    plot_eta_paths,
    run_bundle_suite,
    run_suite,
    write_tex_table,
)

OUT = ROOT / "results" / "typed_shift_stepsize"
DOCS = ROOT / "docs" / "method"


def _strip_traces(suite):
    return {k: v for k, v in suite.items() if k != "traces"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-bundle", action="store_true")
    ap.add_argument("--zip", type=str, default=None)
    ap.add_argument("--seeds", type=int, default=8)
    args = ap.parse_args()

    n_batches = 8 if args.quick else 12
    n_per = 48 if args.quick else 64
    n_seeds = 3 if args.quick else args.seeds
    methods = list(METHODS) if not args.quick else ["constant", "cosine", "plateau_m", "restart_m", "fsds_pi", "tss", "oracle_tss"]
    seeds = list(range(2026, 2026 + n_seeds))

    print("suite", n_seeds, "seeds", n_batches, "batches", methods, flush=True)
    suite = run_suite(seeds=seeds, n_batches=n_batches, n_per=n_per, methods=methods)
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "tss_vs_schedulers.json", _strip_traces(suite))
    plot_comparison(suite, OUT / "tss_vs_schedulers.png")
    plot_eta_paths(suite, OUT / "tss_eta_paths.png")
    write_tex_table(suite, OUT / "TSS_vs_schedulers.tex")
    write_tex_table(suite, DOCS / "TSS_vs_schedulers.tex")

    print("identification", suite["identification"], flush=True)
    for regime, rec in suite["table"].items():
        print("==", regime, flush=True)
        for method, cell in rec.items():
            print(
                "  %s  acc=%.3f  post=%.3f  bwt=%.3f  eta=(v %.4f a %.4f t %.4f)"
                % (
                    method,
                    cell["online_acc"]["mean"],
                    cell["post_acc"]["mean"],
                    cell["bwt"]["mean"],
                    cell["mean_lr_video"]["mean"],
                    cell.get("mean_lr_audio", {"mean": float("nan")})["mean"],
                    cell.get("mean_lr_text", {"mean": float("nan")})["mean"],
                ),
                flush=True,
            )

    if not args.no_bundle:
        zip_path = Path(args.zip) if args.zip else find_feature_zip(ROOT)
        if zip_path is None:
            print("no zip; skip bundle illustration", flush=True)
        else:
            print("bundle", zip_path, flush=True)
            bundle = load_window_bundle(zip_path, root=ROOT)
            bsuite = run_bundle_suite(
                bundle,
                n_batches=8 if args.quick else 10,
                methods=[m for m in methods if m != "oracle_tss"],
            )
            slim = {
                m: {
                    **{k: v for k, v in rec.items() if k != "history"},
                    "last_lr": rec["history"][-1]["lr"],
                    "last_c": rec["history"][-1]["c"],
                    "last_delta": rec["history"][-1]["delta"],
                }
                for m, rec in bsuite.items()
            }
            write_json(OUT / "tss_msrvtt_illustration.json", slim)
            print("bundle last_acc", {m: slim[m]["last_acc"] for m in slim}, flush=True)
            print("bundle bwt", {m: slim[m]["bwt"] for m in slim}, flush=True)

            y_graft = graft_midclip_concept(bundle)
            gsuite = run_bundle_suite(
                bundle,
                n_batches=8 if args.quick else 10,
                methods=[m for m in methods if m != "oracle_tss"],
                y=y_graft,
                label="msrvtt_grafted_concept",
            )
            gslim = {
                m: {
                    **{k: v for k, v in rec.items() if k != "history"},
                    "last_lr": rec["history"][-1]["lr"],
                    "last_c": rec["history"][-1]["c"],
                    "last_delta": rec["history"][-1]["delta"],
                }
                for m, rec in gsuite.items()
            }
            write_json(OUT / "tss_msrvtt_grafted_concept.json", gslim)
            print("grafted last_acc", {m: gslim[m]["last_acc"] for m in gslim}, flush=True)
            print("grafted bwt", {m: gslim[m]["bwt"] for m in gslim}, flush=True)
            print("grafted post_acc", {m: gslim[m]["post_acc"] for m in gslim}, flush=True)

    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
