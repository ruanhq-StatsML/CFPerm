#!/usr/bin/env python3
"""Run the continuous trainer prototype on MSR-VTT window features.

  python3 scripts/run_msrvtt_continuous_trainer.py
  python3 scripts/run_msrvtt_continuous_trainer.py --synthetic
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_continuous_trainer import plot_continuous_trainer, run_continuous_trainer
from msrvtt_multimodal_attribution import find_feature_zip, load_window_bundle, make_synthetic_bundle, write_json

OUT = ROOT / "results" / "msrvtt_continuous_trainer"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", type=str, default=None)
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--n-batches", type=int, default=10)
    ap.add_argument("--eta0", type=float, default=0.05)
    ap.add_argument("--steps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    if args.synthetic:
        bundle = make_synthetic_bundle(seed=args.seed)
        print("synthetic", bundle.X.shape, flush=True)
    else:
        zip_path = Path(args.zip) if args.zip else find_feature_zip(ROOT)
        if zip_path is None:
            raise SystemExit("feature_video_audio.zip not found")
        print("loading", zip_path, flush=True)
        bundle = load_window_bundle(zip_path, root=ROOT)

    summary, _ = run_continuous_trainer(
        bundle,
        n_batches=args.n_batches,
        eta0=args.eta0,
        steps_per_batch=args.steps,
        seed=args.seed,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "msrvtt_continuous_trainer.json", {k: v for k, v in summary.items()})
    plot_continuous_trainer(summary, OUT / "msrvtt_continuous_trainer.png")
    print("mean_pi", summary["mean_pi"], "first_batch_acc", summary["first_batch_acc"], flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
