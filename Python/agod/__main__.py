"""python -m agod --quick   or   python -m agod --out paper/results"""

from __future__ import annotations

import argparse
from pathlib import Path

from .chronoberg import SyntheticChronoBergConfig
from .experiment import plot_result, run_baseline_comparison, save_result
from .online import AGODConfig


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Run AGOD B1/B2/B3 on a ChronoBerg-style stream.")
    p.add_argument("--out", type=Path, default=Path("paper/results"))
    p.add_argument("--n-per-window", type=int, default=80)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quick", action="store_true", help="Tiny CPU run for smoke tests.")
    args = p.parse_args(argv)

    n = 48 if args.quick else args.n_per_window
    steps = 4 if args.quick else 8
    pretrain = 18 if args.quick else 55
    cfg = SyntheticChronoBergConfig(n_per_window=n, seed=args.seed)
    tcfg = AGODConfig(
        seed=args.seed,
        steps_per_window=steps,
        pretrain_steps=pretrain,
        teacher_dim=10 if args.quick else 12,
    )
    result = run_baseline_comparison(config=cfg, trainer_config=tcfg)
    args.out.mkdir(parents=True, exist_ok=True)
    save_result(result, args.out / "agod_metrics.json")
    paths = plot_result(result, args.out)
    print("ranking:", " > ".join(result.ranking))
    for name, metrics in result.metrics.items():
        print(f"{name}: drift-recall={metrics['drift_subgroup_recall']:.3f} "
              f"audio-alpha={metrics['audio_alpha_on_drift']:.3f} "
              f"audio-align={metrics['audio_alignment_on_drift']:.3f}")
    for label, path in paths.items():
        print(f"wrote {label}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
