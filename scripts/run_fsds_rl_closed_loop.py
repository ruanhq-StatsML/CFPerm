#!/usr/bin/env python3
"""FSDS→RL closed loop on typed-drift proxy boards.

  python3 scripts/run_fsds_rl_closed_loop.py
  python3 scripts/run_fsds_rl_closed_loop.py --quick
  python3 scripts/run_fsds_rl_closed_loop.py --boards bandit_mixed,causal_proxy,procgen_proxy
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fsds_rl_closed_loop import (  # noqa: E402
    BOARDS,
    METHODS,
    plot_suite,
    run_suite,
    write_json,
    write_readme,
)

OUT = ROOT / "results" / "fsds_rl_closed_loop"
DOCS = ROOT / "docs" / "method"


def _opportunity_lifts(suite):
    """Δ online_value of fsds_full vs static / naive_boost per board."""
    lifts = {}
    for board in suite["boards"]:
        tab = suite["table"][board]
        base = tab["static"]["online_value"]["mean"]
        naive = tab["naive_boost"]["online_value"]["mean"]
        full = tab["fsds_full"]["online_value"]["mean"]
        lifts[board] = {
            "vs_static": float(full - base),
            "vs_naive_boost": float(full - naive),
            "fsds_full": float(full),
            "static": float(base),
            "naive_boost": float(naive),
            "best_module": max(
                (
                    (m, tab[m]["online_value"]["mean"])
                    for m in ("lr_boost_only", "reward_only", "meta_only", "fsds_full")
                ),
                key=lambda kv: kv[1],
            )[0],
        }
    ranked = sorted(lifts.items(), key=lambda kv: kv[1]["vs_static"], reverse=True)
    return {"per_board": lifts, "ranked_vs_static": [b for b, _ in ranked]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-per", type=int, default=200)
    ap.add_argument("--n-batches", type=int, default=8)
    ap.add_argument("--boards", type=str, default=None, help="comma-separated subset")
    ap.add_argument("--methods", type=str, default=None, help="comma-separated subset")
    args = ap.parse_args()

    n_per = 80 if args.quick else int(args.n_per)
    n_batches = 6 if args.quick else int(args.n_batches)
    n_seeds = 2 if args.quick else int(args.seeds)
    boards = [b.strip() for b in args.boards.split(",")] if args.boards else list(BOARDS)
    methods = [m.strip() for m in args.methods.split(",")] if args.methods else list(METHODS)
    seeds = list(range(2026, 2026 + n_seeds))

    print(
        "fsds-rl", boards, methods, n_seeds, "seeds", n_per, "per", n_batches, "batches",
        flush=True,
    )
    suite = run_suite(
        boards=boards, methods=methods, seeds=seeds, n_per=n_per, n_batches=n_batches
    )
    lifts = _opportunity_lifts(suite)
    suite["lifts"] = lifts

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "fsds_rl_closed_loop.json", suite)
    (OUT / "opportunity_lifts.json").write_text(json.dumps(lifts, indent=2))
    plot_suite(suite, OUT / "fsds_rl_closed_loop.png")
    write_readme(suite, OUT / "README.md")
    write_readme(suite, DOCS / "FSDS_RL_closed_loop.md")

    for board in boards:
        print("==", board, flush=True)
        for method in methods:
            cell = suite["table"][board][method]
            print(
                "  %s  value=%.4f (%.4f)  w_err=%.4f  mix_mse=%.4f"
                % (
                    method,
                    cell["online_value"]["mean"],
                    cell["online_value"]["sd"],
                    cell["mean_w_err"]["mean"],
                    cell["mean_mix_mse"]["mean"],
                ),
                flush=True,
            )
        L = lifts["per_board"][board]
        print(
            "  lift fsds_full vs static=%.4f  vs naive=%.4f  best=%s"
            % (L["vs_static"], L["vs_naive_boost"], L["best_module"]),
            flush=True,
        )
    print("ranked opportunity:", lifts["ranked_vs_static"], flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
