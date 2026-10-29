#!/usr/bin/env python3
"""FSDS→RL closed loop on typed-drift proxy boards.

  python3 scripts/run_fsds_rl_closed_loop.py
  python3 scripts/run_fsds_rl_closed_loop.py --quick
  python3 scripts/run_fsds_rl_closed_loop.py --compare-meta-iv
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
    CORE_BOARDS,
    DIVERSE_BOARDS,
    METHODS,
    plot_suite,
    run_suite,
    write_json,
    write_readme,
)

OUT = ROOT / "results" / "fsds_rl_closed_loop"
DOCS = ROOT / "docs" / "method"


def _opportunity_lifts(suite):
    """Δ online_value of fsds_full / fsds_meta_iv vs static / naive_boost."""
    lifts = {}
    for board in suite["boards"]:
        tab = suite["table"][board]
        base = tab["static"]["online_value"]["mean"]
        naive = tab["naive_boost"]["online_value"]["mean"]
        full = tab["fsds_full"]["online_value"]["mean"]
        meta = tab.get("fsds_meta_iv", {}).get("online_value", {}).get("mean")
        candidates = [
            (m, tab[m]["online_value"]["mean"])
            for m in (
                "lr_boost_only",
                "reward_only",
                "meta_only",
                "fsds_full",
                "fsds_meta_iv",
            )
            if m in tab
        ]
        cell = {
            "vs_static": float(full - base),
            "vs_naive_boost": float(full - naive),
            "fsds_full": float(full),
            "static": float(base),
            "naive_boost": float(naive),
            "best_module": max(candidates, key=lambda kv: kv[1])[0],
        }
        if meta is not None:
            cell["fsds_meta_iv"] = float(meta)
            cell["meta_iv_vs_static"] = float(meta - base)
            cell["meta_iv_vs_full"] = float(meta - full)
            cell["meta_iv_vs_naive"] = float(meta - naive)
            if "fsds_meta_iv" in tab:
                cell["meta_iv_mix_mse"] = float(tab["fsds_meta_iv"]["mean_mix_mse"]["mean"])
                cell["full_mix_mse"] = float(tab["fsds_full"]["mean_mix_mse"]["mean"])
                cell["meta_iv_n_iv"] = float(tab["fsds_meta_iv"]["mean_n_iv"]["mean"])
                cell["meta_iv_r_po"] = float(tab["fsds_meta_iv"]["mean_uplift_r_po"]["mean"])
                cell["meta_iv_ate"] = float(
                    tab["fsds_meta_iv"].get("mean_uplift_ate", {}).get("mean", float("nan"))
                )
                cell["meta_iv_ess"] = float(
                    tab["fsds_meta_iv"].get("mean_ess_frac", {}).get("mean", float("nan"))
                )
                cell["meta_iv_regret"] = float(
                    tab["fsds_meta_iv"].get("mean_regret", {}).get("mean", float("nan"))
                )
                cell["meta_iv_jaccard"] = float(
                    tab["fsds_meta_iv"].get("mean_iv_jaccard", {}).get("mean", float("nan"))
                )
        lifts[board] = cell
    ranked = sorted(lifts.items(), key=lambda kv: kv[1]["vs_static"], reverse=True)
    ranked_meta = sorted(
        ((b, c) for b, c in lifts.items() if "meta_iv_vs_full" in c),
        key=lambda kv: kv[1]["meta_iv_vs_full"],
        reverse=True,
    )
    return {
        "per_board": lifts,
        "ranked_vs_static": [b for b, _ in ranked],
        "ranked_meta_iv_vs_full": [b for b, _ in ranked_meta],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--n-per", type=int, default=200)
    ap.add_argument("--n-batches", type=int, default=8)
    ap.add_argument(
        "--boards",
        type=str,
        default=None,
        help="comma-separated subset, or 'all' / 'core' / 'diverse'",
    )
    ap.add_argument("--methods", type=str, default=None, help="comma-separated subset")
    ap.add_argument(
        "--diverse",
        action="store_true",
        help="run RecSys / trading / dialogue / multimodal boards",
    )
    ap.add_argument(
        "--compare-meta-iv",
        action="store_true",
        help="focus ablation: static / fsds_full / fsds_meta_iv / naive_boost",
    )
    args = ap.parse_args()

    n_per = 80 if args.quick else int(args.n_per)
    n_batches = 6 if args.quick else int(args.n_batches)
    n_seeds = 2 if args.quick else int(args.seeds)
    if args.boards:
        key = args.boards.strip().lower()
        if key == "all":
            boards = list(BOARDS)
        elif key == "core":
            boards = list(CORE_BOARDS)
        elif key == "diverse":
            boards = list(DIVERSE_BOARDS)
        else:
            boards = [b.strip() for b in args.boards.split(",")]
    elif args.diverse:
        boards = list(DIVERSE_BOARDS)
    else:
        boards = list(BOARDS)
    if args.methods:
        methods = [m.strip() for m in args.methods.split(",")]
    elif args.compare_meta_iv:
        methods = ["static", "fsds_full", "fsds_meta_iv", "naive_boost"]
    else:
        methods = list(METHODS)
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
                "  %s  value=%.4f (%.4f)  w_err=%.4f  mix_mse=%.4f  "
                "n_iv=%.1f  R_PO=%.4f  ATE=%.4f  ESS=%.3f  regret=%.4f  IV-J=%.3f"
                % (
                    method,
                    cell["online_value"]["mean"],
                    cell["online_value"]["sd"],
                    cell["mean_w_err"]["mean"],
                    cell["mean_mix_mse"]["mean"],
                    cell.get("mean_n_iv", {}).get("mean", float("nan")),
                    cell.get("mean_uplift_r_po", {}).get("mean", float("nan")),
                    cell.get("mean_uplift_ate", {}).get("mean", float("nan")),
                    cell.get("mean_ess_frac", {}).get("mean", float("nan")),
                    cell.get("mean_regret", {}).get("mean", float("nan")),
                    cell.get("mean_iv_jaccard", {}).get("mean", float("nan")),
                ),
                flush=True,
            )
        L = lifts["per_board"][board]
        print(
            "  lift fsds_full vs static=%.4f  vs naive=%.4f  best=%s"
            % (L["vs_static"], L["vs_naive_boost"], L["best_module"]),
            flush=True,
        )
        if "meta_iv_vs_full" in L:
            print(
                "  meta_iv vs full=%+.4f  vs static=%+.4f  n_iv=%.1f  R_PO=%.4f  mixΔ=%+.4f"
                % (
                    L["meta_iv_vs_full"],
                    L["meta_iv_vs_static"],
                    L.get("meta_iv_n_iv", float("nan")),
                    L.get("meta_iv_r_po", float("nan")),
                    L.get("meta_iv_mix_mse", float("nan")) - L.get("full_mix_mse", float("nan")),
                ),
                flush=True,
            )
    print("ranked opportunity:", lifts["ranked_vs_static"], flush=True)
    if lifts.get("ranked_meta_iv_vs_full"):
        print("ranked meta_iv vs full:", lifts["ranked_meta_iv_vs_full"], flush=True)
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
