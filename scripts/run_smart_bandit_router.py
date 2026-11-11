#!/usr/bin/env python3
"""smart_bandit_router_v1 go-live evaluation.

  python3 scripts/run_smart_bandit_router.py
  python3 scripts/run_smart_bandit_router.py --quick
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from smart_bandit_router import plot_go_live, run_suite, write_readme  # noqa: E402

OUT = ROOT / "results" / "smart_bandit_router"
DOCS = ROOT / "docs" / "method"


def _jsonable(obj):
    import numpy as np

    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    print("smart_bandit_router_v1", "quick" if args.quick else "full", flush=True)
    suite = run_suite(quick=bool(args.quick), seed=int(args.seed))
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    # drop bulky paths from json except synthetic
    slim = {}
    for name, rec in suite.items():
        if "error" in rec:
            slim[name] = rec
            continue
        slim[name] = {
            "dataset": rec["dataset"],
            "meta": rec["meta"],
            "policies": rec["policies"],
            "go_live": rec["go_live"],
        }
        if name == "synthetic_credit":
            slim[name]["paths"] = rec["paths"]

    payload = _jsonable(slim)
    (OUT / "smart_bandit_router.json").write_text(json.dumps(payload, indent=2))
    plot_go_live(suite, OUT / "smart_bandit_router.png")
    write_readme(suite, OUT / "README.md")
    # Method note for evaluators is hand-authored (FSDS/CFPerm granularity):
    # docs/method/Smart_Bandit_FSDS_GoLive.md — do not overwrite from the runner.

    for name, rec in suite.items():
        if "error" in rec:
            print(" ", name, "ERROR", rec["error"], flush=True)
            continue
        g = rec["go_live"]
        s = g["summary_reward"]
        bayes = g.get("bayes_subgroups", {})
        print(
            " ",
            name,
            "verdict=%s" % g.get("verdict", "?"),
            "bandit=%.0f" % s["bandit"],
            "keep=%.0f" % s["keep_baseline"],
            "lift=%+.0f" % s["lift"],
            "bayes_conf=%.0f%%" % (100.0 * bayes.get("share_mass_confident", 0.0)),
            "fair=%s" % g["gates"]["fairness_floor"]["pass"],
            flush=True,
        )
        if g.get("verdict") != "READY":
            fails = [k for k, v in g["gates"].items() if not v["pass"]]
            print("    FAIL gates:", ",".join(fails), flush=True)
            if g.get("verdict_reason"):
                print("    reason:", g["verdict_reason"], flush=True)


if __name__ == "__main__":
    main()
