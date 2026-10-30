#!/usr/bin/env python3
"""CLI for po-risk-var.

  po-risk-var justify --quick
  po-risk-var logic
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path


def main(argv=None):
    # editable / src layout
    pkg_root = Path(__file__).resolve().parents[1]
    if str(pkg_root) not in sys.path:
        sys.path.insert(0, str(pkg_root))

    from po_risk_var._path import REPO_ROOT
    from po_risk_var.product import estimator_priority
    from po_risk_var.suite import run_justify_suite

    ap = argparse.ArgumentParser(prog="po-risk-var")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_log = sub.add_parser("logic", help="print product next-step logic")
    p_jus = sub.add_parser("justify", help="synthetic + real-world apps")
    p_jus.add_argument("--quick", action="store_true")
    p_jus.add_argument("--seed", type=int, default=2026)
    p_jus.add_argument("--method", default="cf", choices=("cf", "bcf", "boot"))
    sub.add_parser("products", help="five product applications suite")

    args = ap.parse_args(argv)

    if args.cmd == "logic":
        steps = [
            "1. fit φ, τ̂ (PermuCATE RF or CausalForestDML)",
            "2. unit risk u = (φ − τ̂)²",
            "3. Mag×Var: μ = E[u], s = Std(sub-env mean u)",
            "4. Var(u): CF leaf-size-debiased + leaf merge (min n=30/50)  ← DEFAULT",
            "5. Y-flip on text integrity boards (attacked rows only)",
            "6. triage → hold / conformal / VAI_freeze",
            "7. honesty: bias_ratio = B_leaf / V_naive",
            "",
            "BCF: rpy2 → R bcf (product Bayesian path)",
            "Bootstrap: intuitive fallback / sanity check",
        ]
        print("PO-risk Var(u) — product next-step logic")
        print("----------------------------------------")
        for s in steps:
            print(s)
        print("\nEstimator priority:")
        for e in estimator_priority():
            print(f"  - {e['name']}  [{e['role']}]  {e['when']}")
        return 0

    if args.cmd == "products":
        from po_risk_var.apps_product import run_product_apps_suite

        out = REPO_ROOT / "results" / "as_dta_po_risk_var_product_apps"
        payload = run_product_apps_suite(out_dir=out)
        docs = REPO_ROOT / "docs" / "method"
        docs.mkdir(parents=True, exist_ok=True)
        shutil.copy(out / "README.md", docs / "AS_DTA_po_risk_var_product_apps.md")
        h = payload["headline"]
        print(
            "[products] routing_save=$%.0f  canary_micro_loss=%s  cascade_blocked=%s"
            % (
                payload["industrial_cost_saved_usd"],
                h["canary_micro_loss"],
                h["agent_cascade_blocked"],
            ),
            flush=True,
        )
        print("[products] wrote", out, flush=True)
        return 0

    if args.cmd == "justify":
        out = REPO_ROOT / "results" / "as_dta_po_risk_var_apps"
        payload = run_justify_suite(
            seed=int(args.seed),
            out_dir=out,
            quick=bool(args.quick),
            method=str(args.method),
        )
        h = payload["headline"]
        print(
            "[po-risk-var] synth=%s  real=%s  meets=%s"
            % (h["synth_n"], h["real_n"], h["meets_expectation"]),
            flush=True,
        )
        docs = REPO_ROOT / "docs" / "method"
        docs.mkdir(parents=True, exist_ok=True)
        shutil.copy(out / "README.md", docs / "AS_DTA_po_risk_var_apps.md")
        art = Path("/opt/cursor/artifacts")
        art.mkdir(parents=True, exist_ok=True)
        if (out / "po_risk_var_apps.png").exists():
            shutil.copy(out / "po_risk_var_apps.png", art / "as_dta_po_risk_var_apps.png")
        print("[po-risk-var] wrote", out, flush=True)
        print(json.dumps(h, indent=2))
        return 0 if h["meets_expectation"] else 1

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
