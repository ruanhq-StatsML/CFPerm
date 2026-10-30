#!/usr/bin/env python3
"""PO-risk Var(u) estimators: CF debiased / BCF / bootstrap.

  python3 scripts/run_as_dta_po_risk_var.py --quick
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_po_risk_var import run_po_risk_var_suite  # noqa: E402

OUT = ROOT / "results" / "as_dta_po_risk_var"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    payload = run_po_risk_var_suite(
        seed=int(args.seed), out_dir=OUT, quick=bool(args.quick)
    )
    print(
        "[po-var] hit_rate=%.2f  (%d/%d)  bcf_backend=%s"
        % (
            payload["hit_rate"],
            payload["n_estimators_separating"],
            payload["n_estimators_total"],
            payload.get("bcf_backend_used"),
        ),
        flush=True,
    )
    for k, v in payload["separation"].items():
        print(
            "[po-var] %-16s  ben=%.4g  mal=%.4g  gap=%.4g  ratio=%.3g  %s"
            % (
                k,
                v["benign"],
                v["malignant"],
                v["gap"],
                v["ratio"] if v["ratio"] == v["ratio"] else float("nan"),
                "OK" if v["separates"] else "MISS",
            ),
            flush=True,
        )
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_po_risk_var_estimators.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "po_risk_var.png", ARTIFACTS / "as_dta_po_risk_var.png")
    print("[po-var] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
