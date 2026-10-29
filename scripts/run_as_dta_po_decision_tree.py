#!/usr/bin/env python3
"""PO Mag×Var product decision tree (bootstrap var → mean / kurtosis → VAI).

  python3 scripts/run_as_dta_po_decision_tree.py
  python3 scripts/run_as_dta_po_decision_tree.py --quick
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_po_decision_tree import run_product_tree  # noqa: E402

OUT = ROOT / "results" / "as_dta_po_decision_tree"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    payload = run_product_tree(quick=bool(args.quick), seed=int(args.seed), out_dir=OUT)
    print(
        "[po-tree] family_acc=%.2f branch_acc=%.2f"
        % (payload["mean_family_acc"], payload["mean_branch_acc"]),
        flush=True,
    )
    for b in payload["boards"]:
        print("[po-tree]", b["board"], "family=%.2f branch=%.2f" % (b["family_accuracy"], b["branch_accuracy"]), flush=True)
        for r in b["rows"]:
            print(
                "         %-18s → %-18s  d_mu=%.3f d_s=%.3f kurt=%.2f vai=%.2f  [%s]"
                % (
                    r["oracle"],
                    r["pred_branch"],
                    r["decision"]["d_mu"],
                    r["decision"]["d_s"],
                    r["live"]["kurtosis"],
                    r["vai"]["top_share"],
                    "OK" if r["correct_branch"] else ("famOK" if r["correct_family"] else "MISS"),
                ),
                flush=True,
            )

    DOCS.mkdir(parents=True, exist_ok=True)
    # product justify doc
    doc = DOCS / "AS_DTA_po_decision_tree_product.md"
    doc.write_text(
        (OUT / "README.md").read_text()
        + "\n## Justification (landing)\n\n"
        + "Bootstrap Mag×Var is the **triage sensor**. Stable variance + mean lift → cheap "
        + "benign medicine (conformal / wait for label delay). Unstable variance + kurtosis>10 "
        + "→ sparse attack → VAI localizes token/feature blocks to freeze. Otherwise dense/"
        + "schema → heavier review. This is the productized form of Magnitude×Variance — "
        + "not a single OOD alarm.\n"
    )
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "po_decision_tree.png", ARTIFACTS / "as_dta_po_decision_tree.png")
    print("[po-tree] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
