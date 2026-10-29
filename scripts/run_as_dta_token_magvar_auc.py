#!/usr/bin/env python3
"""Token perturb → Mag×Var AUC: benign vs malignant (no multimodal).

  python3 scripts/run_as_dta_token_magvar_auc.py
  python3 scripts/run_as_dta_token_magvar_auc.py --quick
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_token_magvar_auc import run_auc_suite  # noqa: E402

OUT = ROOT / "results" / "as_dta_token_magvar_auc"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    payload = run_auc_suite(quick=bool(args.quick), seed=int(args.seed), out_dir=OUT)
    print("[tok-auc] pooled AUC:", payload["pooled_auc"], flush=True)
    for b in payload["boards"]:
        a = b["auc"]
        print(
            "[tok-auc] %-14s  AUC(sev)=%.3f  AUC(J)=%.3f  AUC(Yflip)=%.3f  AUC(tok)=%.3f"
            % (
                b["board"],
                a.get("severity", float("nan")),
                a["J"],
                a.get("label_flip_frac", float("nan")),
                a["token_top_share"],
            ),
            flush=True,
        )
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_token_magvar_auc.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "token_magvar_auc.png", ARTIFACTS / "as_dta_token_magvar_auc.png")
    print("[tok-auc] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
