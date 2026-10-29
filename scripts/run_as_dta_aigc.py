#!/usr/bin/env python3
"""Run AS-DTA on HC3: fine-tune → LIVE shifts → PO-risk Peaks → diagnosis board.

  python3 scripts/run_as_dta_aigc.py --quick
  python3 scripts/run_as_dta_aigc.py
  python3 scripts/run_as_dta_aigc.py --max-pairs 3000 --attack-max 220
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_aigc import run_as_dta  # noqa: E402

OUT = ROOT / "results" / "as_dta_aigc"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--max-pairs", type=int, default=2500)
    ap.add_argument("--K", type=int, default=30, help="bootstrap sub-batches")
    ap.add_argument("--synonym-rate", type=float, default=0.30)
    ap.add_argument("--attack-rate", type=float, default=0.40)
    ap.add_argument("--attack-max", type=int, default=220)
    ap.add_argument("--use-marian", action="store_true", help="real en↔zh back-translate")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    payload = run_as_dta(
        max_pairs=int(args.max_pairs),
        K_subbatches=int(args.K),
        synonym_rate=float(args.synonym_rate),
        attack_rate=float(args.attack_rate),
        attack_max=int(args.attack_max),
        use_marian=bool(args.use_marian),
        seed=int(args.seed),
        out_dir=OUT,
        quick=bool(args.quick),
    )

    # mirror docs + artifacts
    DOCS.mkdir(parents=True, exist_ok=True)
    readme = OUT / "README.md"
    if readme.exists():
        shutil.copy(readme, DOCS / "AS_DTA_AIGC_diagnosis.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    png = OUT / "as_dta_diagnosis_dashboard.png"
    html = OUT / "as_dta_diagnosis_board.html"
    if png.exists():
        shutil.copy(png, ARTIFACTS / "as_dta_diagnosis_dashboard.png")
    if html.exists():
        shutil.copy(html, ARTIFACTS / "as_dta_diagnosis_board.html")

    print("routing:", payload["routing"], flush=True)
    print("dashboard:", png, flush=True)


if __name__ == "__main__":
    main()
