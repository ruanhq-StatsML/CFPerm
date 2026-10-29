#!/usr/bin/env python3
"""客群/region Mag×Var on PermuCATE PO.

  python3 scripts/run_as_dta_segment_magvar.py --quick
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_segment_magvar import run_segment_suite  # noqa: E402

OUT = ROOT / "results" / "as_dta_segment_magvar"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    payload = run_segment_suite(quick=bool(args.quick), seed=int(args.seed), out_dir=OUT)
    print(
        "[seg-magvar] hit_rate=%.2f  (%d/%d)"
        % (payload["segment_hit_rate"], payload["n_separated"], payload["n_segments"]),
        flush=True,
    )
    for s in payload["segments"]:
        print(
            "[seg-magvar] %-10s  ben=%.3f  mal=%.3f  gap=%.3f  %s"
            % (
                s["segment"],
                s["benign"]["severity"],
                s["malignant"]["severity"],
                s["severity_gap"],
                "OK" if s["separates_benign_malignant"] else "MISS",
            ),
            flush=True,
        )
    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_segment_magvar.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "segment_magvar.png", ARTIFACTS / "as_dta_segment_magvar.png")
    print("[seg-magvar] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
