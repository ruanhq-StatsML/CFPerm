#!/usr/bin/env python3
"""MVP: six real-data watermark kinds → strength → AS-DTA action + viz.

  python3 scripts/run_as_dta_watermark_mvp.py
  python3 scripts/run_as_dta_watermark_mvp.py --n-per-kind 100
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_watermark_mvp import KINDS, run_mvp  # noqa: E402

OUT = ROOT / "results" / "as_dta_watermark_mvp"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-per-kind", type=int, default=80)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    payload = run_mvp(n_per_kind=int(args.n_per_kind), seed=int(args.seed), out_dir=OUT)
    summary = payload["summary"]
    print("[wm-mvp] kinds:", list(KINDS), flush=True)
    for k in KINDS:
        b = summary["by_kind"][k]
        print(
            "  %-14s n=%d  pos=%d gray=%d neg=%d miss=%d  vague=%.2f  actions=%s"
            % (
                k,
                b["n"],
                b["strength"].get("strong_pos", 0),
                b["strength"].get("gray", 0),
                b["strength"].get("strong_neg", 0),
                b["strength"].get("missing", 0),
                b["vague_rate"],
                b["actions"],
            ),
            flush=True,
        )
    print("[wm-mvp] same:", *summary["same_diff"]["same"], sep="\n  - ", flush=True)
    print("[wm-mvp] diff:", *summary["same_diff"]["diff"], sep="\n  - ", flush=True)

    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_watermark_mvp.md")
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "watermark_mvp_kinds.png", ARTIFACTS / "as_dta_watermark_mvp_kinds.png")
    print("[wm-mvp] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
