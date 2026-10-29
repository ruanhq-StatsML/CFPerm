#!/usr/bin/env python3
"""Mag×Var multi-board eval — 2×2 contingency + metrics.

  python3 scripts/run_as_dta_mag_var_boards.py
  python3 scripts/run_as_dta_mag_var_boards.py --quick
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_mag_var_eval import METRICS_LIST, run_eval  # noqa: E402
from as_dta_shift_data import BOARD_NAMES  # noqa: E402

OUT = ROOT / "results" / "as_dta_mag_var_boards"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--boards", nargs="*", default=None)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    live_n = 100 if args.quick else 180
    K = 16 if args.quick else 28
    max_pairs = 600 if args.quick else 1000
    names = list(args.boards) if args.boards else list(BOARD_NAMES)

    payload = run_eval(
        board_names=names, live_n=live_n, K=K, seed=int(args.seed), max_pairs=max_pairs, out_dir=OUT
    )
    agg = payload["aggregate"]
    c = agg["contingency_2x2"]
    print("[magvar] metrics:", ", ".join(METRICS_LIST), flush=True)
    print(
        "[magvar] pooled 2x2  [[%d,%d],[%d,%d]]  acc=%.2f bal=%.2f mal_rec=%.2f"
        % (c["tn"], c["fp"], c["fn"], c["tp"], agg["pooled"]["accuracy"],
           agg["pooled"]["balanced_accuracy"], agg["pooled"]["malignant_recall"]),
        flush=True,
    )
    for b in payload["boards"]:
        m = b["metrics"]
        print(
            "[magvar] %-14s acc=%.2f J_gap=%.3f  d_s=%.4f/%.4f  mag_only=%.2f"
            % (b["board"], m["accuracy"], m["J_gap"], m["d_s_benign"], m["d_s_malignant"],
               m["magnitude_only_acc"]),
            flush=True,
        )

    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_mag_var_boards_eval.md")
    method = DOCS / "AS_DTA_mag_var_boards_eval.md"
    method.write_text(
        method.read_text()
        + "\n## Formulate\n\n"
        + "- Data: `Python/src/as_dta_shift_data.py` (`load_board` / `BOARD_NAMES`)\n"
        + "- Eval: `Python/src/as_dta_mag_var_eval.py` (2×2 + metrics list)\n"
        + "- Logic: `docs/method/AS_DTA_magnitude_variance_logic.md`\n"
    )
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "mag_var_boards_eval.png", ARTIFACTS / "as_dta_mag_var_boards_eval.png")
    print("[magvar] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
