#!/usr/bin/env python3
"""Attach / repair ``direction`` on localize summary.json (no Drill change).

Closes the gap: stale summaries lack tip_signs → cards always show sign=0 / flat
story with blank tip polarity. Fills tip_signs from feature_shift_diagnostics.csv;
leaves Dy missing→flat unless already present (does not invent S1 from window rates).

  PYTHONPATH=. python3 scripts/tencent_gr/attach_direction_to_summary.py \\
    --summary results/tencent_gr_w1w2_mmd_po_fsds/summary.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from direction_report import ensure_direction, scenario_from_direction  # noqa: E402
from feature_methods_panel import panel_from_dir  # noqa: E402


def attach(
    summary_path: Path,
    *,
    feat_diag_path: Path | None = None,
    inplace: bool = True,
) -> dict:
    blob = json.loads(summary_path.read_text())
    diag = feat_diag_path or (summary_path.parent / "feature_shift_diagnostics.csv")
    direction = ensure_direction(blob, feat_diag_path=diag if diag.exists() else None)
    scenario = scenario_from_direction(direction)
    tips = list((direction.get("tip_signs") or {}).keys())
    if not tips:
        for key in ("fsds_W1_holdout", "fsds_W1", "fsds_W2_temporal", "fsds_W2"):
            tops = (blob.get(key) or {}).get("top_features")
            if isinstance(tops, list) and tops:
                tips = [str(x) for x in tops]
                break
    panel = panel_from_dir(
        summary_path.parent,
        tip_features=tips,
        tip_signs=direction.get("tip_signs") or {},
    )
    blob["direction"] = direction
    blob["scenario"] = scenario
    blob["feature_methods"] = panel
    if inplace:
        summary_path.write_text(json.dumps(blob, indent=2, ensure_ascii=False, default=str) + "\n")
    return {
        "summary": str(summary_path),
        "sign_Dy": direction.get("sign_Dy"),
        "dy_missing": direction.get("dy_missing"),
        "n_tip_signed": sum(
            1 for v in (direction.get("tip_signs") or {}).values() if v in ("+", "-")
        ),
        "scenario": scenario,
        "feature_methods_read": panel.get("read"),
        "consensus_top": panel.get("consensus_top"),
        "report": direction.get("report"),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Attach direction (+ scenario) to summary.json")
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--feat-diag", type=Path, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    info = attach(args.summary, feat_diag_path=args.feat_diag, inplace=not args.dry_run)
    print(json.dumps(info, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
