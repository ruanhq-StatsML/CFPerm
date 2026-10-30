#!/usr/bin/env python3
"""Synthetic + real-world justify for po-risk-var product package.

  python3 scripts/run_as_dta_po_risk_var_apps.py --quick
  # same as: po-risk-var justify --quick
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "po_risk_var"))
sys.path.insert(0, str(ROOT / "Python" / "src"))

from po_risk_var.cli import main

if __name__ == "__main__":
    raise SystemExit(main(["justify", "--quick"] if len(sys.argv) == 1 else None))
