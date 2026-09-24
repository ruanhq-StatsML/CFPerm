"""Smoke test for tomorrow demo plot script."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]


def test_plot_rsi_tomorrow_demo_smoke(tmp_path):
    from scripts.plot_rsi_tomorrow_demo import main
    import sys

    out = tmp_path / "demo.png"
    argv = sys.argv
    try:
        sys.argv = [
            "plot_rsi_tomorrow_demo.py",
            "--eff",
            str(ROOT / "results/agod_po_eff/po_eff_scorecard.json"),
            "--power",
            str(ROOT / "results/agod_po_power_eff/po_power_eff.json"),
            "--freeze",
            str(ROOT / "results/agod_freeze_eff/freeze_eff.json"),
            "--out",
            str(out),
        ]
        main()
    finally:
        sys.argv = argv
    assert out.exists() and out.stat().st_size > 1000
