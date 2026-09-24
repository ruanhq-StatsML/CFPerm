"""Tests for sqrt vs cbrt IPTW power efficiency."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from agod.po_power_eff import power_card_from_dataset, power_scorecard_from_summary

ROOT = Path(__file__).resolve().parents[1]


def test_soft_win_and_best_uniform():
    block = {
        "results": {
            "uniform": {"mse_mean_sig": 100.0, "gate_duty": 0.2},
            "sqrt": {"mse_mean_sig": 120.0},
            "cbrt": {"mse_mean_sig": 110.0},
            "gated_sqrt": {"mse_mean_sig": 105.0},
            "gated_cbrt": {"mse_mean_sig": 102.0},
            "dre": {"mse_mean_sig": 130.0},
        }
    }
    c = power_card_from_dataset("toy", block, n_batches=40, batch_size=100)
    assert c["soft_win_cbrt_le_sqrt"] is True
    assert c["best_sig_mode"] == "uniform"
    assert c["soft_mse_eff_gap_cbrt_minus_sqrt"] > 0  # ∛ hurts less per FLOP


def test_repo_cbrt_summary_smoke():
    path = ROOT / "results/agod_po_cbrt/summary.json"
    if not path.exists():
        return
    card = power_scorecard_from_summary(json.loads(path.read_text()))
    assert card["n_datasets"] >= 1
    assert 0.0 <= card["soft_win_rate_cbrt_le_sqrt"] <= 1.0
    assert "median_soft_mse_eff_gap_cbrt_minus_sqrt" in card
    # median should be more robust than mean when one pack dominates
    assert np.isfinite(card["median_soft_mse_eff_gap_cbrt_minus_sqrt"])
