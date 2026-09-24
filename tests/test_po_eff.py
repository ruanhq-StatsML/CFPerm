"""Tests for PO-risk adaptation efficiency scorecard."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from agod.po_eff import (
    fit_flops,
    mode_relative_flops,
    mse_efficiency,
    rank_efficiency,
    scorecard_from_dataset_block,
    scorecard_from_summary,
)

ROOT = Path(__file__).resolve().parents[1]


def test_flops_ordering_probe_gt_refit_when_few_rejects():
    # Always-on probe over 40 batches vs refit on 5 rejects
    probe = mode_relative_flops("probe", n_batches=40, n_reject=5, batch_size=100)
    refit = mode_relative_flops("refit", n_batches=40, n_reject=5, batch_size=100, n_control=1)
    ref = mode_relative_flops("ref", n_batches=40, n_reject=5, batch_size=100)
    assert probe > refit > ref
    assert fit_flops(100) == 100 * 30 * 6


def test_rank_eff_prefers_cheap_delta():
    # Same Δρ, cheaper flops → higher eff
    e_cheap = rank_efficiency(0.40, 0.30, relative_flops=1e5)
    e_dear = rank_efficiency(0.40, 0.30, relative_flops=1e7)
    assert e_cheap > e_dear


def test_mse_eff_sign():
    assert mse_efficiency(80.0, 100.0, 1e6) > 0
    assert mse_efficiency(120.0, 100.0, 1e6) < 0


def test_scorecard_from_mini_block():
    block = {
        "results": {
            "uniform": {
                "mse_mean_sig": 100.0,
                "gate_duty": 0.2,
                "po_quality": {
                    "n_reject": 4,
                    "ref": {"spearman": 0.2},
                    "probe": {"spearman": 0.35},
                    "refit": {"spearman": 0.33},
                },
            },
            "ref_po": {"mse_mean_sig": 95.0},
            "probe_po": {"mse_mean_sig": 110.0},
            "refit_po": {"mse_mean_sig": 90.0},
        }
    }
    c = scorecard_from_dataset_block(
        "toy", block, n_batches=20, batch_size=50, n_control=1
    )
    assert c["ok"]
    assert abs(c["gate_duty"] - 0.2) < 1e-9
    assert abs(c["budget_ratio_refit_vs_probe"] - 0.2) < 1e-9
    assert c["best_rank_eff_mode"] in ("probe", "refit")
    assert c["best_mse_eff_mode"] == "refit"  # only positive mse_eff spender
    refit = next(m for m in c["modes"] if m["mode"] == "refit")
    assert np.isfinite(refit["rank_eff_expected"])


def test_expected_flops_scales_with_duty():
    from agod.po_eff import expected_adapt_flops, budget_ratio_refit_vs_probe

    low = expected_adapt_flops(
        "refit", gate_duty=0.1, n_batches=40, batch_size=100, n_control=1
    )
    high = expected_adapt_flops(
        "refit", gate_duty=0.5, n_batches=40, batch_size=100, n_control=1
    )
    probe = expected_adapt_flops(
        "probe", gate_duty=0.1, n_batches=40, batch_size=100, n_control=1
    )
    assert high > low
    assert abs(probe - expected_adapt_flops(
        "probe", gate_duty=0.9, n_batches=40, batch_size=100
    )) < 1e-9  # probe duty-invariant
    assert abs(budget_ratio_refit_vs_probe(0.25, n_control=2) - 0.5) < 1e-9
    assert abs(low / probe - 0.1) < 1e-9


def test_scorecard_from_repo_summary_smoke():
    path = ROOT / "results/agod_po_ref_vs_refit/summary.json"
    if not path.exists():
        return
    summary = json.loads(path.read_text())
    card = scorecard_from_summary(summary)
    assert card["n_datasets"] >= 1
    assert "headline" in card
