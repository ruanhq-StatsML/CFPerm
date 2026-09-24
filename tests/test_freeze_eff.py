"""Tests for freeze MSE–FLOPs efficiency."""
from __future__ import annotations

from agod.freeze_eff import (
    dominates_always,
    freeze_efficiency,
    scorecard_canonical,
    scorecard_from_freeze_bundle,
)


def test_freeze_eff_and_pareto():
    assert freeze_efficiency(0.86, 0.7) > 0
    assert freeze_efficiency(1.15, 0.57) < 0
    assert dominates_always(0.86, 0.7) is True
    assert dominates_always(1.15, 0.57) is False


def test_canonical_scorecard():
    c = scorecard_canonical()
    assert c["n_datasets"] == 2
    assert c["n_pareto_dominances"] >= 1
    assert "freeze_low_share" in c["best_freeze_eff_counts"] or "freeze_early" in c[
        "best_freeze_eff_counts"
    ]


def test_bundle_parser():
    bundle = {
        "electricity": {
            "_agg": {
                "always_adapt": {"mse_post_mean": 1.0, "flops_mean": 100.0},
                "freeze_early": {"mse_post_mean": 0.9, "flops_mean": 60.0},
            }
        }
    }
    c = scorecard_from_freeze_bundle(bundle)
    assert c["cards"][0]["dominates_always"] == ["freeze_early"]
