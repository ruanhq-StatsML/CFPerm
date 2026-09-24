"""Tests for stats-empowered Tree-of-Thoughts over efficiency policies."""
from __future__ import annotations

from agod.tot_eff import beam_search_policy_tot, tot_from_pack_metrics


def test_tot_prefers_refit_under_low_duty_and_soften():
    # Low duty → probe pruned; soft IPTW hurts → SOFTEN_ONLY ∛ leaf exists
    rep = beam_search_policy_tot(
        gate_duty=0.2,
        rank_effs={"ref": 0.0, "probe": 0.3, "refit": 2.0},
        mse_effs={"ref": 0.0, "probe": -0.1, "refit": -0.05},
        rel_sqrt=1.1,
        rel_cbrt=1.05,
        beam_k=3,
    )
    assert "probe" not in " ".join(rep["beam_adapt"]) or "refit" in " ".join(
        rep["beam_adapt"]
    )
    burns = {L["burn"] for L in rep["leaves"]}
    assert "SOFTEN_ONLY" in burns or "KEEP_UNIFORM" in burns
    assert "BURN_SQRT" not in burns


def test_tot_burns_sqrt_when_beats_uniform():
    rep = beam_search_policy_tot(
        gate_duty=0.18,
        rank_effs={"ref": 0.0, "probe": 0.1, "refit": 0.5},
        mse_effs={"ref": 0.0, "probe": 0.0, "refit": 1.0},
        rel_sqrt=0.93,
        rel_cbrt=0.95,
        mse_eff_sqrt=1.0,
        mse_eff_cbrt=0.8,
        beam_k=3,
    )
    burns = {L["burn"] for L in rep["leaves"]}
    assert "BURN_SQRT" in burns


def test_tot_from_pack_smoke():
    rep = tot_from_pack_metrics(
        gate_duty=0.24,
        rank_refit=2.0,
        rank_probe=0.5,
        mse_refit=-0.01,
        mse_probe=-0.02,
        rel_sqrt=1.08,
        rel_cbrt=1.04,
    )
    assert rep["best"] is not None
    assert "reading" in rep
