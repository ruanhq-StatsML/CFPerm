"""Tests for 180-agent continuous-learning fleet."""
from __future__ import annotations

import numpy as np

from agod.agent_fleet import (
    DEFAULT_COHORTS,
    antifraud_useful_rate,
    build_fleet,
    compare_even_vs_pile,
    continuous_learn_tick,
    cuped_variance_ratio,
    even_force_assign_wip,
    fairness_report,
    gini,
    hhi,
    pile_on_baseline,
)


def test_default_fleet_sums_to_180():
    assert sum(DEFAULT_COHORTS.values()) == 180
    st = build_fleet()
    assert st.n == 180
    assert st.agents[0].cohort == "antifraud"
    assert sum(1 for a in st.agents if a.cohort == "cuped") == 10


def test_gini_even_vs_uneven():
    assert gini([1, 1, 1, 1]) < 0.05
    assert gini([10, 0, 0, 0]) > 0.5


def test_hhi_uniform():
    assert abs(hhi([1, 1, 1, 1]) - 0.25) < 1e-9


def test_cuped_reduces_variance_when_x_predicts_y():
    rng = np.random.default_rng(1)
    x = rng.normal(size=800)
    y = 0.7 * x + rng.normal(scale=0.5, size=800)
    rep = cuped_variance_ratio(y, x)
    assert rep["ratio"] < 1.0
    assert rep["ratio"] < 0.85


def test_antifraud_useful_lift():
    u = [1] * 40 + [0] * 10
    rep = antifraud_useful_rate(u, baseline=0.35)
    assert rep["rate"] == 0.8
    assert rep["lift_vs_baseline"] > 0


def test_even_force_beats_pile_on_fairness():
    even = build_fleet()
    even_force_assign_wip(even, seed=0)
    pile = pile_on_baseline(180, seed=0)
    fe, fp = fairness_report(even), fairness_report(pile)
    # pile-on concentrates headcount on magnet → higher HHI
    assert fp["hhi_headcount"] > fe["hhi_headcount"]


def test_compare_scorecard_smoke():
    rep = compare_even_vs_pile(seed=0)
    assert "verdict" in rep
    assert rep["even"]["fairness"]["n_agents"] == 180
    assert rep["even"]["tick"]["cuped"]["ratio"] < 1.0
    assert rep["even"]["mean_score"] >= 0


def test_cl_tick_reflection_keys():
    st = build_fleet()
    even_force_assign_wip(st, seed=2)
    rng = np.random.default_rng(2)
    snap = continuous_learn_tick(
        st,
        fraud_useful=(rng.random(50) < 0.6).astype(int),
        cuped_y=0.5 * rng.normal(size=200) + rng.normal(scale=0.4, size=200),
        cuped_x=rng.normal(size=200),
        seed=2,
    )
    assert "reflection" in snap
    assert snap["reflection"]["action"] in (
        "keep_specialty_mix",
        "rehome_flex_to_underserved",
    )


def test_multi_tick_cl_demo_smoke():
    from agod.agent_fleet import multi_tick_cl_demo

    demo = multi_tick_cl_demo(n_ticks=5, seed=0)
    assert demo["n_ticks"] == 5
    assert len(demo["actions"]) == 5
    assert "reading" in demo
    # later ticks should show stronger CUPED (lower ratio) than early weak ones
    assert demo["cuped_ratios"][-1] < demo["cuped_ratios"][0]
