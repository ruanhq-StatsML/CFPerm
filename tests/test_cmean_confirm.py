"""Tests for confirm + tip-cmean joint formulation prototype."""
from __future__ import annotations

from agod.cmean_confirm import (
    calibrate_cmean_confirm_thresholds,
    run_cmean_confirm_prototype,
    tip_cmean_delta,
)
from agod.loco_po_monitor import make_tip_shift_stream


def test_tip_cmean_delta_rises_at_shift():
    stream, tips = make_tip_shift_stream(
        n_batches=24, batch_size=64, d=6, tip_idx=(0, 1), shift_at=12, seed=0, shift_mean=3.0
    )
    pre = tip_cmean_delta(stream[10][0], stream[11][0], tips)
    at = tip_cmean_delta(stream[11][0], stream[12][0], tips)
    assert at > pre


def test_calibrate_thresholds_positive():
    stream, tips = make_tip_shift_stream(n_batches=20, batch_size=64, tip_idx=(0, 1), seed=1)
    thr = calibrate_cmean_confirm_thresholds(
        stream, tips, burn_in=8, hard_k=5.0, seed=1, n_estimators=15
    )
    assert thr.eps_delta > 0
    assert thr.thr_po >= 0
    assert thr.recipe == "burn_mean_plus_k_sd"
    assert thr.eps_delta >= thr.delta_burn_mean


def test_prototype_hits_known_shift():
    stream, tips = make_tip_shift_stream(
        n_batches=36,
        batch_size=96,
        d=8,
        tip_idx=(0, 1),
        shift_at=18,
        seed=2,
        shift_mean=3.0,
    )
    out = run_cmean_confirm_prototype(
        stream,
        tips,
        burn_in=8,
        hard_k=5.0,
        seed=2,
        n_estimators=15,
        known_shift=18,
    )
    assert out["t_star"] is not None
    assert abs(out["t_star"] - 18) <= 2
    assert out["detection_delay"] is not None
    assert abs(out["detection_delay"]) <= 2
    # thresholds documented
    assert "eps_delta" in out["thresholds"]
    assert out["thresholds"]["recipe"] == "burn_mean_plus_k_sd"


def test_wrong_tips_do_not_fire_at_tip_shift():
    stream, _ = make_tip_shift_stream(
        n_batches=36, batch_size=96, d=8, tip_idx=(0, 1), shift_at=18, seed=3, shift_mean=3.0
    )
    wrong = run_cmean_confirm_prototype(
        stream, (6, 7), burn_in=8, hard_k=5.0, seed=3, n_estimators=15, known_shift=18
    )
    # wrong tips: either miss or not within ±1 of the tip shift
    if wrong["t_star"] is not None:
        assert abs(wrong["t_star"] - 18) > 1
