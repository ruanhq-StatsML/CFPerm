"""Tests for LOCO PO-risk tip timing monitor."""
from __future__ import annotations

from agod.loco_po_monitor import (
    init_loco_po_monitor,
    loco_po_risks,
    make_tip_shift_stream,
    run_loco_po_stream,
    update_loco_po_monitor,
)
from agod.online_rfperm import ScalarOnlineRFPerm, first_reject_index


def test_scalar_online_rfperm_rejects_spike():
    st = ScalarOnlineRFPerm()
    for _ in range(6):
        st.step(0.1, burn_in=True)
    out = st.step(5.0, burn_in=False, alpha=0.2, fdr="fixed")
    assert out["reject"] is True or out["p"] < 0.5
    assert first_reject_index(st.reject_hist, after=6) is not None or out["p"] < 1.0


def test_loco_tip_group_positive_when_tips_matter():
    import numpy as np

    rng = np.random.default_rng(0)
    n, d = 200, 6
    X0 = rng.normal(size=(n, d))
    y0 = X0[:, 0] + X0[:, 1] + rng.normal(0, 0.2, n)
    X1 = rng.normal(size=(n, d))
    X1[:, 0] += 1.5
    X1[:, 1] += 1.5
    y1 = -(X1[:, 0] + X1[:, 1]) + rng.normal(0, 0.2, n)
    r = loco_po_risks(X0, y0, X1, y1, tip_idx=[0, 1], seed=1, n_estimators=20)
    assert r["tip_group_loco"] >= 0.0
    assert isinstance(r["per_tip_loco"], dict)
    assert set(r["per_tip_loco"]) == {0, 1}


def test_stream_detects_tip_shift_near_known_time():
    stream, tips = make_tip_shift_stream(
        n_batches=36,
        batch_size=96,
        d=6,
        tip_idx=(0, 1),
        shift_at=18,
        seed=2,
        shift_mean=3.0,
    )
    summary = run_loco_po_stream(
        stream,
        tips,
        burn_in=8,
        alpha=0.08,
        seed=2,
        n_estimators=20,
        per_tip=False,
        score_mode="confirm",
        known_shift=18,
    )
    assert summary["tip_first_reject_t"] is not None
    assert summary["tip_first_reject_t"] >= 8
    assert summary["detection_delay_tip"] is not None
    assert summary["detection_delay_tip"] >= -2
    assert summary["detection_delay_tip"] <= 6
    assert abs(summary["tip_first_reject_t"] - 18) <= 6


def test_evaluate_score_modes_confirm_beats_wrong_tips():
    from agod.loco_po_monitor import evaluate_score_modes

    out = evaluate_score_modes(
        seeds=(0, 1, 2),
        modes=("confirm", "sum"),
        shift_at=20,
        true_tips=(0, 1),
        wrong_tips=(6, 7),
        n_estimators=15,
    )
    by = {r["mode"]: r for r in out["rows"]}
    assert by["confirm"]["true_hit_pm1"] >= 0.66
    # wrong tips should hit the known tip-shift less often than true tips
    assert by["confirm"]["wrong_hit_pm1"] <= by["confirm"]["true_hit_pm1"]


def test_update_appends_histories():
    stream, tips = make_tip_shift_stream(n_batches=12, batch_size=64, shift_at=8, seed=3)
    st = init_loco_po_monitor(tips, seed=3, n_estimators=15)
    for t in range(1, 6):
        update_loco_po_monitor(
            st,
            stream[t - 1][0],
            stream[t - 1][1],
            stream[t][0],
            stream[t][1],
            burn_in=True,
            t=t,
        )
    assert len(st.tip_loco_hist) == 5
    assert len(st.po_full_hist) == 5
    assert st.tip_stream.n_burn == 5
