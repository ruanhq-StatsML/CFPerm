"""Tests for DPO + label pollution + decision-path replay."""
from __future__ import annotations

import numpy as np

from sandbox.flywheel_dpo_replay import (
    dpo_preference_loss,
    dpo_update_from_pairs,
    pollute_labels_region,
    replay_with_shuffled_paths,
    sample_period_indices,
    shuffle_decision_path,
)


def test_pollute_labels_permutation_preserves_multiset():
    y = np.arange(20, dtype=float)
    idx = np.array([2, 3, 4, 5, 10, 11])
    y2, meta = pollute_labels_region(y, idx, seed=0)
    assert meta["n_polluted"] == 6
    # outside region unchanged
    mask = np.ones(20, dtype=bool)
    mask[idx] = False
    np.testing.assert_array_equal(y2[mask], y[mask])
    # inside region is a permutation of original values
    assert sorted(y2[idx].tolist()) == sorted(y[idx].tolist())
    # seed 0 on 6 distinct values should move at least one entry
    assert not np.array_equal(y2[idx], y[idx])


def test_sample_periods_and_shuffle_path():
    idx = sample_period_indices(200, n_periods=2, period_len=30, seed=1)
    assert len(idx) == 60
    idx2 = sample_period_indices(
        500, n_periods=2, period_len=40, seed=2, prefer_range=(200, 400)
    )
    assert len(idx2) == 80
    assert idx2.min() >= 200 and idx2.max() < 400
    path = ["observe", "forecast:hgb", "decide:retrain", "refit", "log"]
    sh = shuffle_decision_path(path, seed=2)
    assert sh[0] == "observe" and sh[-1] == "log"
    assert sorted(sh[1:-1]) == sorted(path[1:-1])


def test_replay_and_dpo_loss():
    steps = [
        {
            "step": 0,
            "residual": 0.1,
            "surprise": 0.2,
            "decision_path": ["observe", "forecast:hgb", "decide:idle", "log"],
        },
        {
            "step": 1,
            "residual": -0.2,
            "surprise": 0.3,
            "decision_path": ["observe", "forecast:hgb", "decide:retrain", "refit", "log"],
        },
    ]
    replay = replay_with_shuffled_paths(steps, seed=3)
    assert replay[0]["replay"] is True
    assert "decision_path_orig" in replay[0]
    loss = dpo_preference_loss(1.0, 0.0, beta=1.0)
    assert 0 < loss < 1  # -log σ(1)
    rep = dpo_update_from_pairs([(1.0, 0.0), (0.5, -0.5)])
    assert rep["n_pairs"] == 2
    assert rep["frac_chosen_better"] == 1.0


def test_pollute_replay_dpo_smoke():
    from sandbox.forecast_bakeoff import Pack
    from sandbox.flywheel_dpo_replay import run_pollute_replay_dpo

    rng = np.random.default_rng(0)
    n = 450
    X = rng.normal(size=(n, 2))
    y = np.zeros(n)
    for i in range(1, n):
        y[i] = 0.55 * y[i - 1] + 0.15 * X[i, 0] + rng.normal(scale=0.3)
    pack = Pack(name="toy_dpo", X=X, y=y)
    rep = run_pollute_replay_dpo(
        pack, model_name="ridge", warm=100, max_steps=80, seed=0, period_len=40
    )
    assert rep["n_step_json"] == 80
    assert rep["dpo"]["n_pairs"] == 2
    assert "reading" in rep
