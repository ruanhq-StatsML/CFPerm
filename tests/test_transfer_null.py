"""Tests for transfer null / excess / probe efficiency."""
from __future__ import annotations

import numpy as np

from agod.transfer_null import (
    excess_auc,
    permute_auc,
    probe_efficiency,
    relative_probe_flops,
    summarize_null_pack,
)
from scripts.run_sample_chunk_adjacent_board import (
    chunk_by_n,
    fs_adjacent,
    interpret_pack,
)


def test_permute_auc_near_half_when_shuffled():
    rng = np.random.default_rng(0)
    y = np.array([0] * 80 + [1] * 20)
    # Strong scores aligned with y → high obs AUC, null ~0.5
    proba = y.astype(float) * 0.9 + rng.uniform(0, 0.05, size=len(y))
    null = permute_auc(y, proba, n_perm=20, seed=1)
    assert 0.35 < null["null_auc_mean"] < 0.65
    obs = float(__import__("sklearn.metrics", fromlist=["roc_auc_score"]).roc_auc_score(y, proba))
    assert excess_auc(obs, null["null_auc_mean"]) > 0.2


def test_probe_eff_scales_with_flops():
    flops = relative_probe_flops(1000, 8)
    assert flops == 1000 * 8 * (60 + 200)
    e1 = probe_efficiency(0.2, flops)
    e2 = probe_efficiency(0.2, flops * 2)
    assert e1 > e2


def test_enrich_and_fs_adjacent_null_fields():
    rng = np.random.default_rng(2)
    n, d = 3000, 6
    X = rng.normal(size=(n, d))
    y = (X[:, 0] + rng.normal(scale=0.3, size=n) > 0.2).astype(float)
    T = chunk_by_n(n, 1000)
    rows = fs_adjacent(
        X,
        y,
        T,
        [f"f{j}" for j in range(d)],
        select_k=3,
        y_quantile=0.7,
        seed=0,
        max_pairs=5,
        min_n=50,
        n_null_perm=4,
    )
    assert len(rows) >= 1
    assert "excess_auc" in rows[0]
    assert "probe_eff" in rows[0]
    s = summarize_null_pack(rows)
    assert s["mean_excess_auc"] is not None


def test_interpret_pack_near_null():
    assert "null" in interpret_pack(0.92, 0.8, mean_excess=0.01)
    assert "stable" in interpret_pack(0.92, 0.8, mean_excess=0.3)


def test_ece_perfect_and_biased():
    from agod.transfer_null import expected_calibration_error

    y = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float)
    p_ok = np.array([0.1, 0.2, 0.2, 0.1, 0.8, 0.9, 0.8, 0.9])
    p_bad = np.ones(8) * 0.9
    assert expected_calibration_error(y, p_ok, n_bins=5) < 0.25
    assert expected_calibration_error(y, p_bad, n_bins=5) > 0.3


def test_fs_adjacent_exposes_ece():
    rng = np.random.default_rng(3)
    n, d = 2500, 5
    X = rng.normal(size=(n, d))
    y = (X[:, 0] > 0).astype(float)
    T = chunk_by_n(n, 500)
    rows = fs_adjacent(
        X,
        y,
        T,
        [f"f{j}" for j in range(d)],
        select_k=2,
        y_quantile=0.7,
        seed=1,
        max_pairs=4,
        min_n=40,
        n_null_perm=3,
    )
    assert any(np.isfinite(r.get("hgb_ece", np.nan)) for r in rows)
    assert any(np.isfinite(r.get("logreg_ece", np.nan)) for r in rows)


def test_blocked_bootstrap_ci_covers_mean():
    from agod.transfer_null import blocked_bootstrap_ci, pack_excess_ci

    rng = np.random.default_rng(0)
    # slowly drifting series → dependence
    vals = list(0.3 + 0.02 * np.arange(20) + rng.normal(0, 0.01, 20))
    ci = blocked_bootstrap_ci(vals, block_size=2, n_boot=300, alpha=0.1, seed=1)
    assert ci["lo"] <= ci["mean"] <= ci["hi"]
    assert ci["hi"] - ci["lo"] > 0
    rows = [{"excess_auc": v} for v in vals]
    pci = pack_excess_ci(rows, seed=2)
    assert np.isfinite(pci["lo"]) and np.isfinite(pci["hi"])


def test_reservoir_sample_uniform_size_and_modes():
    from agod.transfer_null import reservoir_sample_indices, select_pair_indices

    idx = reservoir_sample_indices(100, 10, seed=0)
    assert len(idx) == 10
    assert idx == sorted(idx)
    assert len(set(idx)) == 10
    assert select_pair_indices(50, 50) == list(range(50))
    lin = select_pair_indices(20, 5, mode="linspace", seed=0)
    assert len(lin) == 5
    res = select_pair_indices(20, 5, mode="reservoir", seed=0)
    assert len(res) == 5
    # different seeds → often different reservoirs (smoke, not flaky hard assert)
    r2 = select_pair_indices(20, 5, mode="reservoir", seed=1)
    assert len(r2) == 5
    assert select_pair_indices(20, 5, mode="head") == [0, 1, 2, 3, 4]
