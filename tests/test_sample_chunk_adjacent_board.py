"""Tests for sample-chunk adjacent board helpers."""
from __future__ import annotations

import numpy as np

from scripts.run_sample_chunk_adjacent_board import (
    adjacent_chunk_pairs,
    annotate_stability,
    chunk_by_n,
    fs_adjacent,
    interpret_pack,
    jaccard,
)


def test_chunk_by_n_1000_2000():
    T = chunk_by_n(4500, 1000)
    assert list(np.unique(T)) == [0, 1, 2, 3, 4]
    assert (T == 0).sum() == 1000
    assert (T == 4).sum() == 500
    T2 = chunk_by_n(4500, 2000)
    assert list(np.unique(T2)) == [0, 1, 2]
    assert (T2 == 0).sum() == 2000


def test_jaccard_and_stability_annotation():
    assert jaccard(["a", "b"], ["b", "c"]) == 1 / 3
    rows = annotate_stability(
        [
            {"top_fsds": ["a", "b", "c"], "hgb_auc": 0.7},
            {"top_fsds": ["b", "c", "d"], "hgb_auc": 0.6},
        ]
    )
    assert rows[0]["jaccard_top"] is None
    assert abs(rows[1]["jaccard_top"] - 0.5) < 1e-9


def test_interpret_pack_buckets():
    assert "weak" in interpret_pack(0.55, 0.8)
    assert "stable" in interpret_pack(0.95, 0.7)
    assert "shifting" in interpret_pack(0.95, 0.2)


def test_adjacent_pairs_and_fs_smoke():
    rng = np.random.default_rng(0)
    n, d = 4000, 8
    X = rng.normal(size=(n, d))
    y = rng.normal(size=n)
    y[2000:] += 0.4
    X[2000:, 0] += 0.8
    T = chunk_by_n(n, 1000)
    assert adjacent_chunk_pairs(T) == [(0, 1), (1, 2), (2, 3)]
    names = [f"f{j}" for j in range(d)]
    rows = fs_adjacent(
        X, y, T, names, select_k=4, y_quantile=0.7, seed=0, max_pairs=10, min_n=50
    )
    assert len(rows) >= 2
    assert all("hgb_auc" in r and "delta_Y" in r for r in rows)
    assert rows[0]["jaccard_top"] is None
    assert rows[1]["jaccard_top"] is not None


def test_binary_y_rare_positive():
    """Sparse binary y (like y_convert) must not go through quantile thr=0."""
    rng = np.random.default_rng(1)
    n, d = 4000, 6
    X = rng.normal(size=(n, d))
    y = np.zeros(n)
    y[::40] = 1.0
    X[y == 1, 0] += 1.5
    T = chunk_by_n(n, 1000)
    rows = fs_adjacent(
        X,
        y,
        T,
        [f"f{j}" for j in range(d)],
        select_k=3,
        y_quantile=0.7,
        seed=0,
        max_pairs=10,
        min_n=50,
    )
    assert len(rows) >= 2
    assert all(r["y_threshold"] == 0.5 for r in rows)
    assert all("logreg_auc" in r and "fsds_cmean_jaccard" in r for r in rows)


def test_ship_gate_never_promotes():
    from scripts.run_sample_chunk_adjacent_board import ship_gate_from_board

    g = ship_gate_from_board([])
    assert g["promote_HGB_to_production"] is False
    assert "true_driver" in g["what_needs_other_tools"]
