"""Tests for sample-chunk adjacent board helpers."""
from __future__ import annotations

import numpy as np

from scripts.run_sample_chunk_adjacent_board import (
    adjacent_chunk_pairs,
    chunk_by_n,
    fs_adjacent,
)


def test_chunk_by_n_1000_2000():
    T = chunk_by_n(4500, 1000)
    assert list(np.unique(T)) == [0, 1, 2, 3, 4]
    assert (T == 0).sum() == 1000
    assert (T == 4).sum() == 500
    T2 = chunk_by_n(4500, 2000)
    assert list(np.unique(T2)) == [0, 1, 2]
    assert (T2 == 0).sum() == 2000


def test_adjacent_pairs_and_fs_smoke():
    rng = np.random.default_rng(0)
    n, d = 4000, 8
    X = rng.normal(size=(n, d))
    # mild drift in y and one feature across chunks
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
        X, y, T, [f"f{j}" for j in range(d)],
        select_k=3, y_quantile=0.7, seed=0, max_pairs=10, min_n=50,
    )
    assert len(rows) >= 2
    assert all(r["y_threshold"] == 0.5 for r in rows)

