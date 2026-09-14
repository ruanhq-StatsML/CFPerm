"""Real consecutive-batch streams (no shuffle, no K-fold)."""
from __future__ import annotations

import numpy as np

from agod.po_refit import stream_from_xy
from agod.real_data import (
    ROOT,
    iter_real_streams,
    load_diabetes_readmit,
    load_interstate,
)


def test_stream_from_xy_keeps_row_order():
    X = np.arange(30, dtype=float).reshape(30, 1)
    y = np.arange(30, dtype=float)
    st = stream_from_xy(X, y, n_per=5, n_batches=6, name="ord", task="mse")
    assert st.batch.tolist()[:6] == [0, 0, 0, 0, 0, 1]
    assert np.array_equal(st.y, y)
    assert np.array_equal(st.X, X)


def test_interstate_is_time_ordered_hourly():
    packed = load_interstate(ROOT, max_n=400, seed=0, pca_d=8)
    assert packed is not None
    X, y, task = packed
    assert task == "mse"
    assert len(y) == 400
    assert X.shape[1] >= 8
    assert np.isfinite(X).all() and np.isfinite(y).all()


def test_diabetes_readmit_concatenates_source_then_target():
    packed = load_diabetes_readmit(ROOT, max_n=300, seed=0, pca_d=8)
    assert packed is not None
    X, y, task = packed
    assert task == "acc"
    assert set(np.unique(y).tolist()) <= {0, 1}
    assert len(y) == 300


def test_iter_real_streams_includes_interstate_and_taxi():
    names = []
    for st in iter_real_streams(n_per=50, n_batches=6, max_n=400, seed=0):
        names.append(st.name)
        assert st.batch.min() == 0
        assert st.X.shape[0] == st.y.shape[0]
        assert st.X.shape[0] >= 6 * 50
    assert "interstate" in names
    assert "nyc_taxi" in names
    assert "diabetes_readmit" in names
    assert "california" in names
