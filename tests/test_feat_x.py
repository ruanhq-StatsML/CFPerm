"""Tests for edge∥stream observation X and ToT attachment."""
from __future__ import annotations

import numpy as np
import pytest

from agod.feat_x import (
    align_rows,
    concat_edge_stream,
    split_edge_stream,
    tot_observation_spec,
)
from agod.tot_eff import tot_with_observation


def test_concat_edge_stream_layout_and_names():
    rng = np.random.default_rng(0)
    Xe = rng.normal(size=(12, 3))
    Xs = rng.normal(size=(12, 2))
    y = rng.normal(size=12)
    blocks = concat_edge_stream(
        Xe,
        Xs,
        y=y,
        edge_names=["e_ctr", "u_ctr", "i_ctr"],
        stream_names=["temp", "hour"],
    )
    assert blocks.X.shape == (12, 5)
    assert blocks.edge_dim == 3 and blocks.stream_dim == 2
    assert blocks.names == [
        "edge:e_ctr",
        "edge:u_ctr",
        "edge:i_ctr",
        "stream:temp",
        "stream:hour",
    ]
    assert blocks.to_meta()["layout"] == "edge||stream"
    assert blocks.to_meta()["intermediate_is"] == "policy_thought_not_yhat"
    e2, s2 = split_edge_stream(blocks)
    np.testing.assert_allclose(e2, Xe)
    np.testing.assert_allclose(s2, Xs)


def test_align_rows_truncates_to_min():
    Xe = np.ones((10, 2))
    Xs = np.ones((7, 1))
    y = np.arange(9)
    xe, xs, yy = align_rows(Xe, Xs, y=y)
    assert xe.shape[0] == xs.shape[0] == yy.shape[0] == 7


def test_concat_rejects_name_dim_mismatch():
    with pytest.raises(ValueError):
        concat_edge_stream(
            np.zeros((4, 2)),
            np.zeros((4, 1)),
            edge_names=["only_one"],
        )


def test_tot_with_observation_attaches_contract():
    blocks = concat_edge_stream(
        np.zeros((20, 2)),
        np.ones((20, 3)),
        y=np.zeros(20),
    )
    rep = tot_with_observation(
        blocks,
        gate_duty=0.2,
        rank_effs={"ref": 0.0, "probe": 0.1, "refit": 1.0},
        mse_effs={"ref": 0.0, "probe": -0.1, "refit": -0.05},
        rel_sqrt=1.1,
        rel_cbrt=1.05,
    )
    obs = rep["observation"]
    assert obs["X"].startswith("edge_features")
    assert "policy Thought" in obs["intermediate"]
    assert obs["edge_dim"] == 2 and obs["stream_dim"] == 3
    assert obs["n"] == 20
    assert rep["best"] is not None
    # contract helper stays stable
    spec = tot_observation_spec()
    assert "||" in spec["X"]
