"""Smoke tests for OOD-weight online eval harness."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from ood_weight_online_eval import (  # noqa: E402
    attach_regret,
    load_board,
    mahal_sample_weights,
    run_board_suite,
    run_weight_method,
)


def test_surrogate_boards_load():
    for board in ("interstate", "stock", "synthetic"):
        stream = load_board(board, n_per=30, n_batches=4, seed=0)
        assert stream.X.ndim == 2
        assert stream.y.shape[0] == stream.X.shape[0]
        assert stream.batch.max() == 3


def test_mahal_weights_peak_near_last_batch():
    stream = load_board("interstate", n_per=40, n_batches=4, seed=1)
    w = mahal_sample_weights(stream.X, stream.batch, t=3, gamma=1.0)
    assert float(w[stream.batch == 2].mean()) >= float(w[stream.batch == 0].mean()) - 1e-9


def test_regret_zero_for_uniform():
    stream = load_board("stock", n_per=40, n_batches=4, seed=2)
    rec = run_weight_method(stream, "uniform")
    rec = attach_regret(rec, rec["cum_mse"])
    assert abs(rec["regret"]) < 1e-9


def test_quick_suite_keys():
    suite = run_board_suite(
        board="interstate",
        seeds=[0, 1],
        methods=["uniform", "hop", "dga"],
        n_per=40,
        n_batches=4,
    )
    assert "table" in suite and "hop" in suite["table"]
    assert "regret" in suite["table"]["hop"]
    assert suite["table"]["uniform"]["regret"]["mean"] == 0.0 or abs(
        suite["table"]["uniform"]["regret"]["mean"]
    ) < 1e-9
