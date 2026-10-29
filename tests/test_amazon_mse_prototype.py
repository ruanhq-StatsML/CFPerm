"""Amazon MSE prototype: shapes, SDC stars, GPM-Ridge (no Hub required)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from amazon_mse_prototype import (  # noqa: E402
    describe_shapes,
    run_bank_mse,
    run_mse_method,
    run_sdc_mse,
    sdc_predict,
    star_prototypes,
)
from amazon_continuous_batches import make_amazon_like_stream  # noqa: E402
from msrvtt_multimodal_attribution import P_X  # noqa: E402


def test_describe_shapes_lists_both_boards():
    stream = make_amazon_like_stream(n_batches=4, n_per=20, p=16, seed=1)
    rec = describe_shapes(stream)
    assert rec["msrvtt"]["X"] == "(n, %d)" % P_X
    live = rec["amazon_live"]
    assert live["X"] == (80, 16)
    assert live["K"] == 4
    assert live["P"] == 16
    assert live["n_per"][0] == 20


def test_sdc_predict_stays_in_star_range():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(30, 8))
    y = np.array([1, 2, 3, 4, 5] * 6, dtype=float)
    mu, _ = star_prototypes(X, y)
    pred = sdc_predict(X, mu)
    assert pred.shape == (30,)
    assert pred.min() >= 1.0 and pred.max() <= 5.0


def test_bank_mse_beats_star_sdc_on_synthetic():
    stream = make_amazon_like_stream(n_batches=6, n_per=40, seed=2, cov=0.4)
    bank = run_bank_mse(stream)
    sdc = run_sdc_mse(stream)
    assert bank["online_mse"] == bank["online_mse"]
    assert bank["online_mse"] < sdc["online_mse"]


def test_sdc_mse_finite_on_synthetic():
    stream = make_amazon_like_stream(n_batches=6, n_per=40, seed=2, cov=0.4)
    rec = run_sdc_mse(stream)
    assert rec["online_mse"] == rec["online_mse"]
    assert rec["online_mse"] < 20.0
    assert len(rec["online_path"]) == 5


def test_gpm_and_plateau_return_mse():
    stream = make_amazon_like_stream(n_batches=6, n_per=36, seed=4, cov=0.35)
    plat = run_mse_method(stream, "plateau", seed=4, steps_per_batch=2, warmup_steps=2)
    gpm = run_mse_method(stream, "gpm_typed", seed=4, steps_per_batch=2, warmup_steps=2)
    assert plat["online_mse"] == plat["online_mse"]
    assert gpm["online_mse"] == gpm["online_mse"]
    assert "gate" in gpm["history"][0]
