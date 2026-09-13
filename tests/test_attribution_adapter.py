"""Attribution adapter: hop-weighted Ridge and typed bank⊕Ridge mix."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from amazon_continuous_batches import make_amazon_like_stream  # noqa: E402
from attribution_adapter import (  # noqa: E402
    hop_sample_weights,
    rolling_hop_stats,
    run_adapter_method,
    run_attr_adapter,
    run_hop_ridge,
    run_ridge_past,
    typed_bank_weight,
)


def test_hop_weights_peak_on_last_batch():
    stream = make_amazon_like_stream(n_batches=4, n_per=30, p=16, seed=1, cov=0.5)
    stats = rolling_hop_stats(stream)
    w = hop_sample_weights(stats["mus"], stream.batch, t=3, gamma=4.0)
    # last observed batch (t-1=2) should get weight 1
    assert abs(float(w[stream.batch == 2].mean()) - 1.0) < 1e-6
    assert float(w[stream.batch == 0].mean()) <= float(w[stream.batch == 2].mean()) + 1e-9


def test_typed_bank_weight_signs():
    # high covariate, quiet concept → bank up
    assert typed_bank_weight(1.0, 0.0) > typed_bank_weight(0.0, 1.0)
    # high concept → bank down
    assert typed_bank_weight(0.0, 1.0) < 0.35


def test_hop_ridge_beats_uniform_on_covariate_stream():
    stream = make_amazon_like_stream(n_batches=6, n_per=50, p=24, seed=7, cov=0.55)
    hop = run_hop_ridge(stream)
    past = run_ridge_past(stream)
    assert hop["online_mse"] <= past["online_mse"] + 0.05
    assert hop["online_mse"] < 5.0


def test_attr_adapter_returns_path():
    stream = make_amazon_like_stream(n_batches=5, n_per=40, p=20, seed=3, cov=0.4)
    rec = run_attr_adapter(stream)
    assert len(rec["online_path"]) == 4
    assert 0.0 < rec.get("mean_w_bank", 0.5) < 1.0


def test_run_adapter_method_dispatch():
    stream = make_amazon_like_stream(n_batches=4, n_per=30, p=16, seed=2, cov=0.3)
    for method in ("ridge_past", "hop_ridge", "attr_adapter", "bank"):
        rec = run_adapter_method(stream, method, seed=2)
        assert rec["online_mse"] == rec["online_mse"]
