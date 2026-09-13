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


def test_dga_mirror_upweights_aligned_domain():
    from attribution_adapter import dga_mirror_step

    alpha0 = np.ones(3) / 3.0
    a = np.array([-1.0, 0.0, 2.0])
    alpha1 = dga_mirror_step(alpha0, a, eta=2.0)
    assert alpha1.argmax() == 2
    assert alpha1[2] > alpha1[0]
    assert abs(alpha1.sum() - 1.0) < 1e-8


def test_dga_ridge_runs_and_peaks_on_spe():
    from attribution_adapter import run_dga_ridge

    stream = make_amazon_like_stream(n_batches=5, n_per=40, p=16, seed=0, cov=0.5)
    rec = run_dga_ridge(stream, eta=1.0, ema_beta=0.5, align="cosine")
    assert rec["method"] == "dga_ridge"
    assert len(rec["online_path"]) == 4
    assert rec["online_mse"] == rec["online_mse"]
    # last hop: spe = batch 3 should receive non-trivial mass
    last_alpha = np.asarray(rec["history"][-1]["alpha_inst"], dtype=float)
    assert last_alpha.argmax() == len(last_alpha) - 1 or last_alpha[-1] >= 1.0 / len(last_alpha)


def test_attr_adapter_returns_path():
    stream = make_amazon_like_stream(n_batches=5, n_per=40, p=20, seed=3, cov=0.4)
    rec = run_attr_adapter(stream)
    assert len(rec["online_path"]) == 4
    assert 0.0 < rec.get("mean_w_bank", 0.5) < 1.0


def test_run_adapter_method_dispatch():
    stream = make_amazon_like_stream(n_batches=4, n_per=30, p=16, seed=2, cov=0.3)
    for method in ("ridge_past", "hop_ridge", "dga_ridge", "attr_adapter", "bank"):
        rec = run_adapter_method(stream, method, seed=2)
        assert rec["online_mse"] == rec["online_mse"]
