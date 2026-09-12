"""Amazon continuous-batch TSS: identification, MSE, and regret (no Hub required)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from amazon_continuous_batches import (  # noqa: E402
    concept_intensity_mse,
    covariate_intensity_text,
    featurize_reviews,
    hindsight_constant,
    make_amazon_like_stream,
    parse_jsonl_records,
    run_amazon_method,
    tss_eta_scalar,
)


def test_parse_jsonl_skips_truncated_and_empty():
    raw = (
        '{"rating": 5.0, "title": "ok", "text": "great gift card"}\n'
        '{"rating": 1.0, "title": "bad", "text":\n'
        '{"rating": 4, "title": "", "text": "works"}\n'
    )
    rows = parse_jsonl_records(raw)
    assert len(rows) == 2
    assert rows[0]["y"] == 5.0
    assert "gift" in rows[0]["text"]


def test_covariate_intensity_grows_with_mean_shift():
    quiet = make_amazon_like_stream(n_batches=6, n_per=60, seed=3, cov=0.0)
    shifted = make_amazon_like_stream(n_batches=6, n_per=60, seed=3, cov=0.6)
    i0 = np.flatnonzero(quiet.batch == 0)
    i5 = np.flatnonzero(quiet.batch == 5)
    c0 = covariate_intensity_text(quiet.X[i0], quiet.X[i5])
    c1 = covariate_intensity_text(shifted.X[i0], shifted.X[i5])
    assert c1 > c0 + 0.5


def test_concept_intensity_jumps_after_rating_flip():
    stream = make_amazon_like_stream(
        n_batches=8, n_per=90, seed=5, cov=0.0, concept=1.0, concept_at=4
    )
    i3 = np.flatnonzero(stream.batch == 3)
    i4 = np.flatnonzero(stream.batch == 4)
    pre = np.flatnonzero(stream.batch == 2)
    d_pre, _ = concept_intensity_mse(stream.X[pre], stream.y[pre], stream.X[i3], stream.y[i3])
    d_post, extras = concept_intensity_mse(stream.X[i3], stream.y[i3], stream.X[i4], stream.y[i4])
    assert d_post > d_pre + 0.05
    assert extras["mse_al"] > extras["mse0"]


def test_tss_eta_scalar_signs_and_quiet_hold():
    held = 0.08
    quiet = tss_eta_scalar(0.02, 0.01, eta0=0.10, prev=held, n_iter=1)
    assert quiet == held
    cov = tss_eta_scalar(1.2, 0.0, eta0=0.10, prev=held, n_iter=1)
    concept = tss_eta_scalar(0.0, 0.8, eta0=0.10, prev=held, n_iter=1)
    assert concept > cov
    assert cov < held
    assert concept > held


def test_tss_shrinks_eta_under_covariate_only():
    stream = make_amazon_like_stream(n_batches=8, n_per=50, seed=9, cov=0.45, concept=0.0)
    tss, _ = run_amazon_method(stream, method="tss", eta0=0.10, steps_per_batch=3, seed=9)
    const, _ = run_amazon_method(stream, method="constant", eta0=0.10, steps_per_batch=3, seed=9)
    assert tss["mean_c"] > 0.2
    assert tss["mean_eta"] < const["mean_eta"]


def test_featurize_keeps_nonnegative_tfidf():
    texts = ["great gift card five stars"] * 24 + ["noisy guitar amp strings"] * 24
    payload = {
        "texts": texts,
        "y": np.r_[np.full(24, 5.0), np.full(24, 3.0)],
        "batch": np.r_[np.zeros(24, dtype=int), np.ones(24, dtype=int)],
        "categories": ("gift", "music"),
        "meta": {},
    }
    stream = featurize_reviews(payload, max_features=16, min_df=1)
    assert stream.X.shape[0] == 48
    assert stream.X.min() >= -1e-12


def test_hindsight_beats_other_constant_etas():
    stream = make_amazon_like_stream(n_batches=6, n_per=40, seed=4, cov=0.2)
    hind = hindsight_constant(stream, etas=(0.04, 0.10, 0.20), steps_per_batch=2, seed=4)
    rec, _ = run_amazon_method(
        stream, method="constant", eta0=hind["hindsight_eta"], steps_per_batch=2, seed=4
    )
    assert rec["cum_mse"] == hind["cum_mse"]
    const, _ = run_amazon_method(stream, method="constant", eta0=0.10, steps_per_batch=2, seed=4)
    assert hind["cum_mse"] <= const["cum_mse"] + 1e-12
    assert hind["hindsight_eta"] in (0.04, 0.10, 0.20)
