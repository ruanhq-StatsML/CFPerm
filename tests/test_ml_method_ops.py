"""Tests for ML method operators (compose with transfer_null; no core rewrite)."""
from __future__ import annotations

import numpy as np

from agod.ml_method_ops import (
    excess_learning_curve,
    ml_method_suite,
    page_hinkley_skill,
    partial_excess_auc,
    residualize_vs_confounder,
)


def test_residualize_recovers_linear_confounder():
    rng = np.random.default_rng(0)
    c = rng.normal(size=400)
    s = 2.0 * c + rng.normal(scale=0.1, size=400)
    out = residualize_vs_confounder(s, c)
    assert out["r2"] > 0.9
    assert abs(out["beta"] - 2.0) < 0.2
    assert float(np.std(out["residual"])) < float(np.std(s))


def test_partial_excess_drops_when_score_is_volume():
    rng = np.random.default_rng(1)
    n = 600
    # rare-ish labels driven by volume confounder
    c = rng.integers(0, 50, size=n).astype(float)
    logits = 0.15 * c - 3.0
    p = 1.0 / (1.0 + np.exp(-logits))
    y = (rng.random(n) < p).astype(int)
    # probe score ≈ volume (the confounder)
    score = c + rng.normal(scale=0.5, size=n)
    if len(np.unique(y)) < 2:
        y[0], y[1] = 0, 1
    rep = partial_excess_auc(y, score, c, n_perm=10, seed=1)
    assert np.isfinite(rep["excess_raw"])
    # after residualizing volume, excess should fall
    assert rep["excess_partial"] < rep["excess_raw"] - 0.02
    assert rep["delta_excess"] > 0.02
    assert rep["confounder_r2"] > 0.5


def test_learning_curve_smoke():
    rng = np.random.default_rng(2)
    n = 300
    y = rng.integers(0, 2, size=n)
    y[0], y[1] = 0, 1
    score = y.astype(float) + rng.normal(scale=0.3, size=n)
    curve = excess_learning_curve(y, score, n_perm=4, seed=2)
    assert len(curve["points"]) == 5
    assert "reading" in curve


def test_page_hinkley_alarms_on_skill_drop():
    # stable high excess then sustained drop
    series = [0.20] * 12 + [0.02] * 12
    ph = page_hinkley_skill(series, delta=0.002, lambda_thresh=0.03)
    assert ph["alarm"] is True
    assert ph["alarm_index"] is not None
    assert ph["alarm_index"] >= 10


def test_page_hinkley_quiet_when_stable():
    ph = page_hinkley_skill([0.12, 0.11, 0.13, 0.12, 0.11, 0.12], lambda_thresh=0.2)
    assert ph["alarm"] is False


def test_ml_method_suite_composes():
    rng = np.random.default_rng(3)
    n = 250
    c = rng.normal(size=n)
    y = (rng.random(n) < 0.4).astype(int)
    y[0], y[1] = 0, 1
    score = 0.5 * c + y + rng.normal(scale=0.5, size=n)
    suite = ml_method_suite(
        y, score, confounder=c, excess_history=[0.1, 0.11, 0.09], n_perm=4, seed=3
    )
    assert "excess_auc" in suite
    assert suite["partial_excess"] is not None
    assert "learning_curve" in suite
    assert "page_hinkley" in suite
    assert "does not modify" in suite["method_note"]
