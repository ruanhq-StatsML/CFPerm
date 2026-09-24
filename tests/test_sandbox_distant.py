"""Tests for sandbox distant jobs (no AGOD method coupling)."""
from __future__ import annotations

import numpy as np

from sandbox.forecast_bakeoff import Pack, make_supervised, run_bakeoff
from sandbox.prompt_themes import theme_clusters


def test_make_supervised_shapes():
    X = np.random.randn(100, 3)
    y = np.random.randn(100)
    Z, t = make_supervised(X, y, n_lags=4)
    assert Z.shape[0] == 96
    assert Z.shape[1] == 4 + 3
    assert t.shape == (96,)


def test_bakeoff_smoke_synthetic():
    rng = np.random.default_rng(0)
    n = 400
    X = rng.normal(size=(n, 2))
    # AR(1)-ish
    y = np.zeros(n)
    for i in range(1, n):
        y[i] = 0.7 * y[i - 1] + 0.1 * X[i, 0] + rng.normal(scale=0.3)
    pack = Pack(name="toy_ar1", X=X, y=y)
    card = run_bakeoff(pack, n_lags=3, seed=0)
    assert card["ok"]
    assert card["best_rmse"] in ("naive_last", "ridge", "hgb")
    assert card["models"]["hgb"]["rmse"] < card["models"]["naive_last"]["rmse"] * 1.2


def test_theme_clusters_smoke():
    prompts = [
        f"a red car on the highway under sunset lighting variant {i}"
        for i in range(40)
    ] + [
        f"portrait of a woman with blue eyes studio photo {i}" for i in range(40)
    ] + [
        f"fantasy dragon flying over mountains digital art {i}" for i in range(40)
    ]
    rep = theme_clusters(prompts, k=3, max_features=500, seed=0)
    assert rep["ok"]
    assert len(rep["top_terms"]) == 3
    assert rep["n_docs"] == 120


def test_opportunity_map_and_flywheel_smoke():
    from sandbox.forecast_bakeoff import Pack, run_bakeoff
    from sandbox.forecast_flywheel import opportunity_map_from_bakeoff, run_flywheel

    rng = np.random.default_rng(0)
    n = 500
    X = rng.normal(size=(n, 2))
    y = np.zeros(n)
    for i in range(1, n):
        y[i] = 0.6 * y[i - 1] + 0.2 * X[i, 0] + rng.normal(scale=0.25)
    pack = Pack(name="toy", X=X, y=y)
    card = run_bakeoff(pack, n_lags=3, seed=0)
    bakeoff = {"cards": [card]}
    opps = opportunity_map_from_bakeoff(bakeoff)
    assert any(o["pack"] == "toy" for o in opps)
    fw = run_flywheel(pack, model_name="hgb", warm=120, max_steps=80, seed=0)
    assert fw["ok"]
    assert fw["n_steps"] == 80
    assert "surprise_rate" in fw
