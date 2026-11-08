"""Tests for weak-5 nonlinear concept-drift FICO design."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fico_heloc_cfperm import clean_heloc  # noqa: E402
from fico_heloc_weak_cd import (  # noqa: E402
    make_weak_cd_batches,
    nonlinear_concept_score,
    sample_weak_drivers,
    _recovery_metrics,
)


def _toy(n=300, seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    cols = {f"F{j}": rng.normal(size=n) for j in range(8)}
    cols["RiskPerformance"] = np.where(rng.random(n) > 0.5, "Good", "Bad")
    return pd.DataFrame(cols)


def test_nonlinear_score_and_sample():
    Z = np.random.default_rng(0).normal(size=(50, 5))
    g = nonlinear_concept_score(Z)
    assert g.shape == (50,)
    names = [f"F{j}" for j in range(8)]
    d = sample_weak_drivers(names, k=5, seed=1)
    assert len(d) == 5 and len(set(d)) == 5


def test_make_batches_pure_cd():
    df = clean_heloc(_toy())
    X, Y, W, names, meta = make_weak_cd_batches(
        df, k_drivers=5, intensity=2.0, max_n_per_arm=80, seed=3
    )
    assert X.shape[0] == len(Y) == len(W)
    assert set(np.unique(W)) == {0, 1}
    assert len(meta.driver_features) == 5
    # P(X) shared: means of drivers should be close across arms (random split)
    for j in meta.driver_idx:
        m0, m1 = X[W == 0, j].mean(), X[W == 1, j].mean()
        assert abs(m0 - m1) < 1.0  # loose; not a median cut


def test_recovery_metrics():
    m = _recovery_metrics(
        ["a", "b", "c", "d", "e"],
        ["b", "z"],
        ["a", "b", "x", "y", "w"],
        k=5,
    )
    assert m["n_drivers_in_top5"] == 2.0
    assert m["any_driver_rejected"] == 1.0
