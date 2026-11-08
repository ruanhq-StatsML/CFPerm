"""Tests for post-hoc subset localization helpers."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from posthoc_subset_localization import (  # noqa: E402
    feature_localize,
    run_posthoc_on_batches,
    subset_on_driver,
    _dist,
)


def test_dist_and_feature_localize():
    rng = np.random.default_rng(0)
    n, p = 200, 6
    X0 = rng.normal(size=(n, p))
    X1 = rng.normal(size=(n, p))
    X1[:, 0] += 1.5  # driver shift
    Y0 = X0[:, 0] + rng.normal(scale=0.3, size=n)
    Y1 = X1[:, 0] + rng.normal(scale=0.3, size=n)
    names = [f"X{j}" for j in range(p)]
    d = _dist(X0[:, 0], X1[:, 0])
    assert d["TVD"] >= 0 and d["Hellinger"] >= 0
    feat, meta, deg = feature_localize(
        X0, Y0, X1, Y1, seed=1, names=names, true_idx=[0], top_k=3, n_perm=5
    )
    assert len(feat) == p
    assert meta["top_j"] == 0 or abs(deg[0]) >= np.median(np.abs(deg))
    sub, contrast, flag = subset_on_driver(
        X0, Y0, X1, Y1, B=4, j_star=int(meta["top_j"]), driver_idx=meta["driver_idx"], seed=2, n_perm=3
    )
    assert set(sub.subset) == {"low", "high"}
    assert flag in ("low", "high", "")


def test_run_posthoc_pipeline():
    rng = np.random.default_rng(2)
    n, p = 160, 5
    X0 = rng.normal(size=(n, p))
    X1 = X0.copy() + rng.normal(scale=0.1, size=(n, p))
    X1[:, 2] += 2.0
    Y0 = rng.normal(size=n)
    Y1 = rng.normal(size=n)
    names = [f"F{j}" for j in range(p)]
    out = run_posthoc_on_batches(
        X0, Y0, X1, Y1, names, split_feature="F2", seed=3, B=3, top_k=3, n_perm=3, true_feature="F2"
    )
    assert "feature_rank" in out and "subset_driver" in out
    assert out["summary"]["split_feature"] == "F2"
