"""Unit tests for FICO HELOC CFPerm helpers (synthetic, no network)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from fico_heloc_cfperm import (  # noqa: E402
    cfperm_test,
    clean_heloc,
    median_split_batches,
    power_summary,
    subset_localize,
    top_lr_features,
    vimp_grf,
    vimp_loco,
    vimp_permucate,
)


def _toy_heloc(n=400, seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    x0 = rng.normal(size=n)
    x1 = rng.normal(size=n)
    x2 = rng.normal(size=n)
    noise = rng.normal(size=n)
    logits = 1.2 * x0 + 0.8 * x1 + 0.2 * x2 + 0.3 * noise
    y = (logits > np.median(logits)).astype(int)
    # inject a few -9 rows that clean_heloc should drop
    x0 = x0.copy()
    x0[:3] = -9
    return pd.DataFrame(
        {
            "RiskPerformance": np.where(y == 1, "Good", "Bad"),
            "ExternalRiskEstimate": x0,
            "MSinceOldestTradeOpen": x1,
            "AverageMInFile": x2,
            "NumSatisfactoryTrades": rng.normal(size=n),
            "NumInqLast6M": rng.normal(size=n),
        }
    )


def test_clean_and_top_features():
    df = _toy_heloc()
    clean = clean_heloc(df)
    assert len(clean) == len(df) - 3
    tops = top_lr_features(clean, top_k=2)
    assert len(tops) == 2
    assert "ExternalRiskEstimate" in tops or "MSinceOldestTradeOpen" in tops


def test_median_split_and_vimps():
    df = clean_heloc(_toy_heloc(n=300, seed=1))
    feat = "ExternalRiskEstimate"
    X, Y, W, names = median_split_batches(df, feat, max_n_per_arm=80, seed=2)
    assert X.shape[0] == len(Y) == len(W)
    assert set(np.unique(W)) == {0, 1}
    loco = vimp_loco(X, Y, W, seed=0, n_estimators=20, n_splits=2)
    assert loco.shape == (X.shape[1],)
    assert np.all(np.isfinite(loco))
    grf = vimp_grf(X, Y, W, seed=0, n_estimators=20)
    assert grf.shape == (X.shape[1],)
    permu = vimp_permucate(X, Y, W, seed=0, n_estimators=20, n_perm_cate=3)
    assert permu.shape == (X.shape[1],)


def test_cfperm_and_subset_power():
    df = clean_heloc(_toy_heloc(n=320, seed=3))
    feat = top_lr_features(df, top_k=1)[0]
    X, Y, W, names = median_split_batches(df, feat, max_n_per_arm=90, seed=4)
    res = cfperm_test(
        X,
        Y,
        W,
        names,
        vimp="loco",
        split_feature=feat,
        n_perm=4,
        n_estimators=20,
        n_perm_cate=2,
        level_feature=0.9,
        level_across_feature=0.8,
        seed=0,
        n_jobs=1,
    )
    assert res.n == len(Y)
    assert len(res.imp) == len(names)
    rows = subset_localize(
        X, Y, W, names, res.top_features[:2], split_feature=feat, vimp="loco", n_bins=3
    )
    assert isinstance(rows, list)
    summary = power_summary([res])
    assert "overall" in summary and "by_vimp" in summary
