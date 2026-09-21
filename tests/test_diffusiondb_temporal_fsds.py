"""Smoke tests for DiffusionDB temporal FSDS helpers (no Hub download)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.run_diffusiondb_temporal_fsds import (
    assign_time_windows,
    build_light_token_X,
    mean_shift_by_feature,
    rf_domain_vimp,
    run_fsds_y,
)


def test_time_windows_balanced():
    df = pd.DataFrame(
        {
            "timestamp": pd.date_range("2022-08-01", periods=100, freq="h", tz="UTC"),
            "prompt_clean": ["a red car"] * 50 + ["photorealistic portrait"] * 50,
        }
    )
    out = assign_time_windows(df, n_windows=2)
    assert set(out["T"].unique()) == {0, 1}
    assert abs((out["T"] == 0).sum() - (out["T"] == 1).sum()) <= 1


def test_light_tfidf_and_fsds_smoke():
    prompts = [
        "artstation trending high quality fantasy",
        "artstation detailed illustration",
        "photorealistic 8k portrait photo",
        "photorealistic cinematic lighting",
    ] * 40
    X, names, _ = build_light_token_X(prompts, method="tfidf", max_features=64)
    assert X.shape[0] == len(prompts)
    assert X.shape[1] == len(names) and X.shape[1] > 0
    T = np.array([0] * 80 + [1] * 80)
    # make Y correlate with "photo" side in late half a bit
    y = np.array([0] * 60 + [1] * 20 + [0] * 20 + [1] * 60)
    vimp = rf_domain_vimp(X, T, seed=0, n_trees=20)
    assert vimp.shape == (X.shape[1],)
    res = run_fsds_y(X[T == 0], y[T == 0], X[T == 1], y[T == 1], names, select_k=10, seed=0)
    assert res["ok"]
    assert len(res["selected"]) >= 1
    cmean = mean_shift_by_feature(X[T == 0], X[T == 1], names)
    assert "delta" in cmean.columns and "sign" in cmean.columns
