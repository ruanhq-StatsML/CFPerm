"""Smoke tests for DiffusionDB temporal FS helpers (no Hub download)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.run_diffusiondb_temporal_fsds import (
    assign_time_windows,
    build_light_token_X,
    concat_hyperparams,
    mean_shift_by_feature,
    rf_domain_vimp,
    run_fsds_y,
)


def test_time_windows_equal_count_balanced():
    df = pd.DataFrame(
        {
            "timestamp": pd.date_range("2022-08-01", periods=100, freq="h", tz="UTC"),
            "prompt_clean": ["a red car"] * 50 + ["photorealistic portrait"] * 50,
        }
    )
    out = assign_time_windows(df, n_windows=2, scheme="equal_count")
    assert set(out["T"].unique()) == {0, 1}
    assert abs((out["T"] == 0).sum() - (out["T"] == 1).sum()) <= 1
    assert out.attrs["window_meta"]["scheme"] == "equal_count"


def test_time_windows_equal_time_and_width_differ():
    # denser early, sparser late → equal_time imbalances n
    ts_early = pd.date_range("2022-08-01", periods=80, freq="h", tz="UTC")
    ts_late = pd.date_range("2022-08-10", periods=20, freq="12h", tz="UTC")
    df = pd.DataFrame(
        {
            "timestamp": ts_early.append(ts_late),
            "prompt_clean": ["x"] * 100,
        }
    )
    et = assign_time_windows(df, n_windows=2, scheme="equal_time")
    ec = assign_time_windows(df, n_windows=2, scheme="equal_count")
    wd = assign_time_windows(df, n_windows=2, scheme="width", width_hours=48.0)
    assert et.attrs["window_meta"]["scheme"] == "equal_time"
    assert wd.attrs["window_meta"]["scheme"] == "width"
    # schemes disagree on counts or binning
    assert not (
        (et["T"].to_numpy() == ec["T"].to_numpy()).all()
        and (ec["T"].to_numpy() == wd["T"].to_numpy()).all()
    )


def test_concat_hyperparams_extends_X():
    prompts = ["artstation fantasy"] * 30 + ["photorealistic portrait"] * 30
    X, names, _ = build_light_token_X(prompts, method="tfidf", max_features=32)
    df = pd.DataFrame(
        {
            "cfg": [7.0] * 30 + [12.0] * 30,
            "step": [20] * 30 + [50] * 30,
            "sampler": ["k_euler"] * 30 + ["k_euler_ancestral"] * 30,
        }
    )
    Xp, names_p = concat_hyperparams(X, names, df)
    assert Xp.shape[1] > X.shape[1]
    assert "hp_cfg" in names_p and "hp_step" in names_p
    assert any(n.startswith("sampler=") for n in names_p)


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
    y = np.array([0] * 60 + [1] * 20 + [0] * 20 + [1] * 60)
    vimp = rf_domain_vimp(X, T, seed=0, n_trees=20)
    assert vimp.shape == (X.shape[1],)
    res = run_fsds_y(
        X[T == 0], y[T == 0], X[T == 1], y[T == 1], names, select_k=10, seed=0
    )
    assert res["ok"]
    assert len(res["selected"]) >= 1
    cmean = mean_shift_by_feature(X[T == 0], X[T == 1], names)
    assert "delta" in cmean.columns and "sign" in cmean.columns
