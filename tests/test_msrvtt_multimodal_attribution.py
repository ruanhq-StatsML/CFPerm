"""Unit tests for MSR-VTT multimodal FSDS attribution (no zip required)."""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from msrvtt_multimodal_attribution import (  # noqa: E402
    P_X,
    _mean_sd_ci,
    bundle_from_s,
    drop_group,
    group_label_permutation_test,
    holm_adjust,
    load_window_bundle,
    make_synthetic_bundle,
    modality_mass,
    rf_domain,
    standardize_columns,
)


def test_s_layout_and_early_late():
    b = make_synthetic_bundle(n_videos=4, n_windows=20, seed=1)
    assert b.X.shape == (80, P_X)
    assert b.s.shape == (80, P_X + 1)
    assert set(np.unique(b.W)) == {0, 1}
    # first half of each video is early
    for v in range(4):
        idx = np.flatnonzero(b.video_id == v)
        order = idx[np.argsort(b.window_idx[idx])]
        assert (b.W[order[:10]] == 0).all()
        assert (b.W[order[10:]] == 1).all()


def test_zip_roundtrip_s(tmp_path):
    b = make_synthetic_bundle(n_videos=3, n_windows=16, seed=2)
    np.save(tmp_path / "s.npy", b.s)
    np.save(tmp_path / "video_labels.npy", b.video_id)
    np.save(tmp_path / "window_idx.npy", b.window_idx)
    z = tmp_path / "feature_video_audio.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.write(tmp_path / "s.npy", "s.npy")
        zf.write(tmp_path / "video_labels.npy", "video_labels.npy")
        zf.write(tmp_path / "window_idx.npy", "window_idx.npy")
    got = load_window_bundle(z, root=tmp_path / "extract")
    assert got.X.shape == b.X.shape
    assert np.allclose(got.X, b.X)
    assert (got.W == b.W).all()


def test_group_mmd_block_recovers_video_shift():
    from msrvtt_multimodal_attribution import group_mmd_loco, _share_from_contrib

    b = make_synthetic_bundle(
        n_videos=4, n_windows=20, seed=3, video_shift=1.6, audio_shift=0.05, text_shift=0.02
    )
    X = standardize_columns(b.X)
    _, loco, block = group_mmd_loco(X[b.W == 0], X[b.W == 1], max_n=40, seed=3)
    share = _share_from_contrib(block)
    assert share["video"] > share["audio"]
    assert share["video"] > share["text"]


def test_rf_recovers_video_block_shift():
    b = make_synthetic_bundle(
        n_videos=6,
        n_windows=24,
        seed=7,
        video_shift=1.8,
        audio_shift=0.05,
        text_shift=0.02,
    )
    vimp, auc = rf_domain(b.X[b.W == 0], b.X[b.W == 1], seed=7, n_estimators=40)
    _, share = modality_mass(vimp)
    assert share["video"] > share["audio"]
    assert share["video"] > share["text"]
    assert share["video"] > 0.45
    assert auc > 0.75
    perm = group_label_permutation_test(vimp, B=79, seed=7)
    assert perm["p"] < 0.05


def test_drop_group_shapes():
    X = np.zeros((5, P_X))
    assert drop_group(X, "video").shape[1] == P_X - 768
    assert drop_group(X, "audio").shape[1] == P_X - 512
    assert drop_group(X, "text").shape[1] == P_X - 768


def test_holm_adjust():
    adj = holm_adjust([0.001, 0.04, 0.02])
    assert adj[0] <= adj[2] <= adj[1] or adj[0] < 0.01
    assert np.all(adj >= np.array([0.001, 0.04, 0.02]) - 1e-12)
    assert np.all(adj <= 1.0)


def test_bootstrap_mean_sd_from_ten_replicates():
    rng = np.random.default_rng(0)
    samples = rng.normal(loc=0.4, scale=0.05, size=10)
    rec = _mean_sd_ci(samples)
    assert rec["mean"] == pytest.approx(float(np.mean(samples)))
    assert rec["sd"] == pytest.approx(float(np.std(samples, ddof=1)))
    assert rec["var"] == pytest.approx(float(np.var(samples, ddof=1)))
    half = 1.96 * rec["sd"]
    assert rec["ci95"][0] == pytest.approx(rec["mean"] - half)
    assert rec["ci95"][1] == pytest.approx(rec["mean"] + half)


def test_standardize_columns_unit_variance():
    rng = np.random.default_rng(1)
    X = rng.normal(loc=3.0, scale=4.0, size=(40, 6))
    Z = standardize_columns(X)
    assert np.allclose(Z.mean(0), 0.0, atol=1e-10)
    assert np.allclose(Z.std(0), 1.0, atol=1e-10)


def test_bundle_from_s_uses_label_as_video_id():
    rng = np.random.default_rng(0)
    n_v, n_w = 5, 10
    X = rng.normal(size=(n_v * n_w, P_X))
    y = np.repeat(np.arange(n_v), n_w).astype(float)
    s = np.hstack([X, y[:, None]])
    b = bundle_from_s(s)
    assert len(np.unique(b.video_id)) == n_v
    assert b.X.shape[1] == P_X
