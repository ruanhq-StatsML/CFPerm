"""Unit tests for modality-gap clever covariates (synthetic DGP, no HF)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from clever_covariate_gap import (  # noqa: E402
    clever_design,
    compare_raw_vs_clever,
    decompose_modality_gap,
    inject_block_shift,
    make_synthetic_shift,
    relu_normalize,
)


def test_relu_normalize_sums_to_one():
    pi = relu_normalize(np.array([0.2, -1.0, 0.6, 0.0]))
    assert pi.shape == (4,)
    assert np.isclose(pi.sum(), 1.0)
    assert pi[1] == 0.0


def test_synthetic_valence_is_top_contributor():
    X, W, Y, spec = make_synthetic_shift(n=420, d_text=16, d_vad=4, gt="valence", mean_shift=1.4, seed=7)
    gap = decompose_modality_gap(X, W, spec, seed=7, n_splits=3, n_estimators=40)
    top = spec.names[int(np.argmax(gap.pi_consensus))]
    assert top == "valence"
    assert gap.pi_consensus[spec.names.index("valence")] >= 0.35
    assert gap.instance_share.shape == (len(W), 4)
    assert np.allclose(gap.instance_share.sum(axis=1), 1.0, atol=1e-5)
    # clever H uses W; Z does not leak a copy of W as a column of X
    Z = clever_design(X, gap)
    assert Z.shape == (len(W), X.shape[1] + 8)
    # instance shares are functions of X (finite, in [0,1])
    assert np.all(np.isfinite(gap.instance_share))
    assert np.all(gap.instance_share >= -1e-9)


def test_clever_improves_or_matches_domain_auc_on_synthetic():
    X, W, Y, spec = make_synthetic_shift(n=520, d_text=28, d_vad=3, gt="valence", mean_shift=0.85, seed=11)
    gap = decompose_modality_gap(X, W, spec, seed=11, n_splits=3, n_estimators=45)
    row = compare_raw_vs_clever(X, W, Y, spec, gap, seed=11, gt="valence")
    # high-d text noise: clever Z should not hurt, and should recover GT
    assert row["domain_auc_clever"] + 1e-6 >= row["domain_auc_raw"] - 0.03
    assert row["pi_consensus_on_gt"] >= 0.3
    assert row["z_share_on_gt"] >= 0.25


def test_inject_block_recovers_arousal():
    X, W, Y, spec = make_synthetic_shift(n=400, d_text=12, d_vad=4, gt="text", mean_shift=0.15, seed=3)
    X2 = inject_block_shift(X, W, spec, "arousal", alpha=1.3)
    gap = decompose_modality_gap(X2, W, spec, seed=3, n_splits=3, n_estimators=40)
    top = spec.names[int(np.argmax(gap.pi_consensus))]
    assert top == "arousal"


def test_h_clever_shape_and_scale():
    X, W, Y, spec = make_synthetic_shift(n=200, seed=0)
    gap = decompose_modality_gap(X, W, spec, seed=0, n_splits=2, n_estimators=25)
    assert gap.H_clever.shape == (200, 4)
    assert np.all(np.isfinite(gap.H_clever))
