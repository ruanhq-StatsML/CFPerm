"""Unit tests for modality-gap clever covariates (synthetic DGP, no HF)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from sklearn.metrics import roc_auc_score  # noqa: E402

from clever_covariate_gap import (  # noqa: E402
    adaptive_group_scale,
    block_aware_importance,
    block_column_probs,
    clever_design,
    clever_z,
    compare_raw_vs_clever,
    cv_block_aware_auc,
    decompose_modality_gap,
    fit_block_aware_forest,
    inject_block_shift,
    make_synthetic_shift,
    opinion_pool_scores,
    pi_column_replicates,
    relu_normalize,
    smoothed_pi,
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
    Z = clever_z(gap)
    assert Z.shape == (len(W), 4)
    XZ = clever_design(X, gap)
    assert XZ.shape == (len(W), X.shape[1] + 4)
    # instance shares are functions of X (finite, in [0,1])
    assert np.all(np.isfinite(gap.instance_share))
    assert np.all(gap.instance_share >= -1e-9)


def test_clever_z_beats_raw_in_small_n_high_p():
    X, W, Y, spec = make_synthetic_shift(n=160, d_text=40, d_vad=3, gt="valence", mean_shift=1.0, seed=11)
    gap = decompose_modality_gap(X, W, spec, seed=11, n_splits=3, n_estimators=35)
    row = compare_raw_vs_clever(X, W, Y, spec, gap, seed=11, gt="valence", with_po=False, with_block_train=False)
    assert row["p_clever"] == 4
    assert row["p_raw"] > row["p_clever"]
    assert row["domain_auc_clever"] >= 0.75
    assert row["pi_consensus_on_gt"] >= max(row["mass_on_gt_raw"] - 0.08, 0.4)
    assert row["z_logit_vimp"]["valence"] >= 0.35


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


def test_block_column_probs_puts_mass_on_high_pi_block():
    X, W, Y, spec = make_synthetic_shift(n=80, d_text=20, d_vad=4, seed=0)
    pi = np.array([0.05, 0.80, 0.10, 0.05])
    p = block_column_probs(spec, pi, n_features=X.shape[1], floor=0.0)
    assert np.isclose(p.sum(), 1.0)
    val_mass = float(p[spec.slices[spec.names.index("valence")]].sum())
    text_mass = float(p[spec.slices[spec.names.index("text")]].sum())
    assert val_mass > 0.75
    assert val_mass > text_mass
    sm = smoothed_pi(np.array([1.0, 0.0, 0.0, 0.0]), floor=0.05)
    assert sm.min() > 0
    assert np.isclose(sm.sum(), 1.0)


def test_opinion_pool_tracks_gt_block():
    X, W, Y, spec = make_synthetic_shift(n=360, d_text=16, d_vad=4, gt="valence", mean_shift=1.3, seed=5)
    gap = decompose_modality_gap(X, W, spec, seed=5, n_splits=3, n_estimators=40, light=True)
    scores = opinion_pool_scores(gap)
    auc = float(roc_auc_score(W, scores))
    assert scores.shape == (len(W),)
    assert auc >= 0.78
    # pool is a π-mixture of already-fit block RFs — no second-stage RF
    row = compare_raw_vs_clever(
        X, W, Y, spec, gap, seed=5, gt="valence", with_po=False,
        n_estimators=40, with_block_train=False,
    )
    assert row["domain_auc_pool"] >= 0.78


def test_block_aware_forest_follows_pi_not_column_scale():
    """π must change feature *sampling*; RF splits ignore monotone column scales."""
    X, W, Y, spec = make_synthetic_shift(
        n=320, d_text=28, d_vad=4, gt="valence", mean_shift=1.25, seed=9,
    )
    pi_gt = np.array([0.04, 0.88, 0.04, 0.04])
    pi_unif = np.full(4, 0.25)
    model_gt = fit_block_aware_forest(X, W, spec, pi_gt, n_estimators=50, seed=9, floor=0.0)
    mass = block_aware_importance(model_gt, spec)
    assert mass["valence"] >= mass["text"]
    assert mass["valence"] >= 0.35
    auc_pi, _, _ = cv_block_aware_auc(
        X, W, spec, pi_gt, seed=9, n_estimators=40, n_splits=4, floor=0.0,
    )
    auc_unif, _, _ = cv_block_aware_auc(
        X, W, spec, pi_unif, seed=9, n_estimators=40, n_splits=4, floor=0.0,
    )
    assert auc_pi >= auc_unif - 0.02
    row = compare_raw_vs_clever(
        X, W, Y, spec,
        decompose_modality_gap(X, W, spec, seed=9, n_splits=3, n_estimators=35, light=True),
        seed=9, gt="valence", with_po=False, n_estimators=35, n_splits=4,
    )
    assert "domain_auc_subspace" in row
    assert row["domain_auc_bawf"] >= row["domain_auc_subspace"] - 0.03
    scales = adaptive_group_scale(X, spec, pi_gt, floor=0.0)
    val = spec.slices[spec.names.index("valence")]
    text = spec.slices[spec.names.index("text")]
    assert float(scales[val].mean()) > float(scales[text].mean())
    idx = pi_column_replicates(spec, pi_gt, X.shape[1], floor=0.0)
    val_idx = set(range(val.start, val.stop))
    n_val = int(sum(int(j) in val_idx for j in idx))
    n_text = int(sum(spec.slices[0].start <= int(j) < spec.slices[0].stop for j in idx))
    assert n_val > spec.slices[1].stop - spec.slices[1].start
    assert n_val / max(n_text, 1) > (spec.slices[1].stop - spec.slices[1].start) / max(
        spec.slices[0].stop - spec.slices[0].start, 1
    )
