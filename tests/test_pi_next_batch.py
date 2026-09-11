"""Next-batch group learning rates from π. Not online learning."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from clever_covariate_gap import inject_block_shift, make_synthetic_shift  # noqa: E402
from pi_next_batch import (  # noqa: E402
    block_learning_rates,
    next_batch_lr_eval,
)


def test_mean_preserving_boost_puts_mass_on_pi():
    from clever_covariate_gap import ModalitySpec

    spec = ModalitySpec(
        names=["text", "valence", "arousal", "dominance"],
        slices=[slice(0, 20), slice(20, 24), slice(24, 28), slice(28, 32)],
    )
    pi = np.array([0.05, 0.80, 0.10, 0.05])
    eta0 = 0.4
    eta = block_learning_rates(spec, pi, eta0=eta0, mode="boost", floor=0.0)
    assert np.isclose(eta.mean(), eta0)
    assert eta[1] > 4 * eta[0]
    uni = block_learning_rates(spec, pi, eta0=eta0, mode="uniform")
    assert np.allclose(uni, eta0)
    damp = block_learning_rates(spec, pi, eta0=eta0, mode="damp", floor=0.0)
    assert damp[0] > damp[1]


def test_pi_boost_beats_uniform_in_few_steps():
    """Wide text + valence shift: reallocating steps to valence helps before convergence."""
    X, W, Y, spec = make_synthetic_shift(
        n=720, d_text=36, d_vad=4, gt="valence", mean_shift=1.15, seed=8,
    )
    row = next_batch_lr_eval(
        X, W, Y, spec, gt="valence", seed=8, n_steps=18, eta0=0.5, n_estimators=35,
    )
    auc = row["domain_auc"]
    assert auc["pi_boost"] >= auc["uniform"] - 0.01
    assert auc["pi_boost"] >= auc["pi_damp"] - 0.01
    assert auc["oracle"] >= auc["uniform"] - 0.02
    assert row["mass_on_gt"]["pi_boost"] >= row["mass_on_gt"]["uniform"] - 0.02


def test_inject_pi_boost_not_vimp():
    X, W, Y, spec = make_synthetic_shift(
        n=700, d_text=36, d_vad=4, gt="text", mean_shift=0.12, seed=2,
    )
    X = inject_block_shift(X, W, spec, "valence", alpha=1.2)
    row = next_batch_lr_eval(
        X, W, Y, spec, gt="valence", seed=2, n_steps=18, eta0=0.5, n_estimators=35,
    )
    assert row["pi"]["valence"] >= row["pi_vimp"]["valence"] - 0.05
    assert row["domain_auc"]["pi_boost"] >= row["domain_auc"]["vimp_boost"] - 0.02
