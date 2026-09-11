"""Unit tests for ChronoBerg AGOD helpers (no Hugging Face download)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_chronoberg_agod as agod  # noqa: E402


def test_dummy_labels_are_arange():
    X0 = np.zeros((4, 3))
    X1 = np.ones((5, 3))
    n = len(X0) + len(X1)
    Y = np.arange(n)
    assert Y.tolist() == list(range(9))
    assert Y.shape[0] == n


def test_softmax_and_msg_route_to_drifting_modality():
    auc = {"text": 0.51, "valence": 0.97, "arousal": 0.52}
    vimp = {"text": 0.1, "valence": 0.8, "arousal": 0.1}
    po = {"text": 0.05, "valence": 0.9, "arousal": 0.05}
    g = agod.msg_from_components(auc, vimp, po, gamma=1.0)
    alpha = agod.softmax(g, tau=0.2)
    assert abs(alpha.sum() - 1.0) < 1e-8
    assert int(np.argmax(alpha)) == agod.MODALITIES.index("valence")
    assert alpha[1] > 0.6


def test_distill_step_reduces_global_loss():
    rng = np.random.default_rng(0)
    d = 3 * agod.D_MOD
    X = rng.normal(size=(40, d))
    student = agod.LinearStudent(d, agod.D_OUT, rng)
    alpha = np.array([0.1, 0.8, 0.1])
    before = agod.mse(student.embed(X, False), student.embed(X, True))
    for _ in range(8):
        student.step(X, alpha, lr=0.2)
    after = agod.mse(student.embed(X, False), student.embed(X, True))
    assert after < before


def test_recall_perfect_when_query_equals_gallery():
    rng = np.random.default_rng(1)
    Z = rng.normal(size=(20, 8))
    assert agod.recall_at_k(Z, Z, k=1) == 1.0
