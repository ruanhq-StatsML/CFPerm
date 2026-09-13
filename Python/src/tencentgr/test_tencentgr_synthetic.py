"""Synthetic checks for the three-tower forward pass and the DR-PO learner."""

from __future__ import annotations

import numpy as np
import torch

from tencentgr.dataset import _as_emb
from tencentgr.dr_po_learner import dr_pseudo_outcome, fit_dr_pseudo_outcome
from tencentgr.three_tower import ThreeTowerModel


def synthetic_batch(n: int = 8, t: int = 6, d_user: int = 20, d_emb: int = 64):
    return {
        "user_features": torch.randn(n, d_user),
        "history_items": torch.randint(1, 50, (n, t)),
        "history_actions": torch.randint(0, 3, (n, t)),
        "history_embs": torch.randn(n, t, d_emb),
        "history_mask": torch.ones(n, t),
        "target_item": torch.randint(1, 50, (n,)),
        "target_action": torch.randint(0, 3, (n,)),
        "target_emb": torch.randn(n, d_emb),
    }


def test_parse_string_embedding():
    vec = _as_emb("[1.0, 2.0, 3.0]", dim=3)
    assert vec.shape == (3,)
    assert abs(float(vec[1]) - 2.0) < 1e-6
    pad = _as_emb(np.array([1.0, 2.0], dtype=np.float32), dim=4)
    assert pad.shape == (4,)
    assert float(pad[3]) == 0.0


def test_three_tower_forward():
    batch = synthetic_batch()
    model = ThreeTowerModel(user_dim=20, emb_dim=64, hidden_dim=32, tower_dim=16)
    out = model(batch)
    assert out["loss"].ndim == 0
    assert out["logits"].shape == (8, 8)
    out["loss"].backward()


def test_dr_po_recovers_mean_shift():
    rng = np.random.default_rng(0)
    n = 600
    p = 6
    W = np.r_[np.zeros(n // 2), np.ones(n - n // 2)].astype(int)
    X = rng.normal(size=(n, p))
    X[:, 0] += 0.8 * W
    Y = 0.3 * X[:, 0] + 0.2 * W + rng.normal(scale=0.2, size=n)
    fit = fit_dr_pseudo_outcome(X, Y, W, n_splits=3, seed=0, ridge_alpha=0.5)
    assert fit["po_risk_phi"] > 0
    assert abs(fit["phi_mean"] - 0.2) < 0.15
    phi = dr_pseudo_outcome(
        Y.astype(float),
        W.astype(float),
        np.zeros(n),
        np.zeros(n),
        np.full(n, 0.5),
    )
    assert phi.shape == (n,)


if __name__ == "__main__":
    test_parse_string_embedding()
    test_three_tower_forward()
    test_dr_po_recovers_mean_shift()
    print("synthetic tests ok")
