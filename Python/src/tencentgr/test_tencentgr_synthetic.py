"""Synthetic checks for the three-tower forward pass and the DR-PO learner."""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch

from tencentgr.behavior_features import explode_seq, standard_scale_behavior, synthesize_user_behavior
from tencentgr.dataset import USER_FEAT_DIM, _as_emb, user_row_to_vec
from tencentgr.dr_po_learner import dr_pseudo_outcome, fit_dr_pseudo_outcome
from tencentgr.three_tower import ThreeTowerModel


def synthetic_batch(n: int = 8, t: int = 6, d_user: int = 20, d_emb: int = 64):
    return {
        "user_features": torch.randn(n, d_user),
        "history_items": torch.randint(1, 50, (n, t)),
        "history_actions": torch.randint(0, 3, (n, t)),
        "history_embs": torch.randn(n, t, d_emb),
        "history_mask": torch.ones(n, t),
        "history_emb_obs": torch.ones(n, t),
        "target_item": torch.randint(1, 50, (n,)),
        "target_action": torch.randint(0, 3, (n,)),
        "target_emb": torch.randn(n, d_emb),
        "target_emb_obs": torch.ones(n),
        "any_click": torch.randint(0, 2, (n,)).float(),
        "any_conversion": torch.randint(0, 2, (n,)).float(),
    }


def test_parse_string_embedding():
    vec = _as_emb("[1.0, 2.0, 3.0]", dim=3)
    assert vec.shape == (3,)
    assert abs(float(vec[1]) - 2.0) < 1e-6
    pad = _as_emb(np.array([1.0, 2.0], dtype=np.float32), dim=4)
    assert pad.shape == (4,)
    assert float(pad[3]) == 0.0


def test_user_missing_policy():
    observed = user_row_to_vec({"103": 1, "104": 2, "105": 3, "109": 4, "106": [1], "107": [2], "108": [3], "110": [4]})
    missing = user_row_to_vec({})
    assert observed.shape == (USER_FEAT_DIM,)
    assert missing.shape == (USER_FEAT_DIM,)
    assert float(missing[1]) == 1.0
    assert float(observed[1]) == 0.0


def test_three_tower_forward():
    batch = synthetic_batch()
    model = ThreeTowerModel(user_dim=20, emb_dim=64, hidden_dim=32, tower_dim=16)
    out = model(batch)
    assert out["loss"].ndim == 0
    assert out["logits"].shape == (8, 8)
    assert "click_loss" in out and "conv_loss" in out
    assert "act_loss" not in out
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


def test_behavior_explode_and_scale():
    seq_df = pd.DataFrame(
        [
            {
                "user_id": 1,
                "seq": [
                    {"item_id": 10, "action_type": 0, "timestamp": 100},
                    {"item_id": 11, "action_type": 1, "timestamp": 200},
                    {"item_id": 11, "action_type": 2, "timestamp": 300},
                ],
            },
            {
                "user_id": 2,
                "seq": [
                    {"item_id": 20, "action_type": 0, "timestamp": 50},
                    {"item_id": 21, "action_type": 0, "timestamp": 60},
                ],
            },
        ]
    )
    events = explode_seq(seq_df)
    assert len(events) == 5
    assert set(events["user_id"]) == {1, 2}
    users = synthesize_user_behavior(events)
    assert len(users) == 2
    u1 = users.set_index("user_id").loc[1]
    assert int(u1["n_click"]) == 1
    assert int(u1["n_conversion"]) == 1
    assert abs(float(u1["engage_rate"]) - 2.0 / 3.0) < 1e-9
    scaled, scaler, cols = standard_scale_behavior(users)
    assert "n_click" in cols
    assert "z_n_click" in scaled.columns
    assert abs(float(np.mean(scaled["z_n_click"]))) < 1e-9
    assert float(scaler.scale_[cols.index("n_click")]) > 0


if __name__ == "__main__":
    test_parse_string_embedding()
    test_user_missing_policy()
    test_three_tower_forward()
    test_dr_po_recovers_mean_shift()
    test_behavior_explode_and_scale()
    print("synthetic tests ok")
