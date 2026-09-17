"""Smoke tests for Grad-OnlineRFPerm helpers (single-stream口径)."""
from __future__ import annotations

import torch
import torch.nn as nn

from agod.grad_rfperm import (
    alarm_rate,
    false_alarm_rate,
    init_grad_rfperm,
    init_grad_rfperm_layers,
    layer_grad_norms,
    lead_time,
    relative_grad_shares,
    unfrozen_grad_l2,
    update_grad_rfperm,
)


class Tiny(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1 = nn.Linear(4, 3)
        self.fc2 = nn.Linear(3, 1)

    def forward(self, x):
        return self.fc2(torch.relu(self.fc1(x))).squeeze(-1)


def test_unfrozen_l2_matches_concat():
    m = Tiny()
    x = torch.randn(8, 4)
    y = torch.randn(8)
    loss = ((m(x) - y) ** 2).mean()
    loss.backward()
    g = unfrozen_grad_l2(m)
    flat = torch.cat([p.grad.detach().float().reshape(-1) for p in m.parameters()])
    assert abs(g - float(torch.linalg.vector_norm(flat).item())) < 1e-5
    # freeze fc1 → global norm drops frozen grads
    for p in m.fc1.parameters():
        p.requires_grad_(False)
        p.grad = None
    # recompute only fc2 grads
    m.zero_grad(set_to_none=True)
    loss = ((m(x) - y) ** 2).mean()
    loss.backward()
    g2 = unfrozen_grad_l2(m)
    norms = layer_grad_norms(m, unfrozen_only=True)
    assert "fc1" not in norms
    assert "fc2" in norms
    assert abs(g2 - norms["fc2"]) < 1e-5


def test_single_stream_update():
    st = init_grad_rfperm("unfrozen_l2")
    assert not isinstance(st, dict)
    for _ in range(5):
        update_grad_rfperm(st, 1.0, burn_in=True)
    out = update_grad_rfperm(st, 5.0, burn_in=False, alpha=0.2, fdr="fixed")
    assert "p" in out and "reject" in out
    assert lead_time(3, 5) == -2
    layers = init_grad_rfperm_layers(["fc1", "fc2"])
    assert set(layers) == {"fc1", "fc2"}
    shares = relative_grad_shares({"fc1": 1.0, "fc2": 3.0})
    assert abs(shares["fc2"] - 0.75) < 1e-9


def test_far_is_alarms_over_batches():
    hist = [0, 0, 1, 0, 1, 1, 0, 0]
    # after burn=2: alarms at t=2,4,5 → 3/6
    assert abs(alarm_rate(hist, after=2) - 3 / 6) < 1e-12
    assert abs(false_alarm_rate(hist, after=2) - 3 / 6) < 1e-12
    assert alarm_rate([0, 0, 0], after=0) == 0.0

