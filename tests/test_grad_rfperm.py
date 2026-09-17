"""Smoke tests for Grad-OnlineRFPerm helpers."""
from __future__ import annotations

import torch
import torch.nn as nn

from agod.grad_rfperm import (
    earliest_layer_reject,
    init_grad_rfperm,
    layer_grad_norms,
    lead_time,
    relative_grad_shares,
    update_grad_rfperm,
)


class Tiny(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc1 = nn.Linear(4, 3)
        self.fc2 = nn.Linear(3, 1)

    def forward(self, x):
        return self.fc2(torch.relu(self.fc1(x))).squeeze(-1)


def test_layer_grad_norms_and_shares():
    m = Tiny()
    x = torch.randn(8, 4)
    y = torch.randn(8)
    loss = ((m(x) - y) ** 2).mean()
    loss.backward()
    norms = layer_grad_norms(m)
    assert set(norms) >= {"fc1", "fc2"}
    assert all(v >= 0 for v in norms.values())
    shares = relative_grad_shares(norms)
    assert abs(sum(shares.values()) - 1.0) < 1e-5


def test_update_and_lead_time():
    states = init_grad_rfperm(["fc1"])
    st = states["fc1"]
    for _ in range(5):
        update_grad_rfperm(st, 1.0, burn_in=True)
    # spike
    out = update_grad_rfperm(st, 5.0, burn_in=False, alpha=0.2, fdr="fixed")
    assert "p" in out and "reject" in out
    assert lead_time(3, 5) == -2
    assert lead_time(None, 5) is None
    er = earliest_layer_reject(states, after=0)
    assert "__any__" in er
