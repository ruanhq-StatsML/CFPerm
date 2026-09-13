"""Unit tests for AGOD controller (no Amazon shards / GPU required)."""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import agod
from agod.policies import POLICIES, select_policy


def test_mmd_positive_under_mean_shift():
    rng = np.random.default_rng(1)
    X0 = rng.normal(size=(50, 6))
    X1 = rng.normal(size=(50, 6)) + 1.0
    assert agod.rbf_mmd2(X0, X1, seed=0) > agod.rbf_mmd2(X0, X0 + 1e-6, seed=0)


def test_b5_raises_lr_on_high_po_modality():
    mods = ["text", "image"]
    rng = np.random.default_rng(2)
    b0 = {m: rng.normal(size=(40, 5)) for m in mods}
    b1 = {m: rng.normal(size=(40, 5)) + 0.3 for m in mods}
    y0 = rng.integers(0, 2, 40).astype(float)
    y1 = y0.copy()
    msg = SimpleNamespace(
        auc={m: 0.65 for m in mods},
        vimp={m: 0.1 for m in mods},
        po={"text": 0.05, "image": 0.95},
    )
    raw, decomp, gain = select_policy(
        "B5",
        msg=msg,
        blocks0=b0,
        blocks1=b1,
        y0=y0,
        y1=y1,
        mods=mods,
        tau=0.3,
        kappa=1.25,
        seed=0,
    )
    lr = agod.alpha_to_lr(raw, mods, gain=gain)
    assert decomp["con"]["image"] > decomp["con"]["text"]
    assert lr["image"] > lr["text"]


def test_all_policies_return_finite_lr():
    mods = ["text", "image"]
    rng = np.random.default_rng(3)
    b0 = {m: rng.normal(size=(30, 4)) for m in mods}
    b1 = {m: rng.normal(size=(30, 4)) for m in mods}
    y0 = rng.integers(0, 2, 30).astype(float)
    y1 = rng.integers(0, 2, 30).astype(float)
    msg = SimpleNamespace(
        auc={m: 0.6 for m in mods},
        vimp={m: 0.2 for m in mods},
        po={m: 0.5 for m in mods},
    )
    for pol in POLICIES:
        raw, decomp, gain = select_policy(
            pol,
            msg=msg,
            blocks0=b0,
            blocks1=b1,
            y0=y0,
            y1=y1,
            mods=mods,
            tau=0.3,
            kappa=1.0,
            seed=1,
        )
        lr = agod.alpha_to_lr(raw, mods, gain=gain)
        assert all(np.isfinite(lr[m]) and lr[m] > 0 for m in mods)


if __name__ == "__main__":
    test_mmd_positive_under_mean_shift()
    test_b5_raises_lr_on_high_po_modality()
    test_all_policies_return_finite_lr()
    print("test_agod_controller: OK")
