"""Smoke tests for representation-first cold-start."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "Python" / "src"))

from category_shift_boards import run_cgcd  # noqa: E402
from repr_sep_cold_start import (  # noqa: E402
    ColdStartState,
    adapt_batch_repr_sep,
    known_prototypes,
    margin_separate,
)
from weighted_cold_start import fit_initial  # noqa: E402


def test_margin_separate_increases_gap():
    rng = np.random.default_rng(0)
    # known blob at 0; novel-like near 0.2 (inside margin)
    X_ref = rng.normal(0, 0.3, size=(40, 4))
    y_ref = np.zeros(40, dtype=int)
    feat = np.array([0, 1])
    proto, _ = known_prototypes(X_ref, y_ref, feat)
    X_live = rng.normal(0.15, 0.2, size=(20, 4))
    mem = np.linspace(0.2, 1.0, 20)
    out = margin_separate(X_live, mem, proto, feat, margin=1.5, lr=0.5, mem_thr=0.5, n_steps=4)
    assert out["n_pushed"] > 0
    assert out["gap_after"] >= out["gap_before"] - 1e-9


def test_adapt_batch_repr_sep_phases():
    rng = np.random.default_rng(1)
    p, n_k, n_n = 12, 60, 20
    X0 = rng.normal(0, 1, size=(n_k, p))
    y0 = rng.integers(0, 3, size=n_k)
    # novel lives on first 4 dims
    X1 = rng.normal(0, 1, size=(n_n + 30, p))
    X1[:n_n, :4] += 3.0
    y1 = np.concatenate([np.full(n_n, 3), rng.integers(0, 3, size=30)])
    classes = np.arange(4)
    clf, sc = fit_initial(X0, y0, seed=0, classes=classes)
    X_ref, X_live = sc.transform(X0), sc.transform(X1)
    state = ColdStartState()
    out = adapt_batch_repr_sep(
        clf, X_ref, y0, X_live, y1, state,
        seed=0, force_adapt=True, expand_min=100,  # force buffer-only first
        min_jaccard=0.05, top_k=4,
    )
    assert out["phase"] == "separate"
    assert out["buffer_size"] > 0
    assert out["sep"]["n_pushed"] is not None
    # second call with low expand_min should expand
    out2 = adapt_batch_repr_sep(
        out["clf"], X_ref, y0, X_live, y1, out["state"],
        seed=1, force_adapt=True, expand_min=8,
        min_jaccard=0.05, top_k=4,
    )
    assert out2["phase"] == "expand"
    assert out2["expand"]["fitted"] is True


def test_cgcd_repr_sep_runs():
    res = run_cgcd(seed=2026, mode="repr_sep", quick=True, expand_min=8)
    assert "mean_novel_acc" in res
    assert res["mode"] == "repr_sep"
    assert len(res["stages"]) >= 1


if __name__ == "__main__":
    test_margin_separate_increases_gap()
    test_adapt_batch_repr_sep_phases()
    test_cgcd_repr_sep_runs()
    print("ok")
