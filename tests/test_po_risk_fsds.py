"""Tests for PO-risk FSDS helpers (period W; no ATE claim)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from po_risk_fsds import (  # noqa: E402
    blend_cmean_po_scores,
    bootstrap_pi_select,
    fit_period_po,
    po_feature_table,
    po_help_select,
    rare_pos_n_splits,
)


def test_rare_pos_n_splits_caps_by_positives():
    assert rare_pos_n_splits(np.array([0, 0, 1]), prefer=5) == 2
    assert rare_pos_n_splits(np.array([1, 1, 1, 0, 0]), prefer=5) == 3
    assert rare_pos_n_splits(np.zeros(10), prefer=5) == 2


def test_fit_period_po_returns_vimp_and_risk():
    rng = np.random.default_rng(0)
    n, d = 200, 6
    X = rng.normal(size=(n, d))
    W = np.concatenate([np.zeros(n // 2, int), np.ones(n - n // 2, int)])
    # plant a period shift on col 0
    X[W == 1, 0] += 1.5
    fit = fit_period_po(X, W, seed=0, n_trees=20)
    assert fit["risk"] >= 0.0
    assert fit["vimp"].shape == (d,)
    assert abs(fit["vimp"].sum() - 1.0) < 1e-6 or fit["vimp"].sum() > 0
    assert "not an ATE" in fit["note"]
    tab = po_feature_table([f"f{i}" for i in range(d)], fit["vimp"])
    assert tab.iloc[0]["feature"] == "f0" or tab["po_share"].iloc[0] > 0


def test_bootstrap_pi_select_rare_pos():
    rng = np.random.default_rng(1)
    n, d = 120, 8
    X = rng.normal(size=(n, d))
    y = np.zeros(n, int)
    y[:3] = 1
    X[y == 1, 2] += 3.0
    cols = [f"c{i}" for i in range(d)]
    selected, tab = bootstrap_pi_select(X, y, cols, k=3, n_boot=25, seed=1)
    assert len(selected) == 3
    assert "c2" in selected
    assert int(tab["n_boot_ok"].iloc[0]) > 0
    assert set(tab["feature"]) == set(cols)


def test_po_help_select_pool():
    rng = np.random.default_rng(2)
    cols = [f"g{i}" for i in range(10)]
    n1, n2 = 80, 90
    g1 = pd.DataFrame({c: rng.normal(size=n1) for c in cols})
    g2 = pd.DataFrame({c: rng.normal(size=n2) for c in cols})
    g2["g0"] = g2["g0"] + 2.0
    g1["y_convert"] = 0
    g2["y_convert"] = 0
    report = po_help_select(g1, g2, cols, k=5, seed=2, max_n=200)
    assert len(report["pool"]) == 8  # k+3
    assert report["risk"] >= 0.0
    assert "not an ATE" in report["note"]
    blend = blend_cmean_po_scores(
        cols,
        np.abs(g2[cols].mean().to_numpy() - g1[cols].mean().to_numpy()),
        report["feature_table"].set_index("feature").loc[cols, "po_vimp"].to_numpy(),
        alpha=0.5,
    )
    assert list(blend.columns)[:3] == ["feature", "cmean_abs_norm", "po_vimp_norm"]
