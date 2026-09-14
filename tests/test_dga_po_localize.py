"""PO-risk tail localizes D_spe; mean-diff discretizes the subset."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from dga_po_localize import (  # noqa: E402
    apply_rules,
    discretize_subset,
    instance_po_risk,
    localize_one_hop,
    make_planted_subgroup_stream,
    po_tail_mask,
    run_dga_po_localize,
    subgroup_mean_diff,
    top_feature_hit,
)


def test_po_tail_keeps_high_scores():
    is_spe = np.array([True, True, True, True, False, False])
    risk = np.array([0.1, 9.0, 0.2, 8.0, 7.0, 0.0])
    mask = po_tail_mask(is_spe, risk, q=0.5, min_n=2)
    assert mask.sum() == 2
    assert bool(mask[1]) and bool(mask[3])
    assert not mask[0] and not mask[2]
    # non-spe high risk must not enter the tail
    assert not mask[4]


def test_mean_diff_ranks_shifted_coordinate():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(80, 6))
    high = np.zeros(80, dtype=bool)
    high[:30] = True
    X[high, 2] += 2.5
    ranked = subgroup_mean_diff(X, high)
    assert ranked[0]["j"] == 2
    assert ranked[0]["cohens_d"] > 1.0


def test_discretize_recovers_midpoint_rule():
    ranked = [
        {
            "j": 0,
            "name": "X0",
            "mean_high": 1.0,
            "mean_low": -1.0,
            "cohens_d": 1.5,
            "n_high": 20,
            "n_low": 20,
        }
    ]
    rules = discretize_subset(ranked, k=1)
    assert rules[0]["op"] == ">="
    assert abs(rules[0]["threshold"] - 0.0) < 1e-9
    X = np.array([[0.5], [-0.5], [0.0]])
    mask = apply_rules(X, rules)
    assert mask.tolist() == [True, False, True]


def test_planted_pocket_has_higher_po_and_recovers_features():
    stream = make_planted_subgroup_stream(
        n_batches=5, n_per=140, p=12, seed=7, plant_from=3
    )
    spe = 3
    loc = localize_one_hop(stream.X, stream.y, stream.batch, spe, q=0.35, k_features=3)
    on = stream.batch == spe
    plant = stream.meta["planted_mask"][on]
    po = loc["risk"][on]
    assert loc["n_tail"] < loc["n_spe"]
    assert loc["n_tail"] >= 4
    if plant.any() and (~plant).any():
        assert float(po[plant].mean()) > float(po[~plant].mean())
    hit = top_feature_hit(loc["ranked_features"], stream.meta["shift_coords"], k=4)
    assert hit["recall"] >= 0.5
    rec = loc["tail_mask"][on]
    # tail should be enriched for the planted pocket
    prec = float(np.sum(plant & rec) / max(int(rec.sum()), 1))
    base = float(plant.mean())
    assert prec > base


def test_dga_po_localize_runs_and_reports_rules():
    stream = make_planted_subgroup_stream(n_batches=5, n_per=90, p=10, seed=1)
    rec = run_dga_po_localize(stream, po_q=0.35, k_features=3)
    assert rec["method"] == "dga_po_ridge"
    assert len(rec["online_path"]) == 4
    assert rec["online_mse"] == rec["online_mse"]
    last = rec["history"][-1]
    assert last["n_tail"] <= last["n_spe"]
    assert len(last["alignments"]) == last["spe_domain"] + 1
    assert abs(sum(last["alpha_inst"]) - 1.0) < 1e-6
    assert last["top_features"]
    last_rules = rec["last_rules"]
    assert isinstance(last_rules, list)


def test_dga_alignments_accept_spe_mask():
    from attribution_adapter import _ridge_fit, dga_alignments

    stream = make_planted_subgroup_stream(n_batches=4, n_per=40, p=8, seed=2)
    tr = stream.batch < 3
    clf = _ridge_fit(stream.X[tr], stream.y[tr])
    full = dga_alignments(clf, stream.X[tr], stream.y[tr], stream.batch[tr], [0, 1, 2], 2)
    mask = np.zeros(int(tr.sum()), dtype=bool)
    mask[-8:] = True
    loc = dga_alignments(
        clf,
        stream.X[tr],
        stream.y[tr],
        stream.batch[tr],
        [0, 1, 2],
        2,
        spe_mask=mask,
    )
    assert full.shape == loc.shape == (3,)
    assert np.all(np.isfinite(loc))
    assert np.all(np.isfinite(full))
