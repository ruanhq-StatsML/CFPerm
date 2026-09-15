"""Smoke tests for RFPerm hallucination prototype."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

_PATH = Path(__file__).resolve().parents[1] / "scripts/agod/rfperm_hallucination_proto.py"
_SPEC = importlib.util.spec_from_file_location("rfperm_hallucination_proto", _PATH)
_MOD = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(_MOD)


def test_concept_cut_rates_in_range():
    _, meta = _MOD.make_hallucination_stream(cut=4, seed=0)
    assert 0.05 < meta["rate_before"] < 0.95
    assert 0.05 < meta["rate_after"] < 0.95


def test_rfperm_fires_near_cut():
    report = _MOD.run_proto(gate=1.35, cut=4, seed=0)
    assert report["first_fire_t"] is not None
    assert report["first_fire_t"] >= 3
    auc = report["hop_at_cut"]["ranking"]["auroc_po_risk0"]
    assert auc > 0.55


def test_instance_probe_ranks():
    inst = _MOD.demo_instance_probe(cut=4, seed=1)
    assert inst["auroc"] > 0.55
    assert inst["precision_at_10"] >= 0.4


def test_auroc_perfect():
    y = np.array([0, 0, 1, 1])
    s = np.array([0.1, 0.2, 0.8, 0.9])
    assert abs(_MOD.auroc(y, s) - 1.0) < 1e-9
