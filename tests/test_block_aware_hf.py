"""Offline tests for labeled-HF block encoding (no Hub download)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from block_aware_hf import (  # noqa: E402
    encode_adult_blocks,
    hashed_text,
    overlap_features,
    polarity_features,
    style_features,
    write_dashboard_latex,
    write_dataset_latex,
)
from clever_covariate_gap import inject_block_shift  # noqa: E402


def test_encode_adult_blocks_no_sex_leak():
    df = pd.DataFrame({
        "age": [20, 40, 55, 30],
        "workclass": ["Private", "Private", "Federal-gov", "Private"],
        "fnlwgt": [1e5, 2e5, 1.2e5, 9e4],
        "education": ["Bachelors", "HS-grad", "Masters", "HS-grad"],
        "education.num": [13, 9, 14, 9],
        "marital.status": ["Never-married", "Married-civ-spouse", "Divorced", "Never-married"],
        "occupation": ["Sales", "Craft-repair", "Exec-managerial", "Sales"],
        "relationship": ["Not-in-family", "Husband", "Unmarried", "Own-child"],
        "race": ["White", "White", "Black", "Asian-Pac-Islander"],
        "sex": ["Female", "Male", "Female", "Male"],
        "capital.gain": [0, 5000, 0, 0],
        "capital.loss": [0, 0, 200, 0],
        "hours.per.week": [20, 40, 50, 35],
        "native.country": ["United-States", "United-States", "India", "United-States"],
        "income": ["<=50K", ">50K", ">50K", "<=50K"],
    })
    X, W, Y, spec = encode_adult_blocks(df)
    assert spec.names == ["demography", "work", "hours", "capital"]
    assert X.shape[0] == 4
    assert set(np.unique(W)) == {0, 1}
    assert Y.sum() == 2
    assert spec.slices[-1].stop == X.shape[1]
    # relationship/sex must not be copied in as raw columns
    assert X.shape[1] < 40


def test_text_style_polarity_shapes():
    texts = ["This movie is excellent and wonderful.", "A terrible boring waste of time!!!"]
    st = style_features(texts)
    pol = polarity_features(texts)
    ht = hashed_text(texts, d_text=12, seed=0)
    assert st.shape == (2, 6)
    assert pol.shape == (2, 4)
    assert ht.shape == (2, 12)
    assert pol[0, 0] > pol[1, 0]
    assert pol[1, 1] > pol[0, 1]


def test_overlap_and_inject_on_toy_blocks():
    prem = ["the cat sat on the mat", "a dog ran"]
    hyp = ["the cat sat", "a cat sat"]
    ov = overlap_features(prem, hyp)
    assert ov.shape == (2, 4)
    X = np.hstack([np.random.default_rng(0).normal(size=(40, 6)), np.zeros((40, 4))])
    W = np.array([0] * 20 + [1] * 20)
    from clever_covariate_gap import ModalitySpec
    spec = ModalitySpec(names=["a", "b"], slices=[slice(0, 6), slice(6, 10)])
    X2 = inject_block_shift(X, W, spec, "b", alpha=1.2)
    assert np.linalg.norm(X2[W == 1][:, 6:] - X[W == 1][:, 6:]) > 0


def test_latex_writers(tmp_path):
    summary = {
        "slug": "toy",
        "title": "Toy board",
        "dataset": "example/ds",
        "label": "y",
        "batch": "W = domain",
        "n0": 10,
        "n1": 10,
        "gap": {
            "pi_vimp": {"a": 0.2, "b": 0.8},
            "pi_auc": {"a": 0.3, "b": 0.7},
            "pi_mmd": {"a": 0.1, "b": 0.9},
            "pi_lomo": {"a": 0.4, "b": 0.6},
            "pi_consensus": {"a": 0.25, "b": 0.75},
        },
        "observational_auc": {
            "domain_auc_raw": 0.7,
            "domain_auc_pi_rf": 0.71,
            "domain_auc_pool": 0.69,
            "domain_auc_stack_xz": 0.72,
            "domain_auc_clever": 0.68,
        },
        "comparisons": [{
            "name": "inject b",
            "gt": "b",
            "domain_auc_raw": 0.80,
            "domain_auc_pi_rf": 0.88,
            "domain_auc_pool": 0.82,
            "domain_auc_stack_xz": 0.85,
            "domain_auc_clever": 0.81,
            "pi_rf_on_gt": 0.6,
        }],
    }
    p = tmp_path / "t.tex"
    write_dataset_latex(summary, p)
    text = p.read_text(encoding="utf-8")
    assert "consensus" in text
    assert "inject b" in text
    dash = tmp_path / "dash.tex"
    write_dashboard_latex([summary], None, dash)
    dtext = dash.read_text(encoding="utf-8")
    assert "\\includegraphics" in dtext
    assert "toy_gap_shares.png" in dtext
    assert "\\begin{document}" in dtext
