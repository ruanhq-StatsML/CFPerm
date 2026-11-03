"""Tests for first-indices-per-category + performance recording."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "Python" / "src"))

from first_category_performance import (  # noqa: E402
    first_indices_per_category,
    record_first_category_performance,
    split_first_indices,
    union_indices,
)


def test_first_indices_earliest():
    y = np.array([0, 0, 1, 0, 1, 1, 2, 2, 2, 1])
    idx = first_indices_per_category(y, n_first=2, categories=[0, 1, 2])
    assert list(idx[0]) == [0, 1]
    assert list(idx[1]) == [2, 4]
    assert list(idx[2]) == [6, 7]


def test_first_indices_frac():
    y = np.array([3] * 10 + [5] * 8)
    idx = first_indices_per_category(y, frac=0.25, categories=[3, 5])
    # round(10*0.25)=2; round(8*0.25)=2
    assert len(idx[3]) == 2
    assert len(idx[5]) == 2
    assert list(idx[3]) == [0, 1]
    assert list(idx[5]) == [10, 11]


def test_record_performance():
    y = np.array([0, 1, 0, 1, 2, 2])
    pred = np.array([0, 1, 1, 1, 2, 0])
    first = first_indices_per_category(y, n_first=1)
    rec = record_first_category_performance(y, pred, first)
    # first hits: 0@0, 1@1, 2@4 — all correct
    assert rec["per_category"][0]["acc"] == 1.0
    assert rec["per_category"][1]["acc"] == 1.0
    assert rec["per_category"][2]["acc"] == 1.0
    assert rec["micro_acc"] == 1.0
    assert set(union_indices(first).tolist()) == {0, 1, 4}


def test_split_stays_inside_first_window():
    y = np.array([7] * 20)
    first = first_indices_per_category(y, n_first=8, categories=[7])
    tr, ev = split_first_indices(first, seed=0, eval_frac=0.35)
    assert set(tr[7]).isdisjoint(ev[7]) or len(first[7]) == 1
    assert set(tr[7]).union(ev[7]) == set(first[7].tolist())


if __name__ == "__main__":
    test_first_indices_earliest()
    test_first_indices_frac()
    test_record_performance()
    test_split_stays_inside_first_window()
    print("ok")
