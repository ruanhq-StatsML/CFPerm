"""Smoke tests for subset pairwise helpers."""
from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.run_subset_pairwise_analysis import _rank_tail_masks, _tertile_masks


def test_tertile_masks_balanced():
    s = pd.Series(np.arange(300, dtype=float))
    m = _tertile_masks(s)
    assert m["low"].sum() > 50 and m["high"].sum() > 50


def test_rank_tail_masks():
    s = pd.Series(np.arange(1, 101, dtype=float))  # 1=best
    m = _rank_tail_masks(s, frac=0.25)
    assert m["high"].sum() == 25
    assert m["low"].sum() == 25
    assert set(np.where(m["high"])[0]).issubset(set(range(25)))
