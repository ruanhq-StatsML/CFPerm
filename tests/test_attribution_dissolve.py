"""Tests for hierarchical attribution dissolve (rollup)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "tencent_gr"))
from attribution_dissolve import dissolve_scores, drill_hierarchy, jaccard_topk  # noqa: E402


def test_dissolve_sum_and_rank():
    leaf = pd.DataFrame(
        {
            "user_id": [1, 1, 2, 2, 3],
            "merchant_id": [10, 10, 10, 20, 20],
            "score": [1.0, 2.0, 3.0, 4.0, 5.0],
        }
    )
    u = dissolve_scores(leaf, score_col="score", by="user_id", how="sum")
    # user2=7, user3=5, user1=3
    assert list(u["user_id"].astype(int)) == [2, 3, 1]
    assert float(u.loc[u["user_id"] == 1, "score"].iloc[0]) == 3.0
    m = dissolve_scores(leaf, score_col="score", by="merchant_id", how="sum")
    assert float(m.loc[m["merchant_id"] == 10, "score"].iloc[0]) == 6.0


def test_drill_hierarchy_topk():
    leaf = pd.DataFrame(
        {
            "user_id": np.repeat(np.arange(5), 2),
            "merchant_id": np.tile([100, 200], 5),
            "shift_l2": np.linspace(1, 10, 10),
        }
    )
    hier = drill_hierarchy(
        leaf, score_col="shift_l2", levels=("user_id", "merchant_id"), how="sum", top_k={"user_id": 2}
    )
    assert hier["user_id"]["selected"].sum() == 2
    # top-2: {1,2} vs {2,3} → intersection 1 / union 3
    assert abs(jaccard_topk([1, 2, 3], [2, 3, 4], 2) - 1 / 3) < 1e-9
