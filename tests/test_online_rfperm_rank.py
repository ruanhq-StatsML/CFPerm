#!/usr/bin/python3
"""One sample: onlinePermOOB_rank. Last column is relevance; groups are slates."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from online_rfperm_rank import make_rank_df, onlinePermOOB_rank  # noqa: E402


class OnlinePermOOBRankTests(unittest.TestCase):
    def test_onlinePermOOB_rank_y_is_last_column(self):
        df, groups = make_rank_df(n_queries=40, slate=8, kind="rank_flip", onset_q=20, seed=3)
        rec = onlinePermOOB_rank(df, groups, k=4, ref_batch_size=160, batch_size=40, seed=3)
        self.assertEqual(df.shape[1], 7)
        self.assertEqual(len(groups), len(df))
        self.assertFalse(rec["y_in_x"])
        self.assertTrue(np.all(np.isfinite(rec["NDCG_list"])))
        self.assertIn("ndcg_1", rec)
        self.assertEqual(len(rec["NDCG_list"]), 4)


if __name__ == "__main__":
    unittest.main()
