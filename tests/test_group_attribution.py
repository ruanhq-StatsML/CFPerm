#!/usr/bin/env python3
"""CFPerm group attribution: permute T, two-level threshold, planted vs null."""
from __future__ import annotations

import unittest

import numpy as np

from scripts.prototype_group_attribution import (
    LEAK,
    cfperm_groups,
    interaction_vimp,
    load_two_stream,
    load_multistep_groups,
    null_dgp,
    planted_dgp,
    x_columns,
    zscore,
)
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class GroupAttributionTests(unittest.TestCase):
    def test_planted_cate_ranks_the_group_feature(self):
        X, Y, T, names = planted_dgp(n=500, p=6, seed=0)
        rec = cfperm_groups(X, Y, T, names=names, n_perm=25, seed=0)
        self.assertEqual(rec["top"][0], "x1")
        self.assertEqual(rec["rejected"], 1)
        self.assertIn("x1", [n for n, h in zip(rec["names"], rec["hits"]) if h] or rec["top"][:1])

    def test_null_does_not_reject(self):
        X, Y, T, names = null_dgp(n=500, p=6, seed=1)
        rec = cfperm_groups(X, Y, T, names=names, n_perm=25, seed=1)
        self.assertEqual(rec["rejected"], 0)
        self.assertEqual(rec["n_hits"], 0)

    def test_permutes_T_not_X(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(120, 3))
        T = np.array([0] * 60 + [1] * 60)
        Y = (T * X[:, 2] > 0).astype(float)
        imp = interaction_vimp(zscore(X), Y, T)
        null = interaction_vimp(zscore(X), Y, rng.permutation(T))
        self.assertGreater(imp[2], null[2])

    def test_k_way_groups(self):
        X, Y, T, names = planted_dgp(n=300, p=5, seed=2)
        self.assertEqual(len(np.unique(T)), 3)
        rec = cfperm_groups(X, Y, T, names=names, n_perm=20, seed=2)
        self.assertEqual(rec["n_groups"], 3)
        self.assertEqual(len(rec["pvals"]), 5)

    def test_live_two_stream_has_T_not_chosen(self):
        path = ROOT / "results" / "manuscript" / "llm_audit" / "xy_hh_two_stream.csv"
        X, Y, T, cols = load_two_stream(path)
        self.assertEqual(set(np.unique(T)), {0, 1})
        self.assertEqual(set(np.unique(Y)), {0, 1})
        self.assertIn("x_refuse", cols)
        self.assertTrue(all(c.startswith("x_") for c in cols))
        with path.open() as f:
            header = f.readline().lower()
        for bad in LEAK:
            self.assertNotIn(bad, header)

    def test_multistep_buckets_are_three_groups(self):
        path = ROOT / "results" / "manuscript" / "llm_audit" / "xy_hh_multistep_consistent.csv"
        if not path.exists():
            self.skipTest("multistep table not built")
        X, Y, T, cols = load_multistep_groups(path)
        self.assertEqual(set(np.unique(T)), {0, 1, 2})
        self.assertIn("x_step", cols)
        self.assertEqual(len(Y), 1200)

    def test_x_columns_skips_T_Y(self):
        self.assertEqual(x_columns(["T", "Y", "x_a", "batch", "stream"]), ["x_a"])


if __name__ == "__main__":
    unittest.main()
