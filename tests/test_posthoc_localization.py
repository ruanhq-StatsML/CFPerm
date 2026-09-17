#!/usr/bin/env python3
"""Post-hoc localization calls MMD and PO-risk on pulled subset indices."""
from __future__ import annotations

import unittest

import numpy as np

from scripts.posthoc_localization import (
    localize,
    mmd_pair,
    pairwise_subset_mmd,
    pairwise_subset_po_risk,
    po_risk,
    subset_indices_from_feature,
    subset_indices_from_labels,
)


class PosthocLocalizationTests(unittest.TestCase):
    def test_pulls_indices_then_mmd_separates_shifted_groups(self):
        rng = np.random.default_rng(0)
        X0 = rng.normal(size=(80, 4))
        X1 = rng.normal(size=(80, 4)) + 1.6
        X = np.vstack([X0, X1])
        T = np.array([0] * 80 + [1] * 80)
        subsets = subset_indices_from_labels(T, prefix="T")
        self.assertEqual(sorted(subsets), ["T0", "T1"])
        self.assertEqual(len(subsets["T0"]), 80)
        rows = pairwise_subset_mmd(X, subsets, n_perm=20, seed=0)
        self.assertEqual(len(rows), 1)
        self.assertGreater(rows[0]["mmd"], 0.05)
        self.assertLessEqual(rows[0]["mmd_p"], 0.1)

    def test_null_labels_do_not_look_like_a_subset_hit(self):
        rng = np.random.default_rng(1)
        X = rng.normal(size=(160, 4))
        T = rng.integers(0, 2, size=160)
        Y = rng.normal(size=160)
        rows = pairwise_subset_mmd(X, subset_indices_from_labels(T), n_perm=20, seed=1)
        self.assertGreater(rows[0]["mmd_p"], 0.05)

    def test_po_risk_on_pair_beats_permuted_w_when_cate_is_there(self):
        rng = np.random.default_rng(2)
        X = rng.normal(size=(200, 5))
        W = np.array([0] * 100 + [1] * 100)
        Y = 0.2 * X[:, 0] + 1.8 * W * X[:, 1] + 0.15 * rng.normal(size=200)
        obs = po_risk(X, Y, W)
        null = po_risk(X, Y, rng.permutation(W))
        self.assertGreater(obs, null)
        subsets = {"T0": np.arange(100), "T1": np.arange(100, 200)}
        rows = pairwise_subset_po_risk(X, Y, subsets, n_perm=20, seed=2)
        self.assertLessEqual(rows[0]["po_p"], 0.1)
        self.assertIn("mean_Y_a", rows[0])

    def test_quartile_indices_and_localize_keeps_means_plus_pairs(self):
        rng = np.random.default_rng(3)
        x = np.concatenate([rng.normal(size=60), rng.normal(size=60) + 3.0])
        X = np.column_stack([x, rng.normal(size=120)])
        Y = (x > np.median(x)).astype(float)
        T = (np.arange(120) >= 60).astype(int)
        bins = subset_indices_from_feature(x, n_bins=4)
        self.assertGreaterEqual(len(bins), 2)
        rec = localize(X, Y, group_labels=T, feature=x, n_perm=15, seed=3)
        self.assertIn("groups", rec)
        self.assertIn("bins", rec)
        self.assertIn("conditional_mean", rec["groups"])
        self.assertIsNotNone(rec["bins"]["conditional_mean"][next(iter(rec["bins"]["conditional_mean"]))]["mean_X"])
        self.assertTrue(rec["groups"]["pairwise"])
        self.assertIn("mmd", rec["groups"]["pairwise"][0])
        self.assertIn("po_risk", rec["groups"]["pairwise"][0])

    def test_mean_table_is_emitted_from_already_computed_values(self):
        from scripts.posthoc_localization import markdown_conditional_means

        md = "\n".join(
            markdown_conditional_means(
                [
                    {
                        "title": "toy",
                        "loc_means": {"T0": {"n": 10, "mean_Y": 0.2}, "T1": {"n": 10, "mean_Y": 0.8}},
                    }
                ]
            )
        )
        self.assertIn("Conditional mean", md)
        self.assertIn("0.200", md)
        self.assertIn("0.800", md)

    def test_mmd_pair_is_the_vendor_call(self):
        rng = np.random.default_rng(4)
        a = rng.normal(size=(40, 3))
        b = rng.normal(size=(40, 3)) + 2.0
        self.assertGreater(mmd_pair(a, b), mmd_pair(a, a + 0.01 * rng.normal(size=a.shape)))


if __name__ == "__main__":
    unittest.main()
