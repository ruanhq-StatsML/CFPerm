#!/usr/bin/env python3
"""Direct localization on the on-disk audit csvs."""
from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from scripts.run_posthoc_on_dataset import DATASETS, load_dataset, top_feature

ROOT = Path(__file__).resolve().parents[1]


class PosthocOnDatasetTests(unittest.TestCase):
    def test_catalog_points_at_audit_tables(self):
        names = {d["name"] for d in DATASETS}
        self.assertIn("hh_two_stream", names)
        self.assertIn("hh_multistep", names)
        for d in DATASETS:
            self.assertTrue(str(d["path"]).endswith(".csv"))

    def test_two_stream_pulls_T_not_chosen(self):
        path = ROOT / "results" / "manuscript" / "llm_audit" / "xy_hh_two_stream.csv"
        if not path.exists():
            self.skipTest("hh two-stream missing")
        X, Y, labels, cols, how = load_dataset(path)
        self.assertEqual(set(np.unique(labels)), {0, 1})
        self.assertEqual(set(np.unique(Y)), {0, 1})
        self.assertIn("x_n_toks", cols)
        self.assertIn("T", how)
        with path.open() as f:
            header = f.readline().lower()
        self.assertNotIn("chosen", header)

    def test_multistep_uses_step_buckets(self):
        path = ROOT / "results" / "manuscript" / "llm_audit" / "xy_hh_multistep_consistent.csv"
        if not path.exists():
            self.skipTest("multistep missing")
        X, Y, labels, cols, how = load_dataset(path)
        self.assertEqual(set(np.unique(labels)), {0, 1, 2})
        self.assertIn("step", how)
        feat, name = top_feature(X, cols)
        self.assertEqual(name, "x_n_toks")
        self.assertEqual(len(feat), len(Y))


if __name__ == "__main__":
    unittest.main()
