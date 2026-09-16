#!/usr/bin/env python3
"""Wiring: generic (y, batch, x_*) loader and live-facet catalog."""
from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from scripts.run_serving_gates import STREAMS, compact, load_xy_table, x_columns
from scripts.llm_audit_online_bootstrap_prototype import freeze_and_score, run_stream

ROOT = Path(__file__).resolve().parents[1]


class ServingGatesRunnerTests(unittest.TestCase):
    def test_catalog_has_the_three_live_facets(self):
        facets = {s["facet"] for s in STREAMS}
        self.assertEqual(
            facets,
            {"审核 / judge", "Graph-RAG 子图", "混合检索"},
        )
        self.assertTrue(any(s["regime"] == "quiet" for s in STREAMS))
        self.assertTrue(any(s["regime"] == "hop" for s in STREAMS))
        names = {s["name"] for s in STREAMS}
        self.assertIn("judge_multistep_quiet", names)
        self.assertIn("judge_multistep_hop", names)

    def test_load_xy_uses_all_x_cols_and_rejects_chosen(self):
        path = ROOT / "results" / "manuscript" / "llm_audit" / "xy_beavertails.csv"
        X, y, batch, cols = load_xy_table(path)
        self.assertEqual(X.shape[0], 1200)
        self.assertEqual(len(cols), X.shape[1])
        self.assertTrue(all(c.startswith("x_") for c in cols))
        self.assertEqual(set(np.unique(y)), {0, 1})
        self.assertGreaterEqual(int(batch.max()), 1)

    def test_graph_query_table_is_binary_y(self):
        path = ROOT / "results" / "manuscript" / "graph_rag_batches" / "xy_graph_query.csv"
        if not path.exists():
            self.skipTest("graph pack tables not built")
        X, y, batch, cols = load_xy_table(path)
        self.assertEqual(set(np.unique(y)), {0, 1})
        self.assertIn("x_n_edges", cols)
        self.assertIn("x_lcc_frac", cols)
        self.assertEqual(int(batch.max()) + 1, 15)

    def test_run_stream_on_generic_table(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(160, 3))
        y = (X[:, 0] > 0).astype(int)
        batch = np.repeat(np.arange(4), 40)
        rec = run_stream(
            "toy",
            "toy quiet",
            X,
            y,
            batch,
            regime="quiet",
            n_ref_batches=2,
            cut_batch=2,
            gate=1.5,
            seed=0,
            labeled_cut=False,
        )
        self.assertEqual(rec["n"], 160)
        self.assertIn("mu_ref", rec)
        self.assertFalse(rec["hops"][0]["fired"])
        row = compact(rec, STREAMS[0])
        self.assertEqual(row["name"], "toy")
        self.assertIn("mean_delta", row)

    def test_x_columns_ignores_y_batch(self):
        self.assertEqual(x_columns(["y", "batch", "x_a", "n_queries"]), ["x_a"])


if __name__ == "__main__":
    unittest.main()
