#!/usr/bin/env python3
"""Catalog columns match the csvs on disk."""
from __future__ import annotations

import unittest
from pathlib import Path

from scripts.serving_features import cols, disk_x_cols
from scripts.list_serving_features import check_disk, render_md

ROOT = Path(__file__).resolve().parents[1]


class ServingFeatureCatalogTests(unittest.TestCase):
    def test_catalog_matches_disk(self):
        errors = check_disk()
        self.assertEqual(errors, [])

    def test_judge_has_refuse_not_chosen(self):
        names = cols("judge")
        self.assertIn("x_refuse", names)
        self.assertEqual(len(names), 13)
        self.assertNotIn("chosen", names)

    def test_multistep_adds_x_step(self):
        names = cols("judge_multistep")
        self.assertEqual(names[-1], "x_step")
        self.assertEqual(len(names), 14)
        self.assertTrue(names[:-1] == cols("judge"))

    def test_graph_batch_is_mean_and_std(self):
        q = cols("graph_query")
        b = cols("graph_batch")
        self.assertEqual(len(q), 7)
        self.assertEqual(len(b), 14)
        self.assertEqual(b[:7], q)
        self.assertEqual(b[7:], [f"{c}_std" for c in q])

    def test_hybrid_has_no_gold_column(self):
        names = cols("hybrid")
        self.assertEqual(len(names), 14)
        blob = " ".join(names).lower()
        for bad in ("gold", "hit", "supporting", "question"):
            self.assertNotIn(bad, blob)

    def test_pair_is_ten_slot_geometry(self):
        names = cols("pair")
        self.assertEqual(
            names,
            [
                "x_bm25",
                "x_dense",
                "x_rrf",
                "x_rank_sp",
                "x_rank_de",
                "x_q_overlap",
                "x_title_seed",
                "x_title_deg",
                "x_slot",
            ],
        )

    def test_markdown_lists_live_gate_tables(self):
        md = render_md()
        self.assertIn("x_n_edges", md)
        self.assertIn("x_js_overlap", md)
        self.assertIn("run_serving_gates.py", md)

    def test_disk_helper_reads_header(self):
        path = ROOT / "results" / "manuscript" / "llm_audit" / "xy_beavertails.csv"
        self.assertEqual(disk_x_cols(path), cols("judge"))


if __name__ == "__main__":
    unittest.main()
