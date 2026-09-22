#!/usr/bin/env python3
"""Schema lock for hybrid-retrieval / Graph-RAG Hotpot tables."""
from __future__ import annotations

import csv
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "manuscript" / "hybrid_retrieval"
HYBRID = [
    "x_q_toks",
    "x_q_chars",
    "x_qmark",
    "x_q_ents",
    "x_n_cand",
    "x_js_overlap",
    "x_rank_corr",
    "x_sparse_margin",
    "x_dense_margin",
    "x_rrf_top1_mass",
    "x_fuse_uniq",
    "x_mean_q_overlap",
    "x_ks",
    "x_kd",
]
GRAPH = [
    "x_n_nodes",
    "x_n_edges",
    "x_mean_deg",
    "x_n_cc",
    "x_n_q_seeds",
    "x_seed_frac",
    "x_lcc_frac",
]
LEAK = ("gold", "supporting", "hit", "question", "answer", "chosen", "rejected")


class HybridRetrievalXyTests(unittest.TestCase):
    def _rows(self, name):
        path = OUT / name
        self.assertTrue(path.exists(), msg=path)
        with path.open() as f:
            return list(csv.DictReader(f))

    def test_hybrid_schema(self):
        for name in ("xy_hotpot_hybrid.csv", "xy_hotpot_hybrid_hop.csv"):
            rows = self._rows(name)
            self.assertEqual(len(rows), 1200, name)
            keys = list(rows[0].keys())
            self.assertEqual(keys[:2], ["y", "batch"])
            for c in HYBRID:
                self.assertIn(c, keys, name)
            low = {k.lower() for k in keys}
            for bad in LEAK:
                self.assertNotIn(bad, low, name)
            ys = {int(r["y"]) for r in rows}
            self.assertTrue(ys <= {0, 1}, name)

    def test_graph_adds_geometry_not_gold(self):
        rows = self._rows("xy_hotpot_graph.csv")
        self.assertEqual(len(rows), 1200)
        for c in GRAPH:
            self.assertIn(c, rows[0])
        self.assertNotIn("question", {k.lower() for k in rows[0]})

    def test_two_stream_is_sparse_vs_dense(self):
        rows = self._rows("xy_hotpot_two_stream.csv")
        self.assertEqual(len(rows), 2400)
        self.assertEqual({int(r["T"]) for r in rows}, {0, 1})
        self.assertEqual({r["stream"] for r in rows}, {"sparse", "dense"})

    def test_pair_shape_is_query_by_ten_paras(self):
        rows = self._rows("xy_hotpot_pairs.csv")
        self.assertEqual(len(rows), 12000)
        batches = [int(r["batch"]) for r in rows]
        self.assertEqual(min(batches), 0)
        self.assertEqual(max(batches), 1199)
        n_per = sum(1 for b in batches if b == 0)
        self.assertEqual(n_per, 10)
        for c in (
            "x_bm25",
            "x_dense",
            "x_rrf",
            "x_rank_sp",
            "x_rank_de",
            "x_q_overlap",
            "x_title_seed",
            "x_title_deg",
            "x_slot",
        ):
            self.assertIn(c, rows[0])
        ys = [int(r["y"]) for r in rows[:10]]
        self.assertEqual(sum(ys), 2)
        low = {k.lower() for k in rows[0]}
        for bad in LEAK:
            self.assertNotIn(bad, low)

    def test_pair_hop_keeps_y_changes_dense(self):
        nat = self._rows("xy_hotpot_pairs.csv")
        hop = self._rows("xy_hotpot_pairs_hop.csv")
        self.assertEqual([r["y"] for r in nat], [r["y"] for r in hop])
        later = [
            (a["x_dense"], b["x_dense"])
            for a, b in zip(nat, hop)
            if int(a["batch"]) >= 320
        ]
        self.assertTrue(any(x != y for x, y in later))
        nat = self._rows("xy_hotpot_hybrid.csv")
        hop = self._rows("xy_hotpot_hybrid_hop.csv")
        pre_n = [int(a["y"]) for a, b in zip(nat, hop) if int(a["batch"]) < 4]
        pre_h = [int(b["y"]) for a, b in zip(nat, hop) if int(a["batch"]) < 4]
        self.assertEqual(pre_n, pre_h)
        y1n = [int(r["y"]) for r in nat if int(r["batch"]) >= 4]
        y1h = [int(r["y"]) for r in hop if int(r["batch"]) >= 4]
        self.assertNotEqual(y1n, y1h)


if __name__ == "__main__":
    unittest.main()
