#!/usr/bin/env python3
"""Graph-RAG batch-aggregate table: schema, aggregation, community rewire."""
from __future__ import annotations

import csv
import unittest
from pathlib import Path

import numpy as np

from scripts.prototype_graph_pack_batch_agg import (
    BATCH_COLS,
    GRAPH_COLS,
    aggregate_batch_rows,
    co_mention_edges,
    graph_features,
    pack_usable,
    rewire_edges,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "manuscript" / "graph_rag_batches"
LEAK = ("gold", "supporting", "hit", "question", "answer", "chosen", "rejected", "title")


class GraphPackBatchAggTests(unittest.TestCase):
    def test_pack_usable_needs_one_hop_when_gold_is_a_neighbour(self):
        seeds = np.array([1.0, 0.0, 0.0])
        gold = {2}
        isolated = pack_usable(3, [], seeds, gold)
        linked = pack_usable(3, [(0, 2)], seeds, gold)
        community = pack_usable(3, [(0, 2)], seeds, gold, community=True)
        self.assertEqual(isolated, 0)
        self.assertEqual(linked, 1)
        self.assertEqual(community, 1)

    def test_community_pack_is_largest_cc_not_seeds(self):
        seeds = np.array([1.0, 0.0, 0.0, 0.0])
        gold = {2, 3}
        local = pack_usable(4, [(2, 3)], seeds, gold, community=False)
        comm = pack_usable(4, [(2, 3)], seeds, gold, community=True)
        self.assertEqual(local, 0)
        self.assertEqual(comm, 1)

    def test_rewire_keeps_edge_count_changes_pairs(self):
        rng = np.random.default_rng(0)
        native = [(0, 1), (1, 2)]
        hopped = rewire_edges(6, len(native), rng)
        self.assertEqual(len(hopped), len(native))
        self.assertNotEqual(set(hopped), set(native))

    def test_aggregate_one_row_per_window(self):
        rows = []
        for i in range(4):
            rows.append(
                {
                    "y": i % 2,
                    "batch": 0 if i < 2 else 1,
                    "x_n_nodes": 1.0,
                    "x_n_edges": float(i),
                    "x_mean_deg": 0.1,
                    "x_n_cc": 0.2,
                    "x_n_q_seeds": 0.3,
                    "x_seed_frac": 0.4,
                    "x_lcc_frac": 0.5,
                }
            )
        agg = aggregate_batch_rows(rows)
        self.assertEqual(len(agg), 2)
        self.assertEqual(agg[0]["n_queries"], 2)
        self.assertAlmostEqual(agg[0]["y"], 0.5)
        self.assertAlmostEqual(agg[0]["x_n_edges"], 0.5)
        self.assertIn("x_n_edges_std", agg[0])
        for c in BATCH_COLS:
            self.assertIn(c, agg[0])

    def test_graph_features_scale_and_cc(self):
        seeds = np.ones(4)
        g0 = graph_features(4, [], seeds)
        g1 = graph_features(4, [(0, 1), (1, 2), (2, 3)], seeds)
        self.assertEqual(g0["x_n_nodes"], 0.4)
        self.assertGreater(g1["x_mean_deg"], g0["x_mean_deg"])
        self.assertGreater(g0["x_n_cc"], g1["x_n_cc"])
        self.assertGreater(g1["x_lcc_frac"], g0["x_lcc_frac"])

    def test_co_mention_is_symmetric_and_token_based(self):
        edges = co_mention_edges(["Red Sox", "White Sox", "Unrelated"])
        self.assertEqual(edges, [(0, 1)])

    def _rows(self, name):
        path = OUT / name
        self.assertTrue(path.exists(), msg=path)
        with path.open() as f:
            return list(csv.DictReader(f))

    def test_written_query_and_batch_schema(self):
        if not (OUT / "xy_graph_query.csv").exists():
            self.skipTest("run scripts/prototype_graph_pack_batch_agg.py first")
        q = self._rows("xy_graph_query.csv")
        b = self._rows("xy_graph_batch.csv")
        self.assertEqual(len(q), 1200)
        self.assertEqual(len(b), 15)
        self.assertEqual(list(q[0].keys())[:2], ["y", "batch"])
        for c in GRAPH_COLS:
            self.assertIn(c, q[0])
        for c in BATCH_COLS:
            self.assertIn(c, b[0])
        low = {k.lower() for k in list(q[0].keys()) + list(b[0].keys())}
        for bad in LEAK:
            self.assertNotIn(bad, low)
        self.assertEqual({int(r["n_queries"]) for r in b}, {80})

    def test_hop_changes_post_cut_graph_not_pre_cut(self):
        if not (OUT / "xy_graph_query_hop.csv").exists():
            self.skipTest("run scripts/prototype_graph_pack_batch_agg.py first")
        nat = self._rows("xy_graph_query.csv")
        hop = self._rows("xy_graph_query_hop.csv")
        pre_n = [r["x_n_edges"] for r in nat if int(r["batch"]) < 4]
        pre_h = [r["x_n_edges"] for r in hop if int(r["batch"]) < 4]
        self.assertEqual(pre_n, pre_h)
        post = [
            (a["x_n_cc"], b["x_n_cc"])
            for a, b in zip(nat, hop)
            if int(a["batch"]) >= 4
        ]
        self.assertTrue(any(x != y for x, y in post))
        y_pre = [a["y"] for a, b in zip(nat, hop) if int(a["batch"]) < 4]
        y_pre_h = [b["y"] for a, b in zip(nat, hop) if int(a["batch"]) < 4]
        self.assertEqual(y_pre, y_pre_h)
        y_post_n = [int(a["y"]) for a in nat if int(a["batch"]) >= 4]
        y_post_h = [int(b["y"]) for b in hop if int(b["batch"]) >= 4]
        self.assertNotEqual(y_post_n, y_post_h)


if __name__ == "__main__":
    unittest.main()
