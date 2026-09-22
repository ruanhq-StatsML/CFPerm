#!/usr/bin/env python3
"""Multimodal CFPerm: permute T, importance is a modality block."""
from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from scripts.prototype_multimodal_attribution import (
    GRAPH_BLOCKS,
    HYBRID_BLOCKS,
    LEAK,
    block_index,
    cfperm_blocks,
    load_graph_community,
    load_hybrid_dense_hop,
    null_multimodal,
    planted_multimodal,
)

ROOT = Path(__file__).resolve().parents[1]


class MultimodalAttributionTests(unittest.TestCase):
    def test_planted_hits_vision_block(self):
        X, Y, T, names, blocks = planted_multimodal(n=700, seed=0)
        rec = cfperm_blocks(X, Y, T, blocks, names=names, n_perm=25, seed=0)
        self.assertEqual(rec["block"]["top"][0], "vision")
        self.assertEqual(rec["block"]["rejected"], 1)
        self.assertIn("vision", rec["block"]["hit_names"] or rec["block"]["top"][:1])
        self.assertNotIn("text", rec["block"]["hit_names"])

    def test_null_does_not_reject_blocks(self):
        X, Y, T, names, blocks = null_multimodal(n=700, seed=1)
        rec = cfperm_blocks(X, Y, T, blocks, names=names, n_perm=25, seed=1)
        self.assertEqual(rec["block"]["rejected"], 0)
        self.assertEqual(rec["block"]["n_hits"], 0)

    def test_block_index_covers_hybrid_columns(self):
        names = [c for cols in HYBRID_BLOCKS.values() for c in cols]
        idx = block_index(names, HYBRID_BLOCKS)
        self.assertEqual(set(idx), {"query", "sparse", "dense", "fusion"})
        covered = sorted(int(i) for arr in idx.values() for i in arr)
        self.assertEqual(covered, list(range(len(names))))

    def test_graph_blocks_partition_seven(self):
        names = [c for cols in GRAPH_BLOCKS.values() for c in cols]
        self.assertEqual(len(names), 7)
        idx = block_index(names, GRAPH_BLOCKS)
        self.assertEqual(len(idx["query"]), 2)
        self.assertEqual(len(idx["graph"]), 5)

    def test_live_hybrid_T_is_pack_version_not_chosen(self):
        path = ROOT / "results" / "manuscript" / "hybrid_retrieval" / "xy_hotpot_hybrid_hop.csv"
        if not path.exists():
            self.skipTest("hybrid hop table missing")
        X, Y, T, cols, blocks = load_hybrid_dense_hop(path)
        self.assertEqual(set(np.unique(T)), {0, 1})
        self.assertEqual(set(np.unique(Y)), {0, 1})
        self.assertIn("x_dense_margin", cols)
        self.assertEqual(set(blocks), {"query", "sparse", "dense", "fusion"})
        with path.open() as f:
            header = f.readline().lower()
        for bad in LEAK:
            self.assertNotIn(bad, header)

    def test_live_graph_stack_two_packs(self):
        quiet = ROOT / "results" / "manuscript" / "graph_rag_batches" / "xy_graph_query.csv"
        hop = ROOT / "results" / "manuscript" / "graph_rag_batches" / "xy_graph_query_hop.csv"
        if not quiet.exists() or not hop.exists():
            self.skipTest("graph pack tables missing")
        X, Y, T, cols, blocks = load_graph_community(quiet, hop)
        self.assertEqual(len(Y), 2400)
        self.assertEqual(set(np.unique(T)), {0, 1})
        self.assertEqual(int(np.sum(T == 0)), 1200)
        self.assertIn("x_n_edges", cols)
        self.assertEqual(set(blocks), {"query", "graph"})


if __name__ == "__main__":
    unittest.main()
