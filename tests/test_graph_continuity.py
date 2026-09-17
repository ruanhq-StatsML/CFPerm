#!/usr/bin/env python3
"""GraphSAGE ⊕ hand concat and user-list aggregation."""
from __future__ import annotations

import unittest

import numpy as np

from scripts.prototype_graph_continuity import (
    ALL_COLS,
    GRAPH_COLS,
    HAND_DIM,
    NODE_DIM,
    SAGE_DIM,
    SAGE_LAYERS,
    X_DIM,
    aggregate_user_lists,
    graphsage_mean,
    node_features,
    rank_pvalues,
    user_of,
)


class GraphContinuityTests(unittest.TestCase):
    def test_concat_dim(self):
        self.assertEqual(SAGE_DIM, NODE_DIM * (2**SAGE_LAYERS))
        self.assertEqual(HAND_DIM, 7)
        self.assertEqual(X_DIM, SAGE_DIM + HAND_DIM)
        self.assertEqual(len(ALL_COLS), X_DIM)
        self.assertTrue(ALL_COLS[-len(GRAPH_COLS) :] == list(GRAPH_COLS))

    def test_graphsage_mean_uses_neighbors(self):
        n = 4
        seeds = np.array([1.0, 0.0, 0.0, 0.0])
        isolated = graphsage_mean(node_features(n, [], seeds), [])
        linked = graphsage_mean(node_features(n, [(0, 1), (1, 2)], seeds), [(0, 1), (1, 2)])
        self.assertEqual(isolated.shape, (SAGE_DIM,))
        self.assertEqual(linked.shape, (SAGE_DIM,))
        self.assertFalse(np.allclose(isolated, linked))

    def test_user_list_mean_pool(self):
        rows = []
        for i in range(6):
            rec = {c: float(i) for c in ALL_COLS}
            rec.update({"y": i % 2, "batch": 0 if i < 4 else 1, "user_id": user_of(i)})
            rows.append(rec)
        agg = aggregate_user_lists(rows)
        self.assertTrue(all("user_id" in r and "n_list" in r for r in agg))
        self.assertEqual(sum(r["n_list"] for r in agg), 6)
        keys = {(r["batch"], r["user_id"]) for r in agg}
        self.assertEqual(len(keys), len(agg))

    def test_rank_p_small_when_T_jumps(self):
        T = np.array([0.0, 0.01, 0.0, 0.02, 1.0])
        p = rank_pvalues(T)
        self.assertAlmostEqual(p[0], 1.0)
        self.assertLess(p[-1], 0.3)

    def test_user_hash_is_stable(self):
        self.assertEqual(user_of(0), user_of(20))
        self.assertNotEqual(user_of(0), user_of(1))


if __name__ == "__main__":
    unittest.main()
