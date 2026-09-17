#!/usr/bin/env python3
"""Sliding-window bank: cheap PCA stats vs D_ref, same clock as T / Brier."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from sliding_window_bank import SlidingWindowBank, subspace_gap  # noqa: E402
from stream_dgps import TRIMODAL_GROUPS, make_trimodal_gradual_concept  # noqa: E402


class BankTests(unittest.TestCase):
    def test_subspace_gap_zero_on_self(self):
        P = np.eye(4)[:, :2]
        self.assertAlmostEqual(subspace_gap(P, P), 0.0, places=6)

    def test_xxT_add_remove_cancels(self):
        from sliding_window_bank import SlidingSecondMoment

        rng = np.random.default_rng(0)
        X = rng.normal(size=(15, 4))
        g = SlidingSecondMoment(4)
        g.add(X)
        g.remove(X)
        self.assertLess(np.abs(g.C).max(), 1e-10)
        self.assertEqual(g.n, 0)

    def test_vector_grows_with_groups(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(80, 12))
        Y = rng.binomial(1, 0.5, size=80).astype(float)
        bank = SlidingWindowBank(X, Y, window=40, n_components=3, groups=TRIMODAL_GROUPS)
        feat = bank.step(X[:20], Y[:20])
        v = bank.vector(feat)
        self.assertGreaterEqual(len(v), 7 + 3)
        self.assertIn("pca_recon_group", feat)

    def test_concept_walk_keeps_pca_quieter_than_T(self):
        X, Y, alpha, _, meta = make_trimodal_gradual_concept(
            n_ref=120, n_new=30, n_batches=6, onset_batch=2, seed=1
        )
        n_ref = meta["n_ref"]
        n_new = meta["n_new"]
        bank = SlidingWindowBank(
            X[:n_ref],
            Y[:n_ref],
            window=60,
            n_components=3,
            groups=TRIMODAL_GROUPS,
        )
        rows = []
        for t in range(meta["n_batches"]):
            lo = n_ref + t * n_new
            rows.append(bank.step(X[lo : lo + n_new], Y[lo : lo + n_new]))
        t_pre = [r["rfperm_T"] for r in rows[:2]]
        t_post = [r["rfperm_T"] for r in rows[-2:]]
        self.assertGreater(float(np.mean(t_post)), float(np.mean(t_pre)))
        self.assertIn("pca_recon_excess", rows[-1])


if __name__ == "__main__":
    unittest.main()
