#!/usr/bin/python3
"""OnlineRFPerm + LLM routing: Y stays out of X; FAR helpers run."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from online_rfperm_with_llm import (  # noqa: E402
    compute_vimp,
    default_adapters,
    empirical_pval_large,
    make_random_noise,
    make_stationary,
    online_rfperm_with_llm,
    split_by_vimp,
)


class OnlineRfpermLlmTests(unittest.TestCase):
    def test_empirical_p_is_small_when_t_is_large(self):
        stream = np.array([0.0, 0.1, 0.0, 0.05, 5.0])
        p = empirical_pval_large(stream, burnin=2)
        self.assertLess(p[-1], 0.3)

    def test_y_never_enters_x_and_addis_keys(self):
        X, Y, _ = make_stationary(80, 20, 6, p=6, seed=1)
        llm, dl = default_adapters(seed=1)
        rec = online_rfperm_with_llm(
            X,
            Y,
            llm_adapter=llm,
            dl_adapter=dl,
            llm_cols=(3, 4, 5),
            dl_cols=(0, 1, 2),
            n_ref=80,
            n_new=20,
            seed=1,
        )
        self.assertEqual(X.shape[1], 6)
        self.assertIn("addis", rec["detectors"])
        self.assertEqual(len(rec["MSE_list"]), 6)
        self.assertTrue(np.all(np.isfinite(rec["MSE_list"])))

    def test_vimp_split_covers_p(self):
        X, Y, _ = make_stationary(60, 15, 4, p=8, seed=2)
        vimp = compute_vimp(X[:60], Y[:60], seed=2)
        high, low = split_by_vimp(vimp, n_high=3)
        self.assertEqual(len(high) + len(low), 8)
        self.assertEqual(len(set(high) | set(low)), 8)

    def test_random_noise_has_no_onset(self):
        _, _, meta = make_random_noise(40, 10, 3, p=4, seed=3)
        self.assertIsNone(meta["onset_batch"])
        self.assertEqual(meta["kind"], "random_noise")


if __name__ == "__main__":
    unittest.main()
