#!/usr/bin/env python3
"""Gradual concept: no hop, latency is excess area."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from gradual_latency import (  # noqa: E402
    delay_obs,
    excess_area,
    first_index,
    pre_onset_false_alarms,
    summarize_detector,
)
from stream_dgps import make_trimodal_gradual_concept  # noqa: E402


class GradualLatencyTests(unittest.TestCase):
    def test_alpha_walks_after_onset_no_jump(self):
        X, Y, alpha, tau, meta = make_trimodal_gradual_concept(
            n_ref=80, n_new=20, n_batches=6, onset_batch=2, seed=0
        )
        t0 = meta["onset_tau"]
        self.assertTrue(np.all(alpha[:t0] == 0.0))
        self.assertGreater(float(alpha[-1]), 0.9)
        d = np.diff(alpha[t0:])
        self.assertTrue(np.all(d >= -1e-12))
        self.assertLess(float(np.max(d)), 0.2)

    def test_first_index_and_infinite_delay(self):
        self.assertEqual(first_index([0, 0, 1, 1], start=0), 2)
        self.assertIsNone(first_index([0, 0, 0], start=0))
        self.assertIsNone(delay_obs(None, 2, 40))
        self.assertEqual(delay_obs(5, 2, 40), 120)

    def test_excess_area_until_end_if_never(self):
        loss = [0.2, 0.2, 0.3, 0.4]
        oracle = [0.2, 0.2, 0.2, 0.2]
        self.assertAlmostEqual(excess_area(loss, oracle, 2, None), 0.3)
        self.assertAlmostEqual(excess_area(loss, oracle, 2, 3), 0.1)

    def test_hop_false_alarms_pre_onset(self):
        flags = [False, True, False, True]
        self.assertEqual(pre_onset_false_alarms(flags, 2), 1)
        s = summarize_detector(
            "hop",
            flags,
            t0=2,
            n_new=10,
            alpha_batch=[0, 0, 0.1, 0.4],
            loss=[1, 1, 2, 3],
            oracle=[1, 1, 1, 1],
        )
        self.assertEqual(s["hat_batch"], 3)
        self.assertEqual(s["false_alarms_pre"], 1)
        self.assertAlmostEqual(s["area_until_hat"], 1.0)


if __name__ == "__main__":
    unittest.main()
