#!/usr/bin/env python3
"""Palm–Nagler OnlineARBootstrap: mixing weight, CI, fire rule."""
from __future__ import annotations

import unittest

import numpy as np

from agod.online_ar_bootstrap import (
    BETA,
    OnlineARBootstrap,
    first_significant,
    run_delta_bootstrap,
)
from agod.rf_probe import error_floor, hop_fires, shift_ratio


class OnlineARBootstrapTests(unittest.TestCase):
    def test_beta_is_sqrt2_minus_1(self):
        self.assertAlmostEqual(BETA, 2.0 ** 0.5 - 1.0, places=12)

    def test_rho_formula(self):
        boot = OnlineARBootstrap(n_boot=8, seed=0)
        boot.update(0.0)
        self.assertAlmostEqual(boot.rho(), 1.0 - 1.0 ** (-BETA), places=12)
        boot.update(0.0)
        self.assertAlmostEqual(boot.rho(), 1.0 - 2.0 ** (-BETA), places=12)

    def test_ci_needs_two_updates(self):
        boot = OnlineARBootstrap(n_boot=32, seed=0)
        boot.update(0.1)
        lo, hi = boot.ci()
        self.assertTrue(np.isnan(lo) and np.isnan(hi))
        boot.update(0.1)
        lo, hi = boot.ci()
        self.assertTrue(np.isfinite(lo) and np.isfinite(hi))
        self.assertLessEqual(lo, hi)

    def test_zero_stream_covers_null(self):
        boot = OnlineARBootstrap(n_boot=200, seed=1)
        for _ in range(40):
            boot.update(0.0)
        lo, hi = boot.ci()
        self.assertLessEqual(lo, 0.0)
        self.assertGreaterEqual(hi, 0.0)
        self.assertAlmostEqual(boot.mean, 0.0, places=12)

    def test_positive_shift_fires(self):
        rows, boot = run_delta_bootstrap(np.full(60, 0.8), n_boot=200, seed=2)
        self.assertGreater(boot.mean, 0.5)
        lo, hi = boot.ci()
        self.assertGreater(lo, 0.0)
        hit = first_significant(rows)
        self.assertIsNotNone(hit)
        self.assertTrue(hit["fire"])
        self.assertGreaterEqual(hit["t"], 2)

    def test_fire_is_ci_lo_positive(self):
        rows, _ = run_delta_bootstrap([0.0, 0.0, 1.0, 1.0, 1.0], n_boot=80, seed=3)
        for r in rows:
            if r["t"] < 2:
                self.assertFalse(r["fire"])
            else:
                self.assertEqual(r["fire"], bool(np.isfinite(r["lo"]) and r["lo"] > 0.0))


class HopFiresTests(unittest.TestCase):
    def test_first_hop_quiet(self):
        self.assertFalse(hop_fires(0.5, None, gate=1.5, e_floor=0.02))

    def test_below_floor_quiet(self):
        self.assertFalse(hop_fires(0.4, 0.01, gate=1.5, e_floor=0.02))

    def test_ratio_gate(self):
        self.assertTrue(hop_fires(0.45, 0.20, gate=1.5, e_floor=0.02))
        self.assertFalse(hop_fires(0.24, 0.20, gate=1.5, e_floor=0.02))

    def test_classification_floor(self):
        self.assertEqual(error_floor("acc", 80), 0.02)
        self.assertGreater(shift_ratio(0.3, 0.1, e_floor=0.02), 1.0)


if __name__ == "__main__":
    unittest.main()
