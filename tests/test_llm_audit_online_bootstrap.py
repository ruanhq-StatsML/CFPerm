#!/usr/bin/env python3
"""Prototype wiring: schema, frozen-ref trail, hop overlay, no chosen-as-Y."""
from __future__ import annotations

import csv
import unittest
from pathlib import Path

import numpy as np

from agod.rf_probe import brier_score, fit_online_rf
from scripts.llm_audit_online_bootstrap_prototype import (
    X_COLS,
    apply_preference_hop,
    freeze_and_score,
    last_two_hops,
    linear_mse,
    load_xy,
)

ROOT = Path(__file__).resolve().parents[1]
XY = ROOT / "results" / "manuscript" / "llm_audit"


class LlmAuditBootstrapPrototypeTests(unittest.TestCase):
    def test_load_rejects_nothing_and_keeps_schema(self):
        X, y, batch = load_xy(XY / "xy_beavertails.csv")
        self.assertEqual(X.shape, (1200, 13))
        self.assertEqual(set(np.unique(y)), {0, 1})
        self.assertEqual(int(batch.max()) + 1, 15)

    def test_all_label_tables_exist(self):
        for name in (
            "xy_hh_helpful_consistent.csv",
            "xy_hh_helpful_hop.csv",
            "xy_hh_harmless_consistent.csv",
            "xy_hh_harmless_hop.csv",
            "xy_beavertails.csv",
            "xy_wildguard.csv",
            "xy_toxicchat.csv",
        ):
            with (XY / name).open() as f:
                keys = {k.lower() for k in next(csv.DictReader(f)).keys()}
            for c in X_COLS:
                self.assertIn(c, keys, msg=name)
            self.assertNotIn("chosen", keys)
            self.assertNotIn("rejected", keys)

    def test_freeze_and_score_splits_ref_and_trail(self):
        X, y, batch = load_xy(XY / "xy_hh_helpful_consistent.csv")
        out = freeze_and_score(X, y, batch, n_ref_batches=4, seed=0)
        self.assertEqual(out["ref_batches"], [0, 1, 2, 3])
        self.assertEqual(out["trail_batches"][0], 4)
        self.assertEqual(len(out["scores"]), 11)
        self.assertTrue(np.isfinite(out["mu_ref"]))

    def test_hop_overlay_flips_only_after_cut(self):
        y = np.zeros(16, dtype=int)
        batch = np.repeat(np.arange(4), 4)
        y2, rate = apply_preference_hop(y, batch, cut_batch=2, flip_rate=1.0, seed=0)
        self.assertTrue(np.all(y2[batch < 2] == 0))
        self.assertTrue(np.all(y2[batch >= 2] == 1))
        self.assertGreater(rate, 0)

    def test_last_two_first_hop_quiet(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(160, 4))
        y = (X[:, 0] > 0).astype(int)
        batch = np.repeat(np.arange(4), 40)
        hops = last_two_hops(X, y, batch, gate=1.5, seed=0)
        self.assertFalse(hops[0]["fired"])
        self.assertIsNone(hops[0]["e_prev"])

    def test_brier_on_constant_probe(self):
        X = np.zeros((10, 2))
        y = np.ones(10, dtype=int)
        probe = fit_online_rf(X, y, seed=0, task="acc")
        self.assertLess(brier_score(probe, X, y), 1e-9)

    def test_linear_probe_recovers_easy_map(self):
        rng = np.random.default_rng(1)
        X = rng.normal(size=(200, 3))
        y = (X[:, 0] > 0).astype(int)
        batch = np.repeat(np.arange(5), 40)
        out = freeze_and_score(X, y, batch, n_ref_batches=2, seed=0)
        self.assertLess(out["mu_ref"], 0.2)
        self.assertTrue(np.isfinite(linear_mse(out["probe"], X, y)))


if __name__ == "__main__":
    unittest.main()
