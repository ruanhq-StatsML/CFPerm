#!/usr/bin/python3
"""Recsys grain slices: no Y in X; FSDS recovers planted order columns."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from recsys_dim_monitor import run_dim_stream, slice_xy  # noqa: E402
from stream_dgps import make_order_graph_stream  # noqa: E402


def _tables(kind, seed=0):
    return make_order_graph_stream(
        n_ref=200,
        n_new=80,
        n_batches=4,
        onset_batch=1,
        n_merchants=12,
        n_users=30,
        kind=kind,
        seed=seed,
    )


class RecsysDimTests(unittest.TestCase):
    def test_slice_drops_y_and_region(self):
        tables = _tables("covariate_south", seed=1)
        from graph_fsds_localize import slice_stream

        win = slice_stream(tables, 0)["new"]
        X, names = slice_xy(win, "all")
        self.assertEqual(X.shape[1], len(names))
        self.assertNotIn("Y", names)
        self.assertNotIn("region", names)
        self.assertIn("amount", names)
        self.assertIn("merchant_gmv", names)
        self.assertIn("user_tenure", names)

    def test_covariate_order_recovers_amount(self):
        tables = _tables("covariate_south", seed=2)
        out = run_dim_stream(tables, "order", seed=2, with_po=False)
        self.assertFalse(out["y_in_X"])
        self.assertIn("amount", out["fsds_recovered"] + [r["feature"] for r in out["fsds"][:2]])
        self.assertGreater(out["rows"][-1]["mmd"], out["rows"][0]["mmd"])

    def test_user_grain_does_not_claim_gmv(self):
        tables = _tables("covariate_south", seed=3)
        out = run_dim_stream(tables, "user", seed=3, with_po=False)
        self.assertNotIn("merchant_gmv", out["names"])
        self.assertNotIn("amount", out["names"])
        self.assertFalse(out["y_in_X"])

    def test_posthoc_south_mmd_beats_north_on_order_covariate(self):
        tables = _tables("covariate_south", seed=4)
        out = run_dim_stream(tables, "order", seed=4, with_po=False)
        s = out["localization"]["south"].get("mmd") or 0.0
        n = out["localization"]["north"].get("mmd") or 0.0
        self.assertGreater(s, n)


if __name__ == "__main__":
    unittest.main()
