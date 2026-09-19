#!/usr/bin/python3
"""Recsys grain slices: no Y in X; FSDS recovers planted order columns."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from recsys_dim_monitor import (  # noqa: E402
    kendall_vs_planted,
    planted_in_grain,
    run_dim_stream,
    slice_xy,
)
from stream_dgps import (  # noqa: E402
    GRAIN_FEATS,
    MERCHANT_FEATS,
    ORDER_FEATS,
    USER_FEATS,
    causal_batch_ema,
    causal_ema_lookup,
    causal_rolling_share,
    make_order_graph_stream,
)


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

    def test_localization_emits_mmd_cmean_po(self):
        tables = _tables("covariate_south", seed=6)
        out = run_dim_stream(tables, "order", seed=6, with_po=True)
        s = out["localization"]["south"]
        for key in ("mmd", "cmean_x", "cmean_y", "po"):
            self.assertIn(key, s)
        self.assertIsNotNone(s["mmd"])
        self.assertIsNotNone(s["cmean_x"])
        self.assertIsNotNone(s["cmean_y"])
        if out["localization"]["n_south"] >= 40:
            self.assertIsNotNone(s["po"])

    def test_planted_stays_inside_grain(self):
        names = ("amount", "hour", "n_items", "channel")
        planted = planted_in_grain("covariate_south", names)
        self.assertIn("amount", planted)
        self.assertIn("channel", planted)
        self.assertNotIn("merchant_gmv", planted)

    def test_kendall_recovers_perfect_order(self):
        rank = [
            {"feature": "amount", "score": 3.0},
            {"feature": "channel", "score": 2.0},
            {"feature": "hour", "score": 0.1},
        ]
        tau = kendall_vs_planted(rank, {"amount": 1.0, "channel": 0.7}, score_key="score")
        self.assertGreater(tau, 0.5)

    def test_rfperm_mse_vimp_keys(self):
        tables = _tables("covariate_south", seed=5)
        out = run_dim_stream(tables, "order", seed=5, with_po=False)
        self.assertIn("amount", out["names"])
        self.assertTrue(out["mse_vimp"]["top"])
        self.assertIn(out["addis_status"], ("hit", "FAR", "miss"))


class GrainComputeLeakageTests(unittest.TestCase):
    def test_grain_catalogs_are_disjoint(self):
        o, m, u = set(ORDER_FEATS), set(MERCHANT_FEATS), set(USER_FEATS)
        self.assertFalse(o & m)
        self.assertFalse(o & u)
        self.assertFalse(m & u)
        self.assertEqual(tuple(GRAIN_FEATS["order"]), ORDER_FEATS)

    def test_ema_excludes_current_and_future(self):
        ids = np.array([0, 0, 1, 0])
        vals = np.array([10.0, 20.0, 99.0, 40.0])
        prior = np.array([1.0, 2.0])
        out = causal_ema_lookup(ids, vals, 2, prior, lam=1.0)
        np.testing.assert_allclose(out, [1.0, 10.0, 2.0, 20.0])
        later = vals.copy()
        later[-1] = 1e6
        out2 = causal_ema_lookup(ids, later, 2, prior, lam=1.0)
        np.testing.assert_allclose(out2[:-1], out[:-1])

    def test_gmv_ignores_y_and_current_amount(self):
        tables = _tables("covariate_south", seed=11)
        gmv = tables["X_merchant"][:, list(tables["names_merchant"]).index("merchant_gmv")]
        amount = tables["X_order"][:, list(tables["names_order"]).index("amount")]
        mid = tables["merchant_id"]
        meta = tables["meta"]
        gmv2 = causal_batch_ema(
            mid,
            amount,
            int(meta["n_merchants"]),
            tables["merchant_table"]["X"][:, list(MERCHANT_FEATS).index("merchant_gmv")],
            int(meta["n_ref"]),
            int(meta["n_new"]),
            int(meta["n_batches"]),
        )
        np.testing.assert_allclose(gmv, gmv2)
        n_ref = int(meta["n_ref"])
        n_new = int(meta["n_new"])
        last = slice(n_ref + n_new, n_ref + 2 * n_new)
        # Same merchant, same snapshot inside a new batch.
        m0 = int(mid[last][0])
        same = gmv[last][mid[last] == m0]
        self.assertGreater(len(same), 1)
        np.testing.assert_allclose(same, same[0])
        # Flipping this batch's amounts must not change this batch's GMV.
        amt2 = amount.copy()
        amt2[last] = amt2[last] + 50.0
        gmv_flip = causal_batch_ema(
            mid,
            amt2,
            int(meta["n_merchants"]),
            tables["merchant_table"]["X"][:, list(MERCHANT_FEATS).index("merchant_gmv")],
            int(meta["n_ref"]),
            int(meta["n_new"]),
            int(meta["n_batches"]),
        )
        np.testing.assert_allclose(gmv_flip[last], gmv[last])

    def test_hist_freq_is_past_share_not_y(self):
        tables = _tables("covariate_south", seed=12)
        freq = tables["X_user"][:, list(tables["names_user"]).index("user_hist_freq")]
        y = tables["Y"]
        self.assertLess(abs(float(np.corrcoef(freq, y)[0, 1])), 0.35)
        uid = tables["user_id"]
        # Same user, same frozen snapshot at every row.
        u0 = int(uid[0])
        same = freq[uid == u0]
        np.testing.assert_allclose(same, same[0])

    def test_slice_drops_cross_grain_columns(self):
        tables = _tables("covariate_south", seed=1)
        from graph_fsds_localize import slice_stream

        win = slice_stream(tables, 0)["new"]
        _, order_n = slice_xy(win, "order")
        _, merch_n = slice_xy(win, "merchant")
        _, user_n = slice_xy(win, "user")
        self.assertEqual(order_n, ORDER_FEATS)
        self.assertNotIn("merchant_gmv", order_n)
        self.assertNotIn("amount", merch_n)
        self.assertNotIn("amount", user_n)
        self.assertNotIn("Y", order_n + merch_n + user_n)


if __name__ == "__main__":
    unittest.main()
