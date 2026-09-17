#!/usr/bin/python3
"""Graph + FSDS localization: no Y leakage; subset vs other on three readouts."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from graph_fsds_localize import (  # noqa: E402
    assert_not_outcome,
    freeze_ref_stats,
    fsds_select,
    gated_share,
    join_profile,
    run_pipeline,
    slice_stream,
    unify_to_order,
)
from stream_dgps import make_order_graph_stream  # noqa: E402


def _stream(kind, seed=0):
    return make_order_graph_stream(
        n_ref=240,
        n_new=80,
        n_batches=3,
        onset_batch=1,
        kind=kind,
        seed=seed,
    )


class LeakageTests(unittest.TestCase):
    def test_outcome_name_rejected(self):
        with self.assertRaises(ValueError):
            assert_not_outcome(["amount", "Y"])

    def test_gated_share_ignores_mmd_jitter(self):
        quiet = gated_share({"south": 0.011, "north": 0.010}, abs_floor=0.02)
        self.assertEqual(quiet["south"], 0.0)
        self.assertEqual(quiet["north"], 0.0)
        loud = gated_share({"south": 0.4, "north": 0.05}, abs_floor=0.02)
        self.assertGreater(loud["south"], loud["north"])

    def test_serve_profile_ignores_new_batch_y_and_x(self):
        tables = _stream("covariate_south", seed=2)
        cut = slice_stream(tables, t=2)
        stats = freeze_ref_stats(cut["ref"], seed=2)
        joined = join_profile(cut["new"]["merchant_id"], stats["merchant_profile"])
        shifted = cut["new"]["X_merchant"]
        # Frozen gmv is not the current (shifted) gmv on south orders.
        south = cut["new"]["region"] == "south"
        gmv_j = list(cut["ref"]["names_merchant"]).index("merchant_gmv")
        self.assertGreater(
            abs(float(shifted[south, gmv_j].mean()) - float(joined[south, gmv_j].mean())),
            0.5,
        )
        self.assertNotIn("Y", stats["merchant_profile"]["names"])


class FsdsUnifyTests(unittest.TestCase):
    def test_covariate_selects_amount_not_y(self):
        tables = _stream("covariate_south", seed=3)
        cut = slice_stream(tables, t=2)
        stats = freeze_ref_stats(cut["ref"], seed=3)
        fsds = fsds_select(cut["ref"], cut["new"], stats, top_k=3, seed=3)
        names = fsds["selected_names"]
        self.assertIn("amount", names)
        self.assertNotIn("Y", names)
        self.assertNotIn("y", names)

    def test_unify_rejects_empty_and_keeps_y_out(self):
        tables = _stream("covariate_south", seed=4)
        cut = slice_stream(tables, t=2)
        stats = freeze_ref_stats(cut["ref"], seed=4)
        fsds = fsds_select(cut["ref"], cut["new"], stats, top_k=3, seed=4)
        uni = unify_to_order(cut["new"], fsds["selected"], stats, mode="localize")
        self.assertNotIn("Y", uni["names"])
        self.assertEqual(len(uni["Y"]), len(uni["Z"]))
        with self.assertRaises(ValueError):
            unify_to_order(cut["new"], {"order": [], "merchant": [], "user": []}, stats)


class PipelineTests(unittest.TestCase):
    def test_covariate_south_leads_mmd_and_cmean_x(self):
        tables = _stream("covariate_south", seed=5)
        out = run_pipeline(
            tables,
            t=2,
            grain="order",
            mode="localize",
            subset_by="region",
            seed=5,
            min_n=15,
            with_logo=False,
            with_po=False,
        )
        self.assertFalse(out["leakage"]["y_in_Z"])
        self.assertIn("amount", out["fsds"]["selected_names"])
        self.assertEqual(out["loud_subset"], "south")
        south = next(p for p in out["portraits"] if p["subset"] == "south")
        north = next(p for p in out["portraits"] if p["subset"] == "north")
        self.assertGreater(south["pi_mmd"], north["pi_mmd"])
        self.assertGreater(south["cmean_x"], north["cmean_x"])
        gap = south["gap_vs_other"]
        self.assertIsNotNone(gap["mmd"])
        self.assertGreater(gap["mmd"], 0.0)

    def test_concept_south_leads_cmean_y(self):
        tables = _stream("concept_south", seed=6)
        out = run_pipeline(
            tables,
            t=2,
            grain="order",
            mode="localize",
            subset_by="region",
            seed=6,
            min_n=15,
            with_logo=False,
            with_po=False,
        )
        south = next(p for p in out["portraits"] if p["subset"] == "south")
        north = next(p for p in out["portraits"] if p["subset"] == "north")
        self.assertGreater(abs(south["cmean_y"]), abs(north["cmean_y"]))
        self.assertEqual(out["loud_subset"], "south")
        # Concept keeps P(X) closer: south MMD share should not dominate like covariate.
        self.assertLess(south["mmd"], 2.5)

    def test_merchant_grain_still_tags_south(self):
        tables = _stream("covariate_south", seed=7)
        out = run_pipeline(
            tables,
            t=2,
            grain="merchant",
            mode="localize",
            subset_by="region",
            seed=7,
            min_n=1,
            with_logo=False,
            with_po=False,
        )
        self.assertEqual(out["loud_subset"], "south")
        self.assertFalse(out["leakage"]["y_in_Z"])


class BundledGraphTests(unittest.TestCase):
    def _pack(self, kind, seed, n_merchants=12):
        tables = make_order_graph_stream(
            n_ref=300,
            n_new=120,
            n_batches=3,
            onset_batch=1,
            n_merchants=n_merchants,
            n_users=40,
            kind=kind,
            seed=seed,
        )
        cut = slice_stream(tables, t=2)
        stats = freeze_ref_stats(cut["ref"], seed=seed)
        from graph_fsds_localize import graph_shift_cuts

        return graph_shift_cuts(cut["ref"], cut["new"], stats, seed=seed, min_n=6)

    def test_package_is_networkx_not_pyg(self):
        from order_graph_nx import GRAPH_PACKAGE, CUT_BUNDLED

        self.assertEqual(GRAPH_PACKAGE, "networkx")
        self.assertIn("louvain", CUT_BUNDLED)

    def test_y_never_enters_graph_attrs(self):
        pack = self._pack("covariate_south", seed=8)
        self.assertFalse(pack["bundled"]["y_in_graph"])
        self.assertFalse(pack["y_in_structural"])

    def test_covariate_bundled_cut_recovers_south(self):
        pack = self._pack("covariate_south", seed=9)
        self.assertIsNotNone(pack["bundled"]["loud_community"])
        self.assertGreaterEqual(pack["bundled"]["loud_south_frac"], 0.75)

    def test_concept_bundled_cut_uses_cmean_not_structure(self):
        pack = self._pack("concept_south", seed=10)
        self.assertGreaterEqual(pack["bundled"]["loud_south_frac"], 0.6)

    def test_own_ref_flags_change_full_ref_flags_heterogeneity(self):
        pack = self._pack("covariate_south", seed=9)
        s = pack["own_vs_full"]
        self.assertGreater(s["south_own_mmd"], s["north_own_mmd"])
        self.assertGreater(s["south_own_cmean_x"], s["north_own_cmean_x"])
        own_gap = s["south_own_mmd"] - s["north_own_mmd"]
        full_gap = s["south_full_mmd"] - s["north_full_mmd"]
        self.assertGreater(own_gap, full_gap)

    def test_layers_evaluated_after_lift_to_orders(self):
        pack = self._pack("covariate_south", seed=9)
        layers = pack["layers"]
        self.assertGreater(layers["jaccard_merchant_vs_south"], layers["jaccard_user_vs_south"])
        self.assertGreaterEqual(layers["jaccard_merchant_vs_south"], 0.4)
        self.assertFalse(layers["merchant"]["y_in_graph"])
        self.assertFalse(layers["user"]["y_in_graph"])


if __name__ == "__main__":
    unittest.main()
