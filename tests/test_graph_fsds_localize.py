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
    assemble_feature_library,
    assert_not_outcome,
    freeze_ref_stats,
    fsds_select,
    gated_share,
    join_profile,
    run_pipeline,
    slice_stream,
    unify_to_order,
)
from order_graph_nx import subset_scan  # noqa: E402
from stream_dgps import FEATURE_LIBRARY, make_order_graph_stream  # noqa: E402


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
        # Frozen gmv is the D_ref serving lookup, not the current (amount-driven) gmv.
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

        return graph_shift_cuts(
            cut["ref"], cut["new"], stats, seed=seed, min_n=6, with_graph_contrast=True
        )

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


class LevelSetTests(unittest.TestCase):
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

        return graph_shift_cuts(
            cut["ref"], cut["new"], stats, seed=seed, min_n=6, with_graph_contrast=False
        )

    def test_cut_name_is_level_set_not_louvain(self):
        from order_graph_nx import CUT_LEVEL_SET

        self.assertIn("level_set", CUT_LEVEL_SET)
        self.assertNotIn("louvain", CUT_LEVEL_SET)
        pack = self._pack("covariate_south", seed=9)
        self.assertEqual(pack["level_set"]["cut_name"], CUT_LEVEL_SET)
        self.assertFalse(pack["with_graph_contrast"])
        self.assertIsNone(pack["bundled"])
        self.assertIsNone(pack["structural"])

    def test_covariate_merchant_level_set_recovers_south(self):
        pack = self._pack("covariate_south", seed=9)
        ls = pack["level_set"]
        self.assertGreater(ls["jaccard_merchant_vs_south"], ls["jaccard_user_vs_south"])
        self.assertGreaterEqual(ls["jaccard_merchant_vs_south"], 0.5)
        self.assertGreaterEqual(ls["merchant"]["south_frac_loud_nodes"], 0.75)
        self.assertGreaterEqual(ls["merchant"]["slices"]["mmd"]["south_frac"], 0.6)

    def test_concept_level_set_uses_y_not_community(self):
        pack = self._pack("concept_south", seed=10)
        mer = pack["level_set"]["merchant"]
        self.assertGreaterEqual(mer["south_frac_loud_nodes"], 0.5)
        self.assertGreaterEqual(mer["slices"]["cmean_y"]["south_frac"], 0.5)
        self.assertGreaterEqual(pack["level_set"]["jaccard_merchant_vs_south"], 0.25)

    def test_coverage_scan_recovers_south_without_louvain(self):
        pack = self._pack("covariate_south", seed=9)
        ls = pack["level_set"]
        mer = ls["merchant"]
        self.assertEqual(mer["coverage"]["cut_name"], "subset_scan/{coverage}")
        self.assertNotIn("louvain", mer["coverage"]["cut_name"])
        self.assertGreaterEqual(ls["jaccard_coverage_merchant_vs_south"], 0.4)
        self.assertGreaterEqual(ls["jaccard_mass_merchant_vs_south"], 0.4)
        self.assertGreater(
            ls["jaccard_coverage_merchant_vs_south"], ls["jaccard_user_vs_south"]
        )

    def test_pipeline_level_set_loud_vs_other(self):
        tables = make_order_graph_stream(
            n_ref=300,
            n_new=120,
            n_batches=3,
            onset_batch=1,
            n_merchants=12,
            n_users=40,
            kind="covariate_south",
            seed=11,
        )
        out = run_pipeline(
            tables,
            t=2,
            grain="order",
            mode="localize",
            subset_by="level_set",
            seed=11,
            min_n=15,
            with_logo=False,
            with_po=False,
        )
        self.assertEqual(out["leakage"]["graph_cut"], "level_set/{phi>=tau}")
        self.assertFalse(out["leakage"]["used_networkx"])
        self.assertFalse(out["graph"]["with_graph_contrast"])
        self.assertIsNone(out["graph"]["bundled"])
        self.assertIn(out["loud_subset"], ("loud", "other"))
        loud = next(p for p in out["portraits"] if p["subset"] == "loud")
        self.assertGreaterEqual(loud["south_frac"], 0.6)
        self.assertGreater(loud["pi_mmd"], 0.5)
        self.assertFalse(out["leakage"]["y_in_Z"])


class SubsetScanTests(unittest.TestCase):
    def _scores(self):
        # Two loud, one tiny-loud, two quiet. Scale is estimated from the stack.
        return {
            0: {"mmd": 0.40, "cmean_x": 1.20, "cmean_y": 0.0, "po": 0.0, "n": 40},
            1: {"mmd": 0.35, "cmean_x": 1.00, "cmean_y": 0.0, "po": 0.0, "n": 40},
            2: {"mmd": 0.30, "cmean_x": 0.90, "cmean_y": 0.0, "po": 0.0, "n": 4},
            3: {"mmd": 0.00, "cmean_x": 0.00, "cmean_y": 0.0, "po": 0.0, "n": 40},
            4: {"mmd": 0.00, "cmean_x": 0.00, "cmean_y": 0.0, "po": 0.0, "n": 40},
        }

    def test_loud_ids_are_a_prefix_not_a_partition(self):
        scores = self._scores()
        out = subset_scan(scores, rule="coverage", floor=0.0, coverage=0.80, weight="phi")
        self.assertNotIn("louvain", out["cut_name"])
        ranked = out["ranked"]
        loud = out["loud_ids"]
        self.assertEqual(loud, ranked[: len(loud)])
        self.assertLess(len(loud), len(ranked))

    def test_mass_weight_prefers_large_bags(self):
        scores = self._scores()
        phi = subset_scan(scores, rule="coverage", floor=0.0, coverage=0.50, weight="phi")
        mass = subset_scan(scores, rule="coverage", floor=0.0, coverage=0.50, weight="mass")
        # Tiny bag 2 can rank high on intensity; mass should put 0 or 1 first.
        self.assertIn(mass["loud_ids"][0], (0, 1))
        self.assertTrue(set(phi["loud_ids"]).issubset(set(phi["ranked"])))


class FeatureLibraryTests(unittest.TestCase):
    def test_catalog_has_no_y_and_covers_three_grains(self):
        names = [r["feature"] for r in FEATURE_LIBRARY]
        grains = {r["grain"] for r in FEATURE_LIBRARY}
        self.assertNotIn("Y", names)
        self.assertNotIn("y", names)
        self.assertEqual(grains, {"order", "merchant", "user"})
        self.assertIn("amount", names)
        self.assertIn("merchant_gmv", names)

    def test_pipeline_library_marks_planted_and_selects_amount(self):
        tables = _stream("covariate_south", seed=12)
        self.assertIn("feature_library", tables["meta"])
        out = run_pipeline(
            tables,
            t=2,
            grain="order",
            mode="localize",
            subset_by="level_set",
            seed=12,
            min_n=15,
            with_logo=False,
            with_po=False,
            with_graph_contrast=False,
        )
        lib = out["library"]
        feats = [r["feature"] for r in lib["rows"]]
        self.assertNotIn("Y", feats)
        amount = next(r for r in lib["rows"] if r["feature"] == "amount")
        hour = next(r for r in lib["rows"] if r["feature"] == "hour")
        self.assertTrue(amount["planted"])
        self.assertEqual(amount["planted_how"], "x")
        self.assertFalse(hour["planted"])
        self.assertIn("amount", lib["selected"])
        self.assertGreater(amount["score"], hour["score"])
        self.assertFalse(out["leakage"].get("used_networkx"))
        rank = out["scan_rank"]
        self.assertGreater(rank["n"], 0)
        self.assertEqual(rank["rows"], sorted(rank["rows"], key=lambda r: (-r["phi"], r["merchant_id"])))

    def test_assemble_rejects_outcome_column(self):
        from stream_dgps import FEATURE_LIBRARY as CATALOG

        fake = {
            "rank": {"order": [{"feature": "Y", "score": 1, "mmd": 0, "cmean_x": 0, "cmean_y": 0, "loud": True}]},
            "selected_names": ["Y"],
        }
        # Catalog does not include Y, so assemble stays clean even if FSDS leaked.
        lib = assemble_feature_library(fake, {"kind": "covariate_south"})
        self.assertEqual([r["feature"] for r in lib["rows"]], [c["feature"] for c in CATALOG])
        self.assertNotIn("Y", lib["selected"])


if __name__ == "__main__":
    unittest.main()
