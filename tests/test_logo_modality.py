#!/usr/bin/env python3
"""LOGO modality shares are localization, not a unique decomp."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from dl_model_registry import AnyMLP, apply_train_stem, apply_train_top_i  # noqa: E402
from logo_modality import (  # noqa: E402
    TOWER_STEM,
    TOWER_TOP,
    as_groups,
    cumulative_regret,
    drop_group,
    logo_batch,
    plan_next_batch,
    subset_excess,
    two_layer_ratios,
)
from stream_dgps import TRIMODAL_GROUPS, make_trimodal_stream  # noqa: E402
from streaming_po_risk import ACTION_FREEZE, ACTION_KEEP, ACTION_XSHIFT  # noqa: E402


class LogoUnitTests(unittest.TestCase):
    def test_drop_group_keeps_other_blocks(self):
        X = np.arange(24, dtype=float).reshape(2, 12)
        g = as_groups(TRIMODAL_GROUPS)
        Xa = drop_group(X, g, "audio")
        self.assertEqual(Xa.shape, (2, 8))
        np.testing.assert_array_equal(Xa[0, :4], X[0, :4])
        np.testing.assert_array_equal(Xa[0, 4:], X[0, 8:])

    def test_two_layer_ratios_sum_and_mix(self):
        r = two_layer_ratios(
            {"video": 0.5, "audio": 0.5, "text": 0.0},
            {"video": 1.0, "audio": 0.0, "text": 0.0},
            {"video": 0.0, "audio": 1.0, "text": 0.0},
        )
        self.assertAlmostEqual(r["video"]["mix_po"], 1.0)
        self.assertAlmostEqual(r["audio"]["mix_mmd"], 1.0)
        self.assertTrue(r["text"]["quiet"])

    def test_plan_freeze_trains_top_on_po_loud_tower(self):
        ratios = two_layer_ratios(
            {"video": 1.0, "audio": 0.0, "text": 0.0},
            {"video": 1.0, "audio": 0.0, "text": 0.0},
            {"video": 0.0, "audio": 0.0, "text": 0.0},
        )
        plan = plan_next_batch(ACTION_FREEZE, ratios)
        self.assertEqual(plan["towers"]["video"]["tower"], TOWER_TOP)
        self.assertEqual(plan["update"], "selective_freeze")

    def test_plan_xshift_stems_mmd_loud_tower(self):
        ratios = two_layer_ratios(
            {"video": 0.0, "audio": 1.0, "text": 0.0},
            {"video": 0.0, "audio": 0.0, "text": 0.0},
            {"video": 0.0, "audio": 1.0, "text": 0.0},
        )
        plan = plan_next_batch(ACTION_XSHIFT, ratios)
        self.assertEqual(plan["towers"]["audio"]["tower"], TOWER_STEM)
        self.assertEqual(plan["fusion"], "fusion_infer")

    def test_keep_is_full_train(self):
        ratios = two_layer_ratios(
            {"video": 0.0, "audio": 0.0, "text": 0.0},
            {"video": 0.0, "audio": 0.0, "text": 0.0},
            {"video": 0.0, "audio": 0.0, "text": 0.0},
        )
        plan = plan_next_batch(ACTION_KEEP, ratios)
        self.assertEqual(plan["update"], "full_train")
        self.assertTrue(all(v["tower"] == "full_train" for v in plan["towers"].values()))

    def test_apply_train_stem_only_bottom(self):
        m = AnyMLP(4, 1, (8, 4), dropout=0.0)
        apply_train_stem(m)
        bottom = all(p.requires_grad for p in m.blocks["fc0"].parameters())
        head = all(not p.requires_grad for p in m.head.parameters())
        self.assertTrue(bottom)
        self.assertTrue(head)
        apply_train_top_i(m, 1)
        self.assertTrue(all(p.requires_grad for p in m.head.parameters()))

    def test_regret_cumsum(self):
        r = cumulative_regret([0.4, 0.5, 0.6], [0.3, 0.3, 0.3])
        np.testing.assert_allclose(r, [0.1, 0.3, 0.6])


class LogoDgpTests(unittest.TestCase):
    def _split(self, kind, seed=0):
        X, Y, sl, meta = make_trimodal_stream(
            n_ref=240,
            n_new=80,
            n_batches=3,
            onset_batch=1,
            kind=kind,
            seed=seed,
        )
        n_ref = meta["n_ref"]
        n_new = meta["n_new"]
        t = 2
        lo = n_ref + t * n_new
        hi = lo + n_new
        return X[:n_ref], Y[:n_ref], X[lo:hi], Y[lo:hi], sl[lo:hi], sl[:n_ref], meta

    def test_concept_video_po_share_not_on_text(self):
        Xr, Yr, Xn, Yn, *_ = self._split("concept_video")
        out = logo_batch(Xr, Yr, Xn, Yn, TRIMODAL_GROUPS, seed=1)
        self.assertLess(out["pi_po"]["text"], out["pi_po"]["video"] + 1e-12)

    def test_covariate_audio_mmd_share_leads(self):
        Xr, Yr, Xn, Yn, *_ = self._split("covariate_audio", seed=3)
        out = logo_batch(Xr, Yr, Xn, Yn, TRIMODAL_GROUPS, seed=3)
        self.assertGreaterEqual(out["pi_mmd"]["audio"], out["pi_mmd"]["video"])
        self.assertGreaterEqual(out["pi_mmd"]["audio"], out["pi_mmd"]["text"])

    def test_subset_mmd_peaks_on_planted_slice(self):
        Xr, Yr, Xn, Yn, sl_n, sl_r, meta = self._split("covariate_audio", seed=4)
        rows = subset_excess(Xr, Yr, Xn, Yn, sl_n, labels_ref=sl_r, seed=4, min_n=8)
        self.assertTrue(rows)
        top = max(rows, key=lambda r: r["mmd"])
        self.assertEqual(int(top["subset"]), int(meta["shifted_slice"]))


if __name__ == "__main__":
    unittest.main()
