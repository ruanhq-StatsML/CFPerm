#!/usr/bin/python3
"""One sample: FSDS+LOGO routing on style-transfer tokens. Y is not a feature."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from style_transfer_fsds import (  # noqa: E402
    Config,
    feature_names,
    fsds_route,
    make_data,
    make_lexicon,
    pack_style_df,
    serving_xy,
    styleTransferFSDS,
    style_groups,
)


class StyleTransferFsdsTests(unittest.TestCase):
    def test_fsds_route_noise_walk_lights_mmd(self):
        cfg = Config(n_attr=64, seed=7)
        lex = make_lexicon(cfg, seed=0)
        exist = make_data(64, cfg, seed=3, lexicon=lex)
        cov = make_data(64, cfg, seed=4, drift=True, lexicon=lex)
        con = make_data(64, cfg, seed=5, concept=True, lexicon=lex)
        X, Y = serving_xy(exist)
        names = feature_names(cfg)

        self.assertEqual(X.shape, (64, cfg.max_len))
        self.assertEqual(Y.shape, (64,))
        self.assertNotIn("y", {n.lower() for n in names})
        self.assertNotIn("s_tgt", names)
        self.assertEqual(X.shape[1], exist["src"].shape[1])

        route = fsds_route(exist, cov, cfg, seed=7, with_logo=True)
        self.assertFalse(route["y_in_x"])
        self.assertEqual(route["groups"]["content"], [0, 1, 2, 3])
        self.assertEqual(route["groups"]["style_src"], [4, 5, 6, 7])
        self.assertGreater(route["metrics"]["mmd"], route["quiet"]["mmd"])
        self.assertEqual(route["loud"]["group"], "noise")
        self.assertEqual(route["loud"]["kind"], "mmd")
        top_names = [r["feature"] for r in route["fsds_top"]]
        self.assertTrue(all(n.startswith("noise_") for n in top_names))
        self.assertNotEqual(route["plan"]["towers"]["noise"]["tower"], "freeze_tower")

        route_c = fsds_route(exist, con, cfg, seed=7, with_logo=True)
        self.assertEqual(route_c["loud"]["group"], "style_src")
        self.assertEqual(route_c["loud"]["kind"], "po")
        self.assertTrue(all(n.startswith("style_src_") for n in [r["feature"] for r in route_c["fsds_top"]]))
        self.assertEqual(route_c["plan"]["towers"]["style_src"]["tower"], "train_top")

        df = np.vstack([pack_style_df(exist), pack_style_df(cov)])
        recs = styleTransferFSDS(
            df,
            style_groups(cfg),
            ref_batch_size=64,
            batch_size=64,
            seed=7,
            names=names,
            cfg=cfg,
        )
        self.assertEqual(df.shape[1], cfg.max_len + 1)
        self.assertEqual(len(recs), 1)
        self.assertFalse(recs[0]["y_in_x"])
        self.assertEqual(recs[0]["loud_group"], "noise")

        quiet = make_data(64, cfg, seed=6, lexicon=lex)
        route_q = fsds_route(exist, quiet, cfg, seed=7, with_logo=True)
        self.assertIsNone(route_q["loud"])
        self.assertIsNone((route_q.get("plan") or {}).get("overlay"))


if __name__ == "__main__":
    unittest.main()
