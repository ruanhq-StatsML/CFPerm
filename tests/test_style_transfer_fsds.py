#!/usr/bin/python3
"""One sample: FSDS+LOGO routing on style-transfer tokens. Y is not a feature."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from style_transfer_fsds import (  # noqa: E402
    Config,
    feature_names,
    fsds_route,
    make_data,
    serving_xy,
    style_groups,
)


class StyleTransferFsdsTests(unittest.TestCase):
    def test_fsds_route_noise_walk_lights_mmd(self):
        cfg = Config(n_attr=48, seed=7)
        exist = make_data(48, cfg, seed=3, drift=False)
        new = make_data(48, cfg, seed=4, drift=True)
        X, Y = serving_xy(exist)
        names = feature_names(cfg)

        self.assertEqual(X.shape, (48, cfg.max_len))
        self.assertEqual(Y.shape, (48,))
        self.assertNotIn("y", {n.lower() for n in names})
        self.assertNotIn("s_tgt", names)
        self.assertEqual(X.shape[1], exist["src"].shape[1])

        route = fsds_route(exist, new, cfg, seed=7, with_logo=True)

        self.assertFalse(route["y_in_x"])
        self.assertEqual(route["groups"]["content"], [0, 4])
        self.assertEqual(route["groups"]["style_src"], [4, 8])
        self.assertIn("alpha", route)
        self.assertIn("plan", route)
        self.assertGreater(route["metrics"]["mmd"], route["quiet"]["mmd"])
        self.assertTrue(route["mmd_broken"] or route["cov"] > 0.0)

        top_names = [r["feature"] for r in route["fsds_top"]]
        self.assertTrue(any(n.startswith("noise_") for n in top_names))
        self.assertTrue(all(not n.lower().startswith("y") for n in top_names))

        plan = route["plan"]
        self.assertIsNotNone(plan)
        self.assertIn("towers", plan)
        for g in style_groups(cfg):
            self.assertIn(g, plan["towers"])
        noise_pi = route["logo"]["pi_mmd"]["noise"]
        self.assertGreaterEqual(noise_pi, route["logo"]["pi_mmd"]["content"])


if __name__ == "__main__":
    unittest.main()
