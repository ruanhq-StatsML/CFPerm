#!/usr/bin/python3
"""One sample: onlinePermOOB_reject. Last column is Y; NaN means unshown."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from online_rfperm_reject import (  # noqa: E402
    make_reject_df,
    onlinePermOOB_reject,
    split_xy_s,
)


class OnlinePermOOBRejectTests(unittest.TestCase):
    def test_onlinePermOOB_reject_nan_y_is_unshown(self):
        df = make_reject_df(n=320, p=5, kind="select", onset=160, seed=2)
        X, Y, S = split_xy_s(df)
        self.assertEqual(df.shape[1], 6)
        self.assertTrue(np.any(~np.isfinite(Y)))
        self.assertEqual(int(S.sum()), int(np.isfinite(Y).sum()))
        rec = onlinePermOOB_reject(df, ref_batch_size=120, batch_size=40, seed=2)
        self.assertFalse(rec["y_in_x"])
        self.assertEqual(len(rec["MSE_cc_list"]), 5)
        self.assertEqual(len(rec["MSE_ips_list"]), 5)
        self.assertEqual(len(rec["Brier_sel_list"]), 5)
        self.assertIn("sel_1", rec)
        self.assertTrue(np.all(np.isfinite(X)))


if __name__ == "__main__":
    unittest.main()
