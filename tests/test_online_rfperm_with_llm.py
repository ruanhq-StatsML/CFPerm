#!/usr/bin/python3
"""One sample: onlinePermOOB_with_LLM on a small df. Y is the last column."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from online_rfperm_with_llm import (  # noqa: E402
    make_rf_adapter,
    make_stationary_df,
    make_tabpfn_adapter,
    onlinePermOOB_with_LLM,
)


class OnlinePermOOBWithLLMTests(unittest.TestCase):
    def test_onlinePermOOB_with_LLM_y_is_last_column(self):
        df = make_stationary_df(n=200, p=6, seed=1)
        rec = onlinePermOOB_with_LLM(
            df,
            llm_adapter=make_tabpfn_adapter(seed=1),
            dl_adapter=make_rf_adapter(),
            llm_cols=[3, 4, 5],
            dl_cols=[0, 1, 2],
            ref_batch_size=80,
            batch_size=20,
            burnin=1,
            seed=1,
        )
        self.assertEqual(df.shape[1], 7)
        self.assertIn("addis_1", rec)
        self.assertEqual(len(rec["MSE_list"]), 6)
        self.assertTrue(np.all(np.isfinite(rec["MSE_list"])))
        self.assertEqual(rec["n_llm_cols"], 3)
        self.assertEqual(rec["n_dl_cols"], 3)


if __name__ == "__main__":
    unittest.main()
