#!/usr/bin/env python3
"""Schema lock for manuscript LLM-audit prediction tables."""
from __future__ import annotations

import csv
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "manuscript" / "llm_audit"
X_COLS = [
    "x_n_toks",
    "x_n_chars",
    "x_avg_word",
    "x_qmark",
    "x_bang",
    "x_hedge",
    "x_formal",
    "x_i_count",
    "x_newlines",
    "x_upper",
    "x_refuse",
    "x_please",
    "x_thank",
]


class LlmAuditXyTests(unittest.TestCase):
    def _rows(self, name: str):
        path = OUT / name
        self.assertTrue(path.exists(), msg=path)
        with path.open() as f:
            rows = list(csv.DictReader(f))
        return rows

    def test_prototype_and_real_label_schema(self):
        for name in (
            "xy_hh_helpful_consistent.csv",
            "xy_hh_harmless_consistent.csv",
            "xy_beavertails.csv",
            "xy_wildguard.csv",
            "xy_toxicchat.csv",
        ):
            rows = self._rows(name)
            self.assertEqual(len(rows), 1200, name)
            self.assertEqual(list(rows[0].keys())[:2], ["y", "batch"])
            for c in X_COLS:
                self.assertIn(c, rows[0], msg=f"{name} missing {c}")
            ys = {int(r["y"]) for r in rows}
            self.assertTrue(ys <= {0, 1}, name)

    def test_hh_chosen_is_not_a_column(self):
        for name in ("xy_hh_helpful_consistent.csv", "xy_hh_two_stream.csv"):
            rows = self._rows(name)
            keys = {k.lower() for k in rows[0]}
            self.assertNotIn("chosen", keys)
            self.assertNotIn("rejected", keys)

    def test_cfperm_two_stream(self):
        for name, n in (("xy_hh_two_stream.csv", 2400), ("xy_real_two_stream.csv", 2400)):
            rows = self._rows(name)
            self.assertEqual(len(rows), n, name)
            ts = {int(r["T"]) for r in rows}
            self.assertEqual(ts, {0, 1}, name)
            self.assertIn("x_refuse", rows[0])


if __name__ == "__main__":
    unittest.main()
