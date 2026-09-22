#!/usr/bin/env python3
"""Multi-step audit: one row per assistant hop, chosen is not Y."""
from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts.build_llm_audit_multistep_xy import (
    MIN_TURNS,
    STEP_COLS,
    assistant_turns,
    explode_conversations,
    overlay_hop,
    pack_hops,
    prototype_auditor,
    write_multistep,
)
from scripts.llm_audit_online_bootstrap_prototype import CUT_BATCH


CONV_TWO = """
Human: hello

Assistant: Sure, I can help with that plan in some detail so we get it right.

Human: more please

Assistant: Here is a longer follow-up with several sentences about the next step and how to proceed carefully.
"""

CONV_ONE = """
Human: hi

Assistant: Ok.
"""


class LlmAuditMultistepTests(unittest.TestCase):
    def test_splits_assistant_hops_and_drops_singletons(self):
        self.assertGreaterEqual(len(assistant_turns(CONV_TWO)), MIN_TURNS)
        self.assertEqual(len(assistant_turns(CONV_ONE)), 1)
        rows = explode_conversations([CONV_TWO, CONV_ONE])
        self.assertTrue(all(int(r["step"]) >= 0 for r in rows))
        episodes = {r["episode"] for r in rows}
        self.assertEqual(len(episodes), 1)
        self.assertGreaterEqual(len(rows), 2)

    def test_auditor_is_not_chosen_and_uses_this_hop(self):
        stub = prototype_auditor("No.")
        long = prototype_auditor(
            "Here is a careful, reasonably long reply that actually answers the question "
            "with enough content to ship as a complete turn rather than a stub."
        )
        self.assertIn(stub, (0, 1))
        self.assertIn(long, (0, 1))
        self.assertGreaterEqual(long, stub)

    def test_pack_and_write_schema_has_no_text_or_chosen(self):
        rows = explode_conversations([CONV_TWO] * 8)
        packed = pack_hops(rows, n=8, n_per=4)
        self.assertEqual(len(packed), 8)
        self.assertEqual({r["batch"] for r in packed}, {0, 1})
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "xy.csv"
            write_multistep(path, packed)
            with path.open() as f:
                header = next(csv.reader(f))
            keys = {k.lower() for k in header}
            self.assertEqual(header[:4], ["y", "batch", "episode", "step"])
            self.assertEqual(header[4:], STEP_COLS)
            self.assertIn("x_step", keys)
            self.assertIn("x_refuse", keys)
            for bad in ("chosen", "rejected", "prompt", "text", "conversation"):
                self.assertNotIn(bad, keys)

    def test_hop_overlay_flips_after_cut_keeps_x(self):
        rows = explode_conversations([CONV_TWO] * 20)
        packed = pack_hops(rows, n=16, n_per=4)
        hopped = overlay_hop(packed, cut_batch=2, seed=0)
        y0 = np.array([r["y"] for r in packed])
        y1 = np.array([r["y"] for r in hopped])
        batch = np.array([r["batch"] for r in packed])
        self.assertTrue(np.array_equal(
            np.array([r["x_n_toks"] for r in packed]),
            np.array([r["x_n_toks"] for r in hopped]),
        ))
        self.assertTrue(np.array_equal(y0[batch < 2], y1[batch < 2]))
        self.assertFalse(np.array_equal(y0[batch >= 2], y1[batch >= 2]))
        self.assertEqual(CUT_BATCH, 4)


if __name__ == "__main__":
    unittest.main()
