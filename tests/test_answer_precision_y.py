"""Unit checks for answer-precision vs Jaccard (Hotpot-style dilution)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "agod"))
import online_rfperm_live_infer as live  # noqa: E402


def test_answer_precision_not_diluted_by_long_knowledge():
    short_ans = "Delhi is the head office city"
    long_k = (
        "The Oberoi family is an Indian family famous for hotels. "
        "The Oberoi Group has its head office in Delhi. "
        + ("Distractor paragraph about unrelated history. " * 80)
    )
    prec = live.answer_precision(short_ans, long_k)
    jac = live.jaccard(short_ans, long_k)
    assert prec > 0.8, prec
    assert jac < prec / 2, (jac, prec)
    assert jac < 0.35, jac
    assert prec > jac


def test_binary_y_from_faith_threshold():
    # quiet-style echo → high faith → Y=0 under thr=0.45
    ans = "Arthur's Magazine was an American literary periodical"
    kn = "Arthur's Magazine was an American literary periodical published in Philadelphia"
    faith = live.answer_precision(ans, kn)
    assert faith >= 0.45
    y = 1 if faith < 0.45 else 0
    assert y == 0
    # invent-style → low faith → Y=1
    invent = "I am certain it was founded in Atlantis-42 by Dr. Fabricatus"
    faith_h = live.answer_precision(invent, kn)
    assert faith_h < 0.45
    assert (1 if faith_h < 0.45 else 0) == 1
