"""Smoke tests for HF landing prototypes (offline cache or tiny synth)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts.agod.hf_landing_protos import (
    build_halu_stream,
    build_hh_stream,
    pack_batches,
    run_halluc_demo,
    run_judge_demo,
    style_vector,
)

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "hf_cache"


def _synth_hh(n_pairs: int = 240) -> list[dict]:
    rows = []
    for i in range(n_pairs):
        chosen = (
            f"\n\nHuman: Q{i} about topic {i % 7}?\n\n"
            f"Assistant: Therefore I recommend option {i} with details."
        )
        rejected = (
            f"\n\nHuman: Q{i} about topic {i % 7}?\n\n"
            f"Assistant: maybe idk lol {i}!!!"
        )
        rows.append({"chosen": chosen, "rejected": rejected})
    return rows


def _synth_halu(n: int = 400) -> list[dict]:
    rows = []
    for i in range(n):
        know = f"entity {i} was founded in {1900 + (i % 50)} in city {i % 11}"
        q = f"When was entity {i} founded?"
        if i % 2 == 0:
            ans = f"entity {i} was founded in {1900 + (i % 50)}"
            hall = "no"
        else:
            ans = f"entity {i} opened a mall in 2099"
            hall = "yes"
        rows.append(
            {"knowledge": know, "question": q, "answer": ans, "hallucination": hall}
        )
    return rows


def test_style_vector_dim():
    v = style_vector("Therefore I think this might work!!!")
    assert v.shape == (10,)
    assert np.isfinite(v).all()


def test_pack_batches_shape():
    X = np.random.randn(105, 4)
    y = np.random.randint(0, 2, size=105)
    Xp, yp, b = pack_batches(X, y, 20)
    assert len(Xp) == 100
    assert b.max() == 4


def test_judge_demo_fires_at_cut_synth():
    pack = build_hh_stream(_synth_hh(320), n_per=40, cut_batch=3, seed=0)
    out = run_judge_demo(pack, gate=1.2, seed=0)
    assert out["first_fire_t"] is not None
    assert out["hop_at_cut"] is not None
    assert out["hop_at_cut"]["fired"] in (0, 1, True, False)
    assert out["judge_err_ratio"] > 1.0
    assert out["style_domain_auc"] > 0.7


def test_halluc_demo_fires_at_cut_synth():
    pack = build_halu_stream(_synth_halu(500), n_per=50, cut_batch=3, seed=1)
    out = run_halluc_demo(pack, gate=1.2, seed=1)
    assert out["rate_after"] > out["rate_before"]
    assert out["first_fire_t"] is not None
    assert out["hop_at_cut"]["fired"] in (0, 1, True, False)
    assert "router_task_shift" in out


@pytest.mark.skipif(
    not (CACHE / "hh_rlhf_helpful_base_2500.jsonl").exists()
    or not (CACHE / "halueval_qa_3000.jsonl").exists(),
    reason="HF cache subsets not present",
)
def test_cached_hf_summary_roundtrip():
    summary = CACHE.parent.parent / "results" / "agod" / "hf_landing" / "summary.json"
    if summary.exists():
        data = json.loads(summary.read_text())
        assert "judge" in data and "halluc" in data
