"""Two-dataset stream form → OnlineRFPerm changepoint."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load_mod():
    spec = importlib.util.spec_from_file_location(
        "llm_stream_cp", ROOT / "scripts" / "agod" / "llm_stream_changepoint.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _synth_halu(n: int = 600) -> list[dict]:
    rows = []
    for i in range(n):
        hall = i % 3 == 0
        rows.append(
            {
                "question": f"Q{i} what is fact {i % 7}?",
                "answer": (
                    f"Wrong invented claim {i}." if hall else f"Correct fact {i % 7}."
                ),
                "hallucination": "yes" if hall else "no",
            }
        )
    return rows


def _synth_hh(n_pairs: int = 200) -> list[dict]:
    rows = []
    for i in range(n_pairs):
        rows.append(
            {
                "chosen": (
                    f"\n\nHuman: Q{i} topic {i % 5}?\n\n"
                    f"Assistant: Therefore I recommend option {i} carefully."
                ),
                "rejected": (
                    f"\n\nHuman: Q{i} topic {i % 5}?\n\n"
                    f"Assistant: maybe idk lol {i}!!!"
                ),
            }
        )
    return rows


def test_stream_schema_and_rfperm_fires():
    mod = _load_mod()
    cut = 3

    halu = mod.stream_from_halueval(
        _synth_halu(600), emb_dim=32, n_per=50, cut_batch=cut, seed=0
    )
    assert set(halu[0]) >= {"t", "question", "answer", "embedding", "score"}
    assert halu[0]["embedding"].shape[1] == 1
    assert halu[0]["embedding"].ndim == 2
    assert all(r["t"] == i for i, r in enumerate(halu))

    hh = mod.stream_from_hh(
        _synth_hh(240), emb_dim=32, n_per=40, cut_batch=cut, seed=1
    )
    assert hh[0]["embedding"].shape[0] == 10 + 32 + 1  # style10 + hash + pref

    for stream, seed in ((halu, 0), (hh, 1)):
        hops = mod.online_rfperm_changepoints(stream, gate=1.15, seed=seed)
        assert len(hops) >= 2
        fires = [h for h in hops if h["fired"]]
        assert fires, "expected at least one OnlineRFPerm fire after cut"
        # first strong fire should not be far before the cut
        assert fires[0]["batch_t"] >= cut - 1


@pytest.mark.skipif(
    not (ROOT / "results" / "agod" / "llm_changepoint" / "summary.json").exists(),
    reason="run llm_stream_changepoint.py first",
)
def test_cached_artifacts_roundtrip():
    summary = json.loads(
        (ROOT / "results" / "agod" / "llm_changepoint" / "summary.json").read_text()
    )
    assert summary["method_changepoint"] == "OnlineRFPerm"
    assert summary["method_po_risk_recommended"].startswith("BOCPD")
    assert summary["form"] == ["t", "question", "answer", "embedding(p,1)", "score"]
    assert summary["halu_fires"] >= 1
    assert summary["hh_fires"] >= 1
    emb = np.load(ROOT / "results" / "agod" / "llm_changepoint" / "halu_stream.npy")
    assert emb.ndim == 2 and emb.shape[0] > 0
