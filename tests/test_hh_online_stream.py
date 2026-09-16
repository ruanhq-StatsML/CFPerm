"""HH tidy stream + continuous OnlineRFPerm detection."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


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


def test_tidy_schema_and_continuous_detect():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "hh_stream", ROOT / "scripts" / "agod" / "hh_online_rfperm_stream.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)

    df, feat, meta = mod.build_tidy_table(
        _synth_hh(240), n_per=40, cut_batch=3, featurizer="style_hash", seed=0
    )
    assert list(df.columns)[:5] == [
        "t_idx",
        "batch",
        "question",
        "answer",
        "y",
    ]
    assert "feature" in df.columns
    assert feat.shape[0] == len(df)
    assert meta["featurizer"].startswith("style10")
    # preference hop after cut
    assert df.loc[df["batch"] >= 3, "flipped"].mean() > 0.8

    hist = mod.continuous_online_rfperm(
        feat,
        df["y"].to_numpy(),
        df["batch"].to_numpy(),
        df["t_idx"].to_numpy(),
        gate=1.15,
        topk=5,
    )
    assert len(hist) >= 3
    assert any(h["fired"] for h in hist)
    # continuous: e_prev of hop t equals e_now of hop t-1
    for i in range(1, len(hist)):
        assert hist[i]["e_prev"] == pytest.approx(hist[i - 1]["e_now"], rel=1e-9)


@pytest.mark.skipif(
    not (ROOT / "data" / "hf_cache" / "hh_rlhf_helpful_base_2500.jsonl").exists(),
    reason="HH cache missing",
)
def test_cached_hh_stream_artifacts_roundtrip():
    summary = ROOT / "results" / "agod" / "hh_online_stream" / "summary.json"
    if not summary.exists():
        pytest.skip("run hh_online_rfperm_stream.py first")
    data = json.loads(summary.read_text())
    assert data["columns"] == ["t_idx", "batch", "question", "answer", "y", "feature"]
    assert data["n_fires"] >= 1
