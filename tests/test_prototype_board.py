"""Prototype board: declared FSDS layout and one inspectable hop."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from prototype_board import inspect_shapes, layout_spec  # noqa: E402
from typed_shift_stepsize import make_typed_stream  # noqa: E402


def test_layout_spec_is_768_512_768():
    spec = layout_spec()
    assert spec["constants"]["P_X"] == 2048
    assert spec["constants"]["P_s"] == 2049
    assert spec["groups"]["video"] == {"start": 0, "stop": 768, "dim": 768}
    assert spec["groups"]["audio"]["dim"] == 512
    assert spec["groups"]["text"]["dim"] == 768


def test_inspect_shapes_one_hop():
    stream = make_typed_stream(n_batches=4, n_per=24, n_classes=4, seed=1, cov={"video": 0.2})
    rec = inspect_shapes(stream, t=1)
    assert rec["stream"]["X"] == [96, 2048]
    assert rec["blocks"]["video"]["X_curr"] == [24, 768]
    assert rec["blocks"]["video"]["mu_old"] == [4, 768]
    assert rec["blocks"]["video"]["instdisc"]["V"] == [24, 32]
    assert rec["blocks"]["video"]["cos_sdc"] > rec["blocks"]["video"]["cos_stale"]
