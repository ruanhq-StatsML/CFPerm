"""Synthetic continuous-gains smoke tests (no torch / Affec)."""
from __future__ import annotations

import importlib.util
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def _load_smoke():
    spec = importlib.util.spec_from_file_location(
        "smoke_po_boost_synthetic", _ROOT / "scripts/smoke_po_boost_synthetic.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def test_synthetic_po_fuse_can_ship():
    smoke = _load_smoke()
    payload = smoke.run_synthetic_compare(seed=0)
    fuse = payload["versions"]["po_fuse"]
    assert fuse["mean_flops_rel"] < 1.0
    assert fuse["t_to_acc_star"] is not None
    assert fuse["ship_pass"] is True
    assert fuse["reject_sources"].get("hop_oos", 0) >= 1
    equal = payload["versions"]["equal"]
    assert equal["mean_flops_rel"] == 1.0
