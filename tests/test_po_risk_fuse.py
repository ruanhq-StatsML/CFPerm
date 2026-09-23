"""Tests for long⊗short PO-risk fusion (concept modality emphasis)."""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def _bootstrap():
    if "agod" not in sys.modules or not getattr(sys.modules["agod"], "__path__", None):
        pkg = types.ModuleType("agod")
        pkg.__path__ = [str(_ROOT / "agod")]
        sys.modules["agod"] = pkg

    def load(name: str, rel: str):
        if name in sys.modules and hasattr(sys.modules[name], "fuse_long_short"):
            return sys.modules[name]
        if name in sys.modules and name != "agod.po_risk_train":
            return sys.modules[name]
        spec = importlib.util.spec_from_file_location(name, _ROOT / rel)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        assert spec.loader is not None
        spec.loader.exec_module(mod)
        return mod

    load("agod.lr_controller", "agod/lr_controller.py")
    load("agod.next_step", "agod/next_step.py")
    load("agod.smooth_router", "agod/smooth_router.py")
    return load("agod.po_risk_train", "agod/po_risk_train.py")


_mod = _bootstrap()
fuse_long_short = _mod.fuse_long_short
metric_to_alpha = _mod.metric_to_alpha
next_step_actuators = _mod.next_step_actuators
next_step_actuators_fused = _mod.next_step_actuators_fused


def test_fuse_long_short_spike_tilts_short_weight():
    mods = ["a", "b", "c"]
    out = fuse_long_short(
        mods,
        po={"a": 1.05, "b": 3.5, "c": 1.0},
        po_prev={"a": 1.0, "b": 1.0, "c": 1.0},
        po_ema={"a": 1.0, "b": 1.0, "c": 1.0},
        proto={"a": 0, "b": 0, "c": 0},
    )
    assert out["top_spike_mod"] == "b"
    assert out["omega_short"] >= out["omega_long"]
    assert out["alpha"]["b"] == max(out["alpha"].values())


def test_fuse_calm_prefers_long():
    out = fuse_long_short(
        ["a", "b"],
        po={"a": 2.0, "b": 0.5},
        po_prev={"a": 2.0, "b": 0.5},
        po_ema={"a": 2.0, "b": 0.5},
    )
    assert out["spike_ratio"] < 0.25
    assert out["top_concept_mod"] == "a"
    assert abs(out["omega_long"] - 0.55) < 1e-6


def test_po_fuse_metric_and_hierarchical_actuators():
    mods = ["img", "txt", "aud"]
    rep = metric_to_alpha(
        "po_fuse",
        mods,
        po={"img": 0.2, "txt": 2.0, "aud": 0.3},
        po_prev={"img": 0.2, "txt": 0.4, "aud": 0.3},
        po_ema={"img": 0.5, "txt": 1.2, "aud": 0.4},
        mmd={"img": 0.0, "txt": 0.0, "aud": 0.0},
        proto={"img": 0.0, "txt": 0.2, "aud": 0.0},
    )
    assert "freeze_hint" in rep["diag"]
    act = next_step_actuators_fused(rep, mods)
    assert sum(act["step_alloc"].values()) >= 1
    assert act["flops_rel"] <= 1.0
    assert act["step_alloc"]["txt"] >= max(act["step_alloc"].values()) * 0.4


def test_compare_path_uses_fused_for_po_fuse_only():
    """Mirror compare runner: po_fuse → fused actuators; others → plain."""
    mods = ["a", "b", "c"]
    po = {"a": 1.0, "b": 2.5, "c": 0.4}
    pack = metric_to_alpha(
        "po_fuse",
        mods,
        po=po,
        po_prev={"a": 1.0, "b": 1.0, "c": 0.4},
        po_ema={"a": 1.0, "b": 1.2, "c": 0.5},
    )
    fused = next_step_actuators_fused(pack, mods)
    assert "step_alloc" in fused and "freeze_mask" in fused
    assert pack["diag"].get("top_spike_mod") == "b"
    soft = metric_to_alpha("po_soft", mods, po=po)
    act_soft = next_step_actuators(soft["alpha"], mods)
    assert sum(act_soft["step_alloc"].values()) >= 1
