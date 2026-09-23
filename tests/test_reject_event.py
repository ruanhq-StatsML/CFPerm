"""Reject-event resolver tests (no torch)."""
from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def _boot():
    if "agod" not in sys.modules or not getattr(sys.modules["agod"], "__path__", None):
        pkg = types.ModuleType("agod")
        pkg.__path__ = [str(_ROOT / "agod")]
        sys.modules["agod"] = pkg

    def load(name: str, rel: str):
        if name in sys.modules:
            if name == "agod.po_risk_train" and hasattr(
                sys.modules[name], "stream_reject_proxy"
            ):
                return sys.modules[name]
            if name == "agod.reject_event" and hasattr(
                sys.modules[name], "resolve_stream_reject"
            ):
                return sys.modules[name]
            if name not in ("agod.po_risk_train", "agod.reject_event"):
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
    load("agod.po_risk_train", "agod/po_risk_train.py")
    return load("agod.reject_event", "agod/reject_event.py")


_mod = _boot()


def test_resolve_priority_external_then_hop_then_proxy():
    resolve = _mod.resolve_stream_reject
    hop_fires = _mod.hop_fires
    assert not hop_fires(0.2, None)
    assert hop_fires(0.3, 0.1, gate=1.5)
    assert not hop_fires(0.12, 0.1, gate=1.5)

    ext = resolve(external_rejected=True)
    assert ext["rejected"] and ext["source"] == "external_rfperm"

    quiet_ext = resolve(external_rejected=False, e_now=0.9, e_prev=0.1)
    assert quiet_ext["rejected"] is False

    hop = resolve(e_now=0.4, e_prev=0.2, oos_gate=1.5, use_proxy_fallback=False)
    assert hop["rejected"] and hop["source"] == "hop_oos"

    prox = resolve(
        e_now=0.1,
        e_prev=0.1,
        po_mods={"a": 0.2, "b": 0.2},
        mods=["a", "b"],
        use_proxy_fallback=True,
    )
    assert prox["source"] == "proxy"
    assert prox["rejected"] is True
