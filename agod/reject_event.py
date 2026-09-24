"""Stream reject event resolution (R3 hook for OnlineRFPerm / hop OOS).

Priority
--------
1. ``external_rejected`` — real OnlineRFPerm / CFPerm flag when wired
2. consecutive OOS hop: ``e_now / e_prev ≥ gate`` (``hop_fires``)
3. ``stream_reject_proxy`` — PO/MMD/ΔPO thresholds (compare fallback)

IPTW still goes through ``po_iptw_weights(..., rejected=)`` either way.
"""
from __future__ import annotations

from typing import Any, Mapping, Sequence

from agod.po_risk_train import stream_reject_proxy


def shift_ratio(e_now: float, e_prev: float) -> float:
    return float(e_now) / (float(e_prev) + 1e-8)


def hop_fires(e_now: float | None, e_prev: float | None, gate: float = 1.5) -> bool:
    """Consecutive OOS gate (OnlineRFPerm-style). First hop quiet."""
    if e_now is None or e_prev is None:
        return False
    ratio = shift_ratio(e_now, e_prev)
    return bool(np_isfinite(ratio) and ratio >= float(gate))


def np_isfinite(x: float) -> bool:
    import math

    return bool(math.isfinite(float(x)))


def resolve_stream_reject(
    *,
    external_rejected: bool | None = None,
    e_now: float | None = None,
    e_prev: float | None = None,
    oos_gate: float = 1.5,
    po_mods: Mapping[str, float] | None = None,
    mmd_mods: Mapping[str, float] | None = None,
    po_prev: Mapping[str, float] | None = None,
    mods: Sequence[str] | None = None,
    use_proxy_fallback: bool = True,
) -> dict[str, Any]:
    """Unified reject decision for the compare / post-train stream."""
    if external_rejected is not None:
        return {
            "rejected": bool(external_rejected),
            "source": "external_rfperm",
            "reasons": ["external"] if external_rejected else [],
            "ratio": None,
            "proxy": False,
        }
    if hop_fires(e_now, e_prev, gate=oos_gate):
        return {
            "rejected": True,
            "source": "hop_oos",
            "reasons": ["e_now/e_prev"],
            "ratio": shift_ratio(float(e_now), float(e_prev)),  # type: ignore[arg-type]
            "e_now": float(e_now),  # type: ignore[arg-type]
            "e_prev": float(e_prev),  # type: ignore[arg-type]
            "proxy": False,
        }
    if use_proxy_fallback and po_mods is not None:
        prox = stream_reject_proxy(
            po_mods=po_mods,
            mmd_mods=mmd_mods,
            po_prev=po_prev,
            mods=mods,
        )
        prox = dict(prox)
        prox["source"] = "proxy"
        prox["ratio"] = None
        return prox
    return {
        "rejected": False,
        "source": "none",
        "reasons": [],
        "ratio": None,
        "proxy": False,
    }
