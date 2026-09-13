"""Named LR policies for the Amazon MMD use-case (and clones)."""
from __future__ import annotations

from typing import Any, Mapping, Sequence

from .lr_controller import intensity_gain, softmax_scores
from .shift import decompose_hybrid, decompose_mmd, decompose_rf

POLICIES = ("B1", "B2", "B3", "B4", "B5")


def select_policy(
    policy: str,
    *,
    msg: Any,
    blocks0: Mapping[str, Any],
    blocks1: Mapping[str, Any],
    y0,
    y1,
    mods: Sequence[str],
    tau: float,
    kappa: float,
    seed: int,
) -> tuple[dict[str, float], dict[str, Any], dict[str, float]]:
    """Return (raw_alpha, decomp, gain) for one online window.

    B1 equal | B2 RF | B3 MMD | B4 MMD+gain | B5 MMD-cov+PO (FSDS)
    """
    rf_d = decompose_rf(msg, mods)
    mmd_d = decompose_mmd(blocks0, blocks1, y0, y1, mods, seed=seed)
    unit = {m: 1.0 for m in mods}

    if policy == "B1":
        raw = {m: 1.0 / len(mods) for m in mods}
        return raw, mmd_d, unit
    if policy == "B2":
        return softmax_scores(rf_d["score"], mods, tau), rf_d, unit
    if policy == "B3":
        return softmax_scores(mmd_d["score"], mods, tau), mmd_d, unit
    if policy == "B4":
        return (
            softmax_scores(mmd_d["score"], mods, tau),
            mmd_d,
            intensity_gain(mmd_d, mods, kappa=kappa),
        )
    if policy == "B5":
        hy = decompose_hybrid(msg, mmd_d, mods)
        return softmax_scores(hy["score"], mods, tau), hy, unit
    raise ValueError(f"unknown policy {policy!r}; expected one of {POLICIES}")
