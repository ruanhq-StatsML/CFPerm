"""Significant-batch-only metrics for OnlineRFPerm-gated adaptation.

Non-reject steps all use uniform weights for gated modes, so including them
dilutes the comparison. Evaluate next-batch MSE only where the gate opens.
"""
from __future__ import annotations

from typing import Dict, Iterable, Mapping, Optional, Sequence

import numpy as np


def align_gate_to_mse(
    gate_on: Sequence[int] | np.ndarray,
    mse_next: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Truncate gate flags to steps that have a held-out next-batch MSE."""
    g = np.asarray(gate_on, dtype=bool)
    m = np.asarray(mse_next, dtype=float)
    return g[: len(m)]


def pick_gate_mask(
    results: Mapping[str, Mapping],
    preferred_modes: Iterable[str] = (),
) -> Optional[np.ndarray]:
    """Return boolean gate mask aligned to ``mse_next`` length, or None."""
    modes = list(preferred_modes) + list(results.keys())
    seen = set()
    for mode in modes:
        if mode in seen or mode not in results:
            continue
        seen.add(mode)
        r = results[mode]
        if "gate_on" not in r:
            continue
        mse = r.get("mse_next")
        if not mse:
            continue
        return align_gate_to_mse(r["gate_on"], mse)
    return None


def significant_only_stats(
    results: Mapping[str, Mapping],
    *,
    gate_mask: Optional[np.ndarray] = None,
    preferred_gate_modes: Iterable[str] = (),
) -> Dict[str, dict]:
    """Per-mode mean/std/cum MSE on OnlineRFPerm-significant batches only."""
    mask = gate_mask
    if mask is None:
        mask = pick_gate_mask(results, preferred_gate_modes)
    if mask is None:
        raise ValueError("no gate_on found in results")

    n_sig = int(mask.sum())
    out: Dict[str, dict] = {
        "_meta": {
            "n_significant": n_sig,
            "n_eval": int(len(mask)),
            "duty": float(mask.mean()) if len(mask) else 0.0,
        }
    }
    base_u = np.asarray(results["uniform"]["mse_next"], float)[: len(mask)]
    base_d = np.asarray(results["dre"]["mse_next"], float)[: len(mask)] if "dre" in results else None

    for mode, r in results.items():
        m = np.asarray(r["mse_next"], float)[: len(mask)]
        if n_sig == 0:
            mean = std = cum = float("nan")
            rel_u = rel_d = float("nan")
        else:
            sel = m[mask]
            mean = float(sel.mean())
            std = float(sel.std())
            cum = float(sel.sum())
            rel_u = float(mean / (base_u[mask].mean() + 1e-12))
            if base_d is not None:
                rel_d = float(mean / (base_d[mask].mean() + 1e-12))
            else:
                rel_d = float("nan")
        out[mode] = {
            "mse_mean_sig": mean,
            "mse_std_sig": std,
            "cum_mse_sig": cum,
            "rel_mse_vs_uniform_sig": rel_u,
            "rel_mse_vs_dre_sig": rel_d,
        }
    return out


def annotate_results_with_sig(
    results: Dict[str, dict],
    *,
    preferred_gate_modes: Iterable[str] = (),
) -> Dict[str, dict]:
    """Mutate/return results dict with ``mse_mean_sig`` fields + shared gate meta."""
    stats = significant_only_stats(results, preferred_gate_modes=preferred_gate_modes)
    meta = stats.pop("_meta")
    for mode, r in results.items():
        r["mse_mean_sig"] = stats[mode]["mse_mean_sig"]
        r["mse_std_sig"] = stats[mode]["mse_std_sig"]
        r["cum_mse_sig"] = stats[mode]["cum_mse_sig"]
        r["rel_mse_vs_uniform_sig"] = stats[mode]["rel_mse_vs_uniform_sig"]
        r["n_significant"] = meta["n_significant"]
        r["sig_duty"] = meta["duty"]
    return stats
