"""OnlineRFPerm on continuous ``param.grad.norm()`` time series.

Same continuous-time skeleton as ``agod.online_rfperm`` (rank / EWMA p +
online FDR), but the scalar ``T_t`` is a *layer gradient energy* deviation
from a reference regime — not serving MSE.

Use case (Sep17 MVP dashboard): earlier layer-local reject to guide freeze /
back-prop depth before MSE / PO / MMD fully break.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch
import torch.nn as nn

from agod.online_rfperm import online_fdr_step, rank_pvalue


@dataclass
class GradRFPermState:
    """Per-layer continuous-time OnlineRFPerm state on grad norms."""

    name: str
    e_ref: float = 0.0
    T_hist: List[float] = field(default_factory=list)
    p_hist: List[float] = field(default_factory=list)
    reject_hist: List[int] = field(default_factory=list)
    g_hist: List[float] = field(default_factory=list)  # raw grad norms
    wealth: float = 1.0
    n_burn: int = 0


def layer_grad_norms(model: nn.Module) -> Dict[str, float]:
    """``{layer_name: ||grad||_2}`` over parameters that currently have grad."""
    out: Dict[str, float] = {}
    for name, module in model.named_modules():
        if name == "":
            continue
        grads = [p.grad.detach().float().reshape(-1) for p in module.parameters(recurse=False) if p.grad is not None]
        if not grads:
            continue
        flat = torch.cat(grads)
        out[name] = float(torch.linalg.vector_norm(flat).item())
    return out


def relative_grad_shares(norms: Dict[str, float], eps: float = 1e-12) -> Dict[str, float]:
    """Share of total grad energy per layer (scale-robust OnlineRFPerm input)."""
    tot = float(sum(norms.values())) + eps
    return {k: v / tot for k, v in norms.items()}


def init_grad_rfperm(
    layer_names: Sequence[str],
    *,
    e_ref: Optional[Dict[str, float]] = None,
) -> Dict[str, GradRFPermState]:
    e_ref = e_ref or {}
    return {
        name: GradRFPermState(name=name, e_ref=float(e_ref.get(name, 0.0)))
        for name in layer_names
    }


def update_grad_rfperm(
    state: GradRFPermState,
    g: float,
    *,
    burn_in: bool = False,
    alpha: float = 0.05,
    ewma: bool = True,
    fdr: str = "alpha_investing",
    use_relative_to_ref: bool = True,
) -> dict:
    """One continuous-time step on scalar grad energy ``g``.

    ``T = g - e_ref`` (default) so large unexpected gradients → large T → small p
    under the same orientation as MSE-OnlineRFPerm.
    """
    g = float(g)
    state.g_hist.append(g)
    if use_relative_to_ref:
        T = g - float(state.e_ref)
    else:
        T = g
    out = {"g": g, "T": T, "p": 1.0, "reject": False, "burn_in": burn_in, "layer": state.name}
    if burn_in:
        state.T_hist.append(T)
        state.n_burn += 1
        state.p_hist.append(1.0)
        state.reject_hist.append(0)
        # update e_ref as running mean of burn-in g (regime baseline)
        state.e_ref = float(np.mean(state.g_hist))
        return out
    p = rank_pvalue(T, state.T_hist, ewma=ewma)
    # online_fdr_step expects OnlineRFPermState-like wealth / reject_hist
    rej = online_fdr_step(state, p, alpha=alpha, procedure=fdr)  # type: ignore[arg-type]
    state.T_hist.append(T)
    state.p_hist.append(p)
    out.update({"p": p, "reject": bool(rej)})
    return out


def first_reject_index(reject_hist: Sequence[int], *, after: int = 0) -> Optional[int]:
    """First reject time index with ``t >= after``; None if never."""
    for t, r in enumerate(reject_hist):
        if t >= after and int(r) == 1:
            return int(t)
    return None


def lead_time(grad_reject_t: Optional[int], mse_break_t: Optional[int]) -> Optional[int]:
    """``t_grad - t_mse``; negative ⇒ Grad earlier than MSE break."""
    if grad_reject_t is None or mse_break_t is None:
        return None
    return int(grad_reject_t - mse_break_t)


def earliest_layer_reject(
    states: Dict[str, GradRFPermState],
    *,
    after: int = 0,
) -> Dict[str, Optional[int]]:
    """Per-layer and global first reject indices (``t >= after``)."""
    per = {name: first_reject_index(st.reject_hist, after=after) for name, st in states.items()}
    times = [t for t in per.values() if t is not None]
    per["__any__"] = int(min(times)) if times else None
    return per
