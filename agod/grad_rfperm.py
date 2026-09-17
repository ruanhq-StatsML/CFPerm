"""OnlineRFPerm on continuous ``param.grad.norm()`` — *one* scalar stream.

Engineering口径 (no cross-parameter multiple testing):

1. Collect grads only on **unfrozen** params (``requires_grad=True`` and
   ``grad is not None``).
2. Form **one** scalar
   ``g_t = ||∇_{θ_U} L||_2 = sqrt(Σ_i ||∇_{θ_i} L||_2²)``
   (ℓ₂ pool of the full unfrozen grad vector — *not* a mean of per-layer
   norms, and *not* a separate OnlineRFPerm per layer).
3. Run a **single** continuous-time OnlineRFPerm on ``T_t = g_t - e_ref``
   (rank / EWMA p + online FDR). One hypothesis stream ⇒ no multiplicity
   across layers/parameters.
4. Per-layer relative shares ``||g_ℓ|| / g_t`` are **diagnostics only**
   (freeze-depth ranking after a global reject). They do **not** get their
   own reject / FDR.

Use case (Sep17 MVP): earlier global Grad reject vs MSE / PO / MMD, then
optional layer shares for back-prop depth.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn as nn

from agod.online_rfperm import online_fdr_step, rank_pvalue


@dataclass
class GradRFPermState:
    """Single-stream continuous-time OnlineRFPerm on unfrozen grad energy."""

    name: str = "unfrozen_l2"
    e_ref: float = 0.0
    T_hist: List[float] = field(default_factory=list)
    p_hist: List[float] = field(default_factory=list)
    reject_hist: List[int] = field(default_factory=list)
    g_hist: List[float] = field(default_factory=list)
    wealth: float = 1.0
    n_burn: int = 0


def iter_unfrozen_grads(model: nn.Module) -> Iterable[Tuple[str, torch.Tensor]]:
    """``(param_name, grad)`` for params that are trainable and have grad."""
    for name, p in model.named_parameters():
        if not p.requires_grad:
            continue
        if p.grad is None:
            continue
        yield name, p.grad.detach().float()


def unfrozen_grad_l2(model: nn.Module) -> float:
    """Global ``||∇_{θ_U} L||_2`` over all unfrozen parameters (one scalar)."""
    chunks = [g.reshape(-1) for _, g in iter_unfrozen_grads(model)]
    if not chunks:
        return 0.0
    return float(torch.linalg.vector_norm(torch.cat(chunks)).item())


def layer_grad_norms(model: nn.Module, *, unfrozen_only: bool = True) -> Dict[str, float]:
    """Per-module ``||grad||_2`` — diagnostic / share only, not for FDR reject.

    When ``unfrozen_only``, skips modules whose direct params are all frozen
    or have no grad.
    """
    out: Dict[str, float] = {}
    for name, module in model.named_modules():
        if name == "":
            continue
        grads = []
        for p in module.parameters(recurse=False):
            if unfrozen_only and not p.requires_grad:
                continue
            if p.grad is None:
                continue
            grads.append(p.grad.detach().float().reshape(-1))
        if not grads:
            continue
        out[name] = float(torch.linalg.vector_norm(torch.cat(grads)).item())
    return out


def relative_grad_shares(norms: Dict[str, float], eps: float = 1e-12) -> Dict[str, float]:
    """Share of total grad energy per layer (diagnostic; not a test statistic)."""
    tot = float(sum(norms.values())) + eps
    return {k: v / tot for k, v in norms.items()}


def param_grad_norms(model: nn.Module, *, unfrozen_only: bool = True) -> Dict[str, float]:
    """Per-parameter ``||grad||_2`` (diagnostic only)."""
    out: Dict[str, float] = {}
    for name, p in model.named_parameters():
        if unfrozen_only and not p.requires_grad:
            continue
        if p.grad is None:
            continue
        out[name] = float(torch.linalg.vector_norm(p.grad.detach().float()).item())
    return out


def mean_normalized_param_norms(
    norms: Dict[str, float],
    *,
    ref: Optional[Dict[str, float]] = None,
    eps: float = 1e-12,
) -> float:
    """Optional alternate scalar: mean of per-param norms / burn-in refs.

    **Not** the default gate. Prefer ``unfrozen_grad_l2``. Documented so we
    do not silently average-and-threshold without saying so.
    """
    if not norms:
        return 0.0
    if ref is None:
        return float(np.mean(list(norms.values())))
    vals = [norms[k] / (float(ref.get(k, 1.0)) + eps) for k in norms]
    return float(np.mean(vals))


def init_grad_rfperm(
    name: str = "unfrozen_l2",
    *,
    e_ref: float = 0.0,
) -> GradRFPermState:
    """One OnlineRFPerm state for the global unfrozen grad scalar."""
    return GradRFPermState(name=name, e_ref=float(e_ref))


def init_grad_rfperm_layers(
    layer_names: Sequence[str],
    *,
    e_ref: Optional[Dict[str, float]] = None,
) -> Dict[str, GradRFPermState]:
    """Deprecated multi-stream init — kept for diagnostics / ablations only.

    Do **not** take ``any(reject)`` as a gate (multiple testing). Use
    ``init_grad_rfperm`` for the production gate.
    """
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
    out = {"g": g, "T": T, "p": 1.0, "reject": False, "burn_in": burn_in, "stream": state.name}
    if burn_in:
        state.T_hist.append(T)
        state.n_burn += 1
        state.p_hist.append(1.0)
        state.reject_hist.append(0)
        state.e_ref = float(np.mean(state.g_hist))
        return out
    p = rank_pvalue(T, state.T_hist, ewma=ewma)
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
    """Per-stream first reject (ablation / diagnostic only — not the gate)."""
    per = {name: first_reject_index(st.reject_hist, after=after) for name, st in states.items()}
    times = [t for t in per.values() if t is not None]
    per["__any__"] = int(min(times)) if times else None
    return per


def top_share_layers(shares: Dict[str, float], *, k: int = 3) -> List[Tuple[str, float]]:
    """Rank layers by relative grad share (freeze-depth hint after global reject)."""
    return sorted(shares.items(), key=lambda kv: -kv[1])[:k]
