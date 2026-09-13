"""Online stacking + erank / modality-balance auxiliaries.

Prediction-side actuators (complement LR decorr):

1. **Online stacking** — per-modality logits fused by learned simplex
   ``w = softmax(ψ)``, updated online with the window stream.

2. **Erank / modality balance loss** on modality hidden prototypes:
     G[i,j] = cos(h̄_i, h̄_j)
     erank(G) = exp(H(λ/Σλ))
     L_erank  = ReLU(erank_floor − erank)     # push rank up when collapsed
     L_align  = mean_{i<j} max(cos_ij, 0)²  # damp positive collinearity
     L_bal    = KL(w ‖ prior)                 # prior=α if collapsed else uniform

Task loss = CE on stacked logits. Evaluation MSE = Brier
``mean((softmax(z) − onehot(y))²)``.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def softmax_np(x: np.ndarray, axis: int = -1) -> np.ndarray:
    z = x - np.max(x, axis=axis, keepdims=True)
    e = np.exp(z)
    return e / np.clip(e.sum(axis=axis, keepdims=True), 1e-12, None)


def probs_mse(logits: np.ndarray, y: np.ndarray, n_class: int | None = None) -> float:
    """Brier / MSE of softmax probs vs one-hot labels."""
    logits = np.asarray(logits, float)
    y = np.asarray(y, int).ravel()
    if logits.ndim != 2 or len(logits) == 0:
        return float("nan")
    n_class = int(n_class or logits.shape[1])
    p = softmax_np(logits, axis=1)
    oh = np.zeros_like(p)
    oh[np.arange(len(y)), np.clip(y, 0, n_class - 1)] = 1.0
    return float(np.mean((p - oh) ** 2))


def probs_mse_torch(logits: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    p = F.softmax(logits, dim=-1)
    oh = F.one_hot(y.long(), num_classes=logits.size(-1)).float()
    return ((p - oh) ** 2).mean()


def gram_erank_torch(
    vectors: Sequence[torch.Tensor],
) -> tuple[torch.Tensor, torch.Tensor]:
    """Effective rank of cosine Gram + mean positive-align² over prototypes."""
    mods = list(vectors)
    m = len(mods)
    if m == 0:
        z = torch.tensor(0.0)
        return z, z
    vs = [F.normalize(v.reshape(-1).float(), dim=0) for v in mods]
    g = torch.eye(m, device=vs[0].device, dtype=vs[0].dtype)
    pos = []
    for i in range(m):
        for j in range(i + 1, m):
            c = torch.clamp(torch.dot(vs[i], vs[j]), -1.0, 1.0)
            g[i, j] = g[j, i] = c
            pos.append(torch.relu(c) ** 2)
    evals = torch.linalg.eigvalsh(0.5 * (g + g.T))
    evals = torch.clamp(evals, min=0.0)
    s = evals.sum().clamp_min(1e-12)
    p = evals / s
    p = p[p > 1e-12]
    h = -(p * torch.log(p)).sum()
    erank = torch.exp(h)
    align = torch.stack(pos).mean() if pos else g.new_zeros(())
    return erank, align


def alpha_stack_kl(
    stack_w: torch.Tensor,
    alpha: Mapping[str, float],
    mods: Sequence[str],
    *,
    lambda_kl: float = 0.50,
) -> dict[str, torch.Tensor]:
    """Attribution-guided stacking: always pull ``stack_w`` toward α.

    This is the weight-socket into online stacking:
      L_kl = KL(stack_w ‖ α)   with α from FSDS / MSG / VIMP routing.
    Unlike erank-triggered balance, the prior is applied every step.
    """
    mods = list(mods)
    w = stack_w.clamp_min(1e-8)
    w = w / w.sum()
    prior = torch.tensor(
        [max(float(alpha.get(m, 0.0)), 1e-8) for m in mods],
        device=w.device,
        dtype=w.dtype,
    )
    prior = prior / prior.sum()
    # KL(w || α) — stack mass tracks attribution proportions
    kl = (w * (torch.log(w) - torch.log(prior))).sum()
    loss = float(lambda_kl) * kl
    return {"loss": loss, "kl": kl, "prior": prior.detach(), "w": w.detach()}


def erank_balance_loss(
    hiddens: Mapping[str, torch.Tensor],
    mods: Sequence[str],
    *,
    stack_w: torch.Tensor | None = None,
    alpha: Mapping[str, float] | None = None,
    erank_floor: float | None = None,
    lambda_align: float = 1.0,
    lambda_erank: float = 1.0,
    lambda_bal: float = 0.25,
) -> dict[str, torch.Tensor]:
    """Erank↑ + positive-align↓ + optional stack-weight KL balance."""
    mods = list(mods)
    protos = [hiddens[m].float().mean(dim=0) for m in mods]
    erank, align = gram_erank_torch(protos)
    floor = float(
        erank_floor if erank_floor is not None else max(len(mods) * 0.85, 1.5)
    )
    floor_t = torch.as_tensor(floor, device=erank.device, dtype=erank.dtype)
    l_erank = torch.relu(floor_t - erank)
    l_align = align
    loss = lambda_erank * l_erank + lambda_align * l_align

    l_bal = erank.new_zeros(())
    if stack_w is not None and lambda_bal > 0:
        w = stack_w.clamp_min(1e-8)
        w = w / w.sum()
        if alpha is not None and float(erank.detach()) < floor:
            prior = torch.tensor(
                [max(float(alpha.get(m, 0.0)), 1e-8) for m in mods],
                device=w.device,
                dtype=w.dtype,
            )
            prior = prior / prior.sum()
        else:
            prior = torch.full_like(w, 1.0 / len(mods))
        l_bal = (w * (torch.log(w) - torch.log(prior))).sum()
        loss = loss + lambda_bal * l_bal

    return {
        "loss": loss,
        "erank": erank,
        "align": align,
        "l_erank": l_erank,
        "l_align": l_align,
        "l_bal": l_bal,
    }


class StackFusion(nn.Module):
    """Per-modality towers + online stacking weights over modality logits."""

    def __init__(
        self,
        dims: Mapping[str, int],
        mods: Sequence[str],
        *,
        fuse: int = 128,
        n_class: int = 2,
        temperature: float = 1.0,
    ):
        super().__init__()
        self.mods = list(mods)
        self.n_class = int(n_class)
        self.temperature = float(temperature)
        self.projs = nn.ModuleDict(
            {
                m: nn.Sequential(
                    nn.Linear(int(dims[m]), fuse),
                    nn.ReLU(),
                    nn.Dropout(0.1),
                )
                for m in self.mods
            }
        )
        self.heads = nn.ModuleDict(
            {m: nn.Linear(fuse, self.n_class) for m in self.mods}
        )
        self.stack_logits = nn.Parameter(torch.zeros(len(self.mods)))

    def encode(self, batch: Mapping[str, torch.Tensor]) -> dict[str, torch.Tensor]:
        return {m: self.projs[m](batch[m]) for m in self.mods}

    def modality_logits(
        self, hiddens: Mapping[str, torch.Tensor]
    ) -> dict[str, torch.Tensor]:
        return {m: self.heads[m](hiddens[m]) for m in self.mods}

    def stack_weights(self, temperature: float | None = None) -> torch.Tensor:
        tau = float(self.temperature if temperature is None else temperature)
        tau = max(tau, 1e-4)
        return F.softmax(self.stack_logits / tau, dim=0)

    def alpha_tensor(
        self,
        alpha: Mapping[str, float],
        *,
        device: torch.device | None = None,
        dtype: torch.dtype | None = None,
    ) -> torch.Tensor:
        device = device or self.stack_logits.device
        dtype = dtype or self.stack_logits.dtype
        w = torch.tensor(
            [max(float(alpha.get(m, 0.0)), 1e-8) for m in self.mods],
            device=device,
            dtype=dtype,
        )
        return w / w.sum()

    def fuse_logits(
        self,
        mod_logits: Mapping[str, torch.Tensor],
        w: torch.Tensor | None = None,
    ) -> torch.Tensor:
        w = self.stack_weights() if w is None else w
        out = None
        for i, m in enumerate(self.mods):
            term = w[i] * mod_logits[m]
            out = term if out is None else out + term
        assert out is not None
        return out

    def forward(
        self,
        batch: Mapping[str, torch.Tensor],
        *,
        return_parts: bool = False,
        fixed_w: Mapping[str, float] | torch.Tensor | None = None,
    ):
        h = self.encode(batch)
        logits_m = self.modality_logits(h)
        if fixed_w is None:
            w = self.stack_weights()
        elif torch.is_tensor(fixed_w):
            w = fixed_w.clamp_min(1e-8)
            w = w / w.sum()
        else:
            w = self.alpha_tensor(fixed_w)
        logits = self.fuse_logits(logits_m, w)
        if return_parts:
            return logits, h, logits_m, w
        return logits


class MeanFusion(nn.Module):
    """Mean-pool fusion baseline (matches gradcos-style Fusion)."""

    def __init__(
        self,
        dims: Mapping[str, int],
        mods: Sequence[str],
        *,
        fuse: int = 128,
        n_class: int = 2,
    ):
        super().__init__()
        self.mods = list(mods)
        self.projs = nn.ModuleDict(
            {
                m: nn.Sequential(
                    nn.Linear(int(dims[m]), fuse),
                    nn.ReLU(),
                    nn.Dropout(0.1),
                )
                for m in self.mods
            }
        )
        self.head = nn.Linear(fuse, n_class)

    def forward(self, batch: Mapping[str, torch.Tensor], *, return_h: bool = False):
        hs = {m: self.projs[m](batch[m]) for m in self.mods}
        logits = self.head(torch.stack([hs[m] for m in self.mods], 0).mean(0))
        if return_h:
            return logits, hs
        return logits


# ---------------------------------------------------------------------------
# Weight socket (prediction-side): plug different priors into online stacking.
# ---------------------------------------------------------------------------

WEIGHT_MODES = (
    "mean_ce",       # mean-pool fusion (no stack_w)
    "stack_ce",      # free w = softmax(psi)
    "stack_alpha",   # CE + λ KL(w || alpha)          attribution prior
    "stack_uniform", # CE + λ KL(w || U)               MoE-style balance
    "stack_fixed",   # freeze w := alpha (towers only)
    "stack_temp",    # sharper softmax(psi / tau), no KL
    "stack_erank",   # erank + align + conditional KL
)

STACK_LEARN_PSI = {
    "stack_ce",
    "stack_alpha",
    "stack_uniform",
    "stack_temp",
    "stack_erank",
}


def uniform_alpha(mods: Sequence[str]) -> dict[str, float]:
    mods = list(mods)
    return {m: 1.0 / len(mods) for m in mods}


def stack_weight_aux(
    mode: str,
    *,
    stack_w: torch.Tensor | None,
    alpha: Mapping[str, float] | None,
    mods: Sequence[str],
    hiddens: Mapping[str, torch.Tensor] | None = None,
    lambda_kl: float = 0.50,
    lambda_bal: float = 0.25,
) -> dict[str, torch.Tensor | float]:
    """Auxiliary loss for a stacking weight mode.

    Returns ``loss`` (0 for modes without aux) and optional diagnostics.
    """
    mods = list(mods)
    zero = (
        stack_w.new_zeros(())
        if isinstance(stack_w, torch.Tensor)
        else torch.tensor(0.0)
    )
    if mode in ("mean_ce", "stack_ce", "stack_fixed", "stack_temp"):
        return {"loss": zero, "kl": float("nan"), "mode": mode}
    if mode == "stack_alpha":
        if stack_w is None or alpha is None:
            return {"loss": zero, "kl": float("nan"), "mode": mode}
        pack = alpha_stack_kl(stack_w, alpha, mods, lambda_kl=lambda_kl)
        return {**pack, "mode": mode}
    if mode == "stack_uniform":
        if stack_w is None:
            return {"loss": zero, "kl": float("nan"), "mode": mode}
        pack = alpha_stack_kl(
            stack_w, uniform_alpha(mods), mods, lambda_kl=lambda_kl
        )
        return {**pack, "mode": mode}
    if mode == "stack_erank":
        if hiddens is None or stack_w is None:
            return {"loss": zero, "kl": float("nan"), "mode": mode}
        pack = erank_balance_loss(
            hiddens,
            mods,
            stack_w=stack_w,
            alpha=alpha,
            lambda_bal=lambda_bal,
        )
        return {
            "loss": pack["loss"],
            "kl": pack["l_bal"],
            "erank": pack["erank"],
            "mode": mode,
        }
    raise ValueError(f"unknown weight mode: {mode}")


def uses_stack_fusion(mode: str) -> bool:
    return mode != "mean_ce"


def freezes_stack_psi(mode: str) -> bool:
    return mode == "stack_fixed"


def stack_temperature_for(mode: str, *, default: float = 1.0, temp: float = 0.5) -> float:
    return float(temp) if mode == "stack_temp" else float(default)
