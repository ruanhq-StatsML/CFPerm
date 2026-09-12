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
    ):
        super().__init__()
        self.mods = list(mods)
        self.n_class = int(n_class)
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

    def stack_weights(self) -> torch.Tensor:
        return F.softmax(self.stack_logits, dim=0)

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

    def forward(self, batch: Mapping[str, torch.Tensor], *, return_parts: bool = False):
        h = self.encode(batch)
        logits_m = self.modality_logits(h)
        w = self.stack_weights()
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
