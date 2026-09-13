"""Three-tower next-item model: user / sequence / item. Modest MLP, in-batch softmax."""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F


def l2norm(x: torch.Tensor, eps: float = 1e-6) -> torch.Tensor:
    return x / (x.norm(dim=-1, keepdim=True) + eps)


class MLP(nn.Module):
    def __init__(self, in_dim: int, hidden: int, out_dim: int, dropout: float = 0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, out_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class UserTower(nn.Module):
    def __init__(self, user_dim: int, hidden: int, out_dim: int):
        super().__init__()
        self.mlp = MLP(user_dim, hidden, out_dim)

    def forward(self, user_features: torch.Tensor) -> torch.Tensor:
        return l2norm(self.mlp(user_features))


class SeqTower(nn.Module):
    """Masked attention over history multimodal embeddings + action ids."""

    def __init__(self, emb_dim: int, hidden: int, out_dim: int, n_actions: int = 4):
        super().__init__()
        self.item_proj = nn.Linear(emb_dim, hidden)
        self.action_emb = nn.Embedding(n_actions, hidden)
        self.score = nn.Linear(hidden, 1)
        self.out = nn.Linear(hidden, out_dim)

    def forward(
        self,
        history_embs: torch.Tensor,
        history_actions: torch.Tensor,
        history_mask: torch.Tensor,
    ) -> torch.Tensor:
        acts = history_actions.clamp(min=0, max=self.action_emb.num_embeddings - 1)
        h = torch.tanh(self.item_proj(history_embs) + self.action_emb(acts))
        logits = self.score(h).squeeze(-1)
        logits = logits.masked_fill(history_mask <= 0, -1e9)
        attn = torch.softmax(logits, dim=-1)
        attn = attn * history_mask
        denom = attn.sum(dim=-1, keepdim=True).clamp_min(1e-6)
        pooled = (h * attn.unsqueeze(-1)).sum(dim=1) / denom
        return l2norm(self.out(pooled))


class ItemTower(nn.Module):
    def __init__(self, emb_dim: int, hidden: int, out_dim: int):
        super().__init__()
        self.mlp = MLP(emb_dim, hidden, out_dim)

    def forward(self, target_emb: torch.Tensor) -> torch.Tensor:
        return l2norm(self.mlp(target_emb))


class ThreeTowerModel(nn.Module):
    """Query = user ⊕ seq (sum then L2); item tower on target multimodal emb.

    In-batch InfoNCE for next-item plus a 3-way action head.
    Not a production GR ranker; a small trainable scaffold.
    """

    def __init__(
        self,
        user_dim: int,
        emb_dim: int,
        hidden_dim: int = 128,
        tower_dim: int = 64,
        n_actions: int = 3,
        temperature: float = 0.07,
    ):
        super().__init__()
        self.temperature = temperature
        self.user_tower = UserTower(user_dim, hidden_dim, tower_dim)
        self.seq_tower = SeqTower(emb_dim, hidden_dim, tower_dim)
        self.item_tower = ItemTower(emb_dim, hidden_dim, tower_dim)
        self.query_mix = nn.Linear(tower_dim * 2, tower_dim)
        self.action_head = nn.Linear(tower_dim * 2, n_actions)

    def encode(self, batch: Dict[str, torch.Tensor]) -> Tuple[torch.Tensor, torch.Tensor]:
        u = self.user_tower(batch["user_features"])
        s = self.seq_tower(batch["history_embs"], batch["history_actions"], batch["history_mask"])
        i = self.item_tower(batch["target_emb"])
        q = l2norm(self.query_mix(torch.cat([u, s], dim=-1)))
        return q, i

    def forward(self, batch: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        q, i = self.encode(batch)
        logits = (q @ i.T) / self.temperature
        labels = torch.arange(q.size(0), device=q.device)
        rec_loss = F.cross_entropy(logits, labels)
        act_logits = self.action_head(torch.cat([q, i], dim=-1))
        act_loss = F.cross_entropy(act_logits, batch["target_action"].clamp(min=0, max=act_logits.size(-1) - 1))
        loss = rec_loss + 0.5 * act_loss
        rec_acc = (logits.argmax(dim=-1) == labels).float().mean()
        return {
            "loss": loss,
            "rec_loss": rec_loss.detach(),
            "act_loss": act_loss.detach(),
            "rec_acc": rec_acc.detach(),
            "logits": logits,
        }


def three_tower_loss(model: ThreeTowerModel, batch: Dict[str, torch.Tensor]) -> torch.Tensor:
    return model(batch)["loss"]
