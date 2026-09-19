"""Disentangled Style Transfer with FSDS Attribution-Driven Routing.

Prototype:

  1. ContentTower (style tokens masked) + StyleTower + prefix adapter
  2. Four losses: style / align / contrast / fluency + MINE
  3. Group-first FSDS (MMD / CMean / PO) + LOGO → alpha routing
  4. Soft-adapter LR: MMD-loud stem / PO-loud prefix+top

Serving table (Y is never a feature):
  X = source token ids as float columns
  Y = target style id (s_tgt)
  groups: content=topic words, style_src=source style markers, noise=pad

User-style entry:
    styleTransferFSDS(df, groups, ref_batch_size, batch_size)
    last column of df is Y

Real-deploy swap points:
  ContentTower.backbone → hfl/chinese-roberta-wwm-ext
  PrefixAdapter         → GPT-2 / Qwen prefix-tuning
  fsds_route            → RF-domain + MMD-LOCO + PO-risk (this file)

Run:
    python Python/src/style_transfer_fsds.py
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from graph_fsds_localize import assert_not_outcome, fsds_rank_columns, three_metrics
from logo_modality import (
    FUSION_FULL,
    FUSION_HEAD,
    TOWER_FREEZE,
    TOWER_FULL,
    TOWER_INFER,
    TOWER_STEM,
    TOWER_TOP,
    as_groups,
    logo_batch,
    plan_next_batch,
)
from streaming_po_risk import ACTION_FREEZE, ACTION_KEEP, ACTION_XSHIFT, large_deviation

QUIET_FLOOR_MMD = 0.015
QUIET_FLOOR_PO = 0.01
QUIET_FLOOR_CMEAN_Y = 0.10

N_TOPIC_TOKS = 4


# ============================================================
# 0. Config
# ============================================================

@dataclass
class Config:
    vocab_size: int = 2000
    max_len: int = 32
    n_topics: int = 20
    n_style_markers: int = 4
    n_prefix: int = 4

    embed_dim: int = 96
    content_dim: int = 128
    num_heads: int = 4
    num_layers: int = 2

    num_styles: int = 6
    style_dim: int = 32
    style_hidden: int = 128

    decoder_hidden: int = 192
    mine_hidden: int = 128

    w_style: float = 1.0
    w_align: float = 1.0
    w_contrast: float = 0.5
    w_fluency: float = 1.0
    w_mi: float = 0.05
    contrast_tau: float = 0.1

    alpha_tau: float = 0.5
    cov_weight: float = 1.0
    con_weight: float = 1.0
    lr_base: float = 1e-3
    lr_lo: float = 0.3
    lr_hi: float = 2.0

    batch_size: int = 32
    num_steps: int = 80
    eval_every: int = 20
    n_train: int = 800
    n_val: int = 200
    n_attr: int = 64
    device: str = "cuda" if torch.cuda.is_available() else "cpu"
    seed: int = 42


def set_seed(s: int) -> None:
    np.random.seed(s)
    torch.manual_seed(s)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(s)


def _as_np(t) -> np.ndarray:
    if isinstance(t, torch.Tensor):
        return t.detach().cpu().numpy()
    return np.asarray(t)


def _group_idx(groups: Mapping, name: str) -> np.ndarray:
    v = groups[name]
    if isinstance(v, slice):
        return np.arange(int(v.start), int(v.stop))
    return np.asarray(v, dtype=int)


# ============================================================
# 1. Synthetic data
# ============================================================

def make_lexicon(cfg: Config, seed: int = 0):
    """Frozen token codebook. Shared across D_ref / new so only drift flags move X."""
    rng = np.random.RandomState(seed)
    topic_words = rng.randint(10, 400, size=(cfg.n_topics, N_TOPIC_TOKS))
    style_tokens = (
        500
        + np.arange(cfg.num_styles)[:, None] * 30
        + rng.randint(0, 30, size=(cfg.num_styles, cfg.n_style_markers))
    )
    return topic_words, style_tokens


def make_data(n, cfg: Config, seed=0, drift=False, concept=False, lexicon=None):
    """Synthetic style-transfer rows.

    drift=True   : noise ~ Unif(0, vocab) instead of 400..vocab  (covariate)
    concept=True : s_tgt = (s_src + 1 + K/2) mod K              (P(Y|X) hop)

    Y for FSDS is s_tgt. It is not written into src.
    """
    rng = np.random.RandomState(seed)
    topics = rng.randint(0, cfg.n_topics, size=n)
    s_src = rng.randint(0, cfg.num_styles, size=n)
    s_tgt = (s_src + 1) % cfg.num_styles
    if concept:
        s_tgt = (s_src + 1 + cfg.num_styles // 2) % cfg.num_styles

    topic_words, style_tokens = lexicon if lexicon is not None else make_lexicon(cfg, seed=0)

    src = np.zeros((n, cfg.max_len), dtype=np.int64)
    tgt = np.zeros((n, cfg.n_style_markers), dtype=np.int64)
    for i in range(n):
        parts = list(topic_words[topics[i]]) + list(style_tokens[s_src[i]])
        while len(parts) < cfg.max_len:
            if drift:
                parts.append(rng.randint(0, cfg.vocab_size))
            else:
                parts.append(rng.randint(400, cfg.vocab_size))
        src[i] = np.array(parts[: cfg.max_len])
        tgt[i] = style_tokens[s_tgt[i]]

    return {
        "src": torch.tensor(src, dtype=torch.long),
        "tgt": torch.tensor(tgt, dtype=torch.long),
        "s_tgt": torch.tensor(s_tgt, dtype=torch.long),
        "topic": torch.tensor(topics, dtype=torch.long),
    }


class PairDataset(Dataset):
    def __init__(self, d):
        self.d = d

    def __len__(self):
        return self.d["src"].size(0)

    def __getitem__(self, i):
        return {k: v[i] for k, v in self.d.items()}


# ============================================================
# 2. Model
# ============================================================

class ContentTower(nn.Module):
    """Token sequence → content vector.

    Style-marker positions are padding-masked so the content pool cannot
    read source style. Real deploy: AutoModel.from_pretrained(
        "hfl/chinese-roberta-wwm-ext") with the same span mask.
    """

    def __init__(self, cfg: Config):
        super().__init__()
        self.style_lo = N_TOPIC_TOKS
        self.style_hi = N_TOPIC_TOKS + int(cfg.n_style_markers)
        self.tok = nn.Embedding(cfg.vocab_size, cfg.embed_dim)
        self.pos = nn.Embedding(cfg.max_len, cfg.embed_dim)
        layer = nn.TransformerEncoderLayer(
            cfg.embed_dim,
            cfg.num_heads,
            cfg.embed_dim * 2,
            batch_first=True,
            dropout=0.1,
            activation="gelu",
        )
        self.encoder = nn.TransformerEncoder(layer, cfg.num_layers)
        self.proj = nn.Sequential(
            nn.Linear(cfg.embed_dim, cfg.embed_dim),
            nn.GELU(),
            nn.LayerNorm(cfg.embed_dim),
            nn.Linear(cfg.embed_dim, cfg.content_dim),
        )

    def forward(self, ids):
        B, L = ids.shape
        pos = torch.arange(L, device=ids.device).unsqueeze(0).expand(B, L)
        x = self.tok(ids) + self.pos(pos)
        pad = torch.zeros(B, L, dtype=torch.bool, device=ids.device)
        lo, hi = self.style_lo, min(self.style_hi, L)
        if hi > lo:
            pad[:, lo:hi] = True
            x = x.masked_fill(pad.unsqueeze(-1), 0.0)
        x = self.encoder(x, src_key_padding_mask=pad)
        keep = (~pad).float().unsqueeze(-1)
        x = (x * keep).sum(dim=1) / keep.sum(dim=1).clamp(min=1.0)
        return self.proj(x)


class StyleTower(nn.Module):
    """Discrete style label → narrow style vector. Bottleneck so it cannot hold content."""

    def __init__(self, cfg: Config):
        super().__init__()
        self.emb = nn.Embedding(cfg.num_styles, cfg.style_hidden)
        self.net = nn.Sequential(
            nn.Linear(cfg.style_hidden, cfg.style_hidden),
            nn.GELU(),
            nn.LayerNorm(cfg.style_hidden),
            nn.Linear(cfg.style_hidden, cfg.style_dim),
        )

    def forward(self, s):
        return self.net(self.emb(s))


class PrefixAdapter(nn.Module):
    """Learnable prefix on the fused hidden.

    Real deploy: GPT-2 / Qwen prefix-tuning (virtual tokens in the decoder).
    Here the prefix is a soft adapter on h_f — same slot, toy decoder.
    """

    def __init__(self, cfg: Config):
        super().__init__()
        self.prefix = nn.Parameter(torch.zeros(int(cfg.n_prefix), cfg.decoder_hidden))
        nn.init.normal_(self.prefix, std=0.02)
        self.proj = nn.Linear(cfg.decoder_hidden, cfg.decoder_hidden)

    def forward(self, h_f):
        p = self.proj(self.prefix.mean(dim=0))
        return h_f + p.unsqueeze(0)


class MINE(nn.Module):
    """Mutual-information lower bound (Belghazi et al. 2018)."""

    def __init__(self, c_dim, s_dim, hidden):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(c_dim + s_dim, hidden),
            nn.GELU(),
            nn.Linear(hidden, hidden),
            nn.GELU(),
            nn.Linear(hidden, 1),
        )

    def forward(self, h_c, h_s):
        pos = self.net(torch.cat([h_c, h_s], dim=-1))
        idx = torch.randperm(h_s.size(0), device=h_s.device)
        neg = self.net(torch.cat([h_c, h_s[idx]], dim=-1))
        mi = pos.mean() - (torch.logsumexp(neg, dim=0) - math.log(neg.size(0)))
        return mi


class DisentangledStyleTransfer(nn.Module):
    def __init__(self, cfg: Config):
        super().__init__()
        self.cfg = cfg
        self.content = ContentTower(cfg)
        self.style = StyleTower(cfg)
        self.fusion = nn.Sequential(
            nn.Linear(cfg.content_dim + cfg.style_dim, cfg.decoder_hidden),
            nn.GELU(),
            nn.LayerNorm(cfg.decoder_hidden),
        )
        self.prefix = PrefixAdapter(cfg)
        self.decoder = nn.Sequential(
            nn.Linear(cfg.decoder_hidden, cfg.decoder_hidden),
            nn.GELU(),
            nn.LayerNorm(cfg.decoder_hidden),
            nn.Linear(cfg.decoder_hidden, cfg.n_style_markers * cfg.vocab_size),
        )
        self.style_cls = nn.Linear(cfg.decoder_hidden, cfg.num_styles)
        self.content_cls = nn.Linear(cfg.content_dim, cfg.n_topics)
        self.mine = MINE(cfg.content_dim, cfg.style_dim, cfg.mine_hidden)

    def forward(self, src_ids, tgt_style_ids):
        h_c = self.content(src_ids)
        h_s = self.style(tgt_style_ids)
        h_f = self.prefix(self.fusion(torch.cat([h_c, h_s], dim=-1)))
        logits = self.decoder(h_f).view(-1, self.cfg.n_style_markers, self.cfg.vocab_size)
        return {
            "h_c": h_c,
            "h_s": h_s,
            "h_f": h_f,
            "logits": logits,
            "style_logits": self.style_cls(h_f),
            "content_logits": self.content_cls(h_c),
            "mi": self.mine(h_c, h_s),
        }


# ============================================================
# 3. Losses
# ============================================================

class MultiLoss(nn.Module):
    def __init__(self, cfg: Config):
        super().__init__()
        self.cfg = cfg

    def style_loss(self, style_logits, s_tgt):
        return F.cross_entropy(style_logits, s_tgt)

    def align_loss(self, content_logits, topic):
        return F.cross_entropy(content_logits, topic)

    def contrast_loss(self, h_s, s_tgt):
        z = F.normalize(h_s, dim=-1)
        sim = z @ z.t() / self.cfg.contrast_tau
        mask = (s_tgt[:, None] == s_tgt[None, :]).float()
        mask.fill_diagonal_(0)
        log_prob = F.log_softmax(sim, dim=-1)
        pos = (log_prob * mask).sum(-1) / mask.sum(-1).clamp(min=1)
        return -pos.mean()

    def fluency_loss(self, logits, tgt_tokens):
        B, K, V = logits.shape
        return F.cross_entropy(logits.reshape(B * K, V), tgt_tokens.reshape(B * K))

    def forward(self, out, batch):
        l_s = self.style_loss(out["style_logits"], batch["s_tgt"])
        l_a = self.align_loss(out["content_logits"], batch["topic"])
        l_c = self.contrast_loss(out["h_s"], batch["s_tgt"])
        l_f = self.fluency_loss(out["logits"], batch["tgt"])
        l_mi = out["mi"]
        total = (
            self.cfg.w_style * l_s
            + self.cfg.w_align * l_a
            + self.cfg.w_contrast * l_c
            + self.cfg.w_fluency * l_f
            + self.cfg.w_mi * torch.clamp(l_mi, min=-2.0, max=5.0)
        )
        parts = {
            "style": l_s.item(),
            "align": l_a.item(),
            "contrast": l_c.item(),
            "fluency": l_f.item(),
            "mi": l_mi.item(),
            "total": total.item(),
        }
        return total, parts


# ============================================================
# 4. FSDS + LOGO → alpha / plan / LR
# ============================================================

def style_groups(cfg: Config) -> dict:
    """Partition of X. Y (s_tgt) is not a group and not a column."""
    style_hi = N_TOPIC_TOKS + int(cfg.n_style_markers)
    return {
        "content": slice(0, N_TOPIC_TOKS),
        "style_src": slice(N_TOPIC_TOKS, style_hi),
        "noise": slice(style_hi, int(cfg.max_len)),
    }


def feature_names(cfg: Config) -> list[str]:
    names = [f"content_{j}" for j in range(N_TOPIC_TOKS)]
    names += [f"style_src_{j}" for j in range(cfg.n_style_markers)]
    names += [f"noise_{j}" for j in range(N_TOPIC_TOKS + cfg.n_style_markers, cfg.max_len)]
    return names


def serving_xy(data: Mapping) -> tuple[np.ndarray, np.ndarray]:
    """X = src tokens as float. Y = s_tgt. Y is not concatenated onto X."""
    X = _as_np(data["src"]).astype(float)
    if X.ndim == 1:
        X = X.reshape(1, -1)
    Y = _as_np(data["s_tgt"]).astype(float).ravel()
    return X, Y


def pack_style_df(data: Mapping) -> np.ndarray:
    """Last column Y. Same contract as onlinePermOOB_with_LLM."""
    X, Y = serving_xy(data)
    return np.column_stack([X, Y.reshape(-1, 1)])


def attribution_to_alpha(cov: float, con: float, cfg: Config) -> float:
    s = cfg.cov_weight * float(cov) + cfg.con_weight * float(con)
    return float(1.0 / (1.0 + math.exp(-s / max(cfg.alpha_tau, 1e-6))))


def lr_scale(alpha: float, cfg: Config) -> float:
    return cfg.lr_lo + (cfg.lr_hi - cfg.lr_lo) * float(alpha)


def _quiet_pair(X, Y, seed: int, n: int = 40):
    """Own-ref pair with n≥40 so PO-risk is defined. Bootstrap, same P(X,Y)."""
    rng = np.random.RandomState(int(seed))
    k = max(int(n), 40)
    i1 = rng.choice(len(X), size=k, replace=True)
    i2 = rng.choice(len(X), size=k, replace=True)
    return X[i1], Y[i1], X[i2], Y[i2]


def _quiet_split(X, Y, seed: int):
    return _quiet_pair(X, Y, seed)


def _excess(stream: Optional[float], quiet: Optional[float], floor: float) -> float:
    if stream is None:
        return 0.0
    q = max(0.0 if quiet is None else float(quiet), float(floor))
    return max(float(stream) / q - 1.0, 0.0)


def group_three_metrics(X_e, Y_e, X_n, Y_n, groups, seed: int) -> dict:
    """FSDS three-metrics on each modality block. Grain first, then columns."""
    groups = as_groups(groups)
    out = {}
    for name, idx in groups.items():
        out[name] = three_metrics(
            X_e[:, idx], Y_e, X_n[:, idx], Y_n, seed=seed, with_po=True
        )
    return out


WIDE_GROUP = 8


def sketch_groups(X, groups, names):
    """Wide pad blocks → 4-stat sketch so LOGO is not dim-dominated."""
    groups = as_groups(groups)
    X = np.asarray(X, dtype=float)
    blocks, new_g, new_names, col = [], {}, [], 0
    for g, idx in groups.items():
        block = X[:, idx]
        gn = [names[int(i)] for i in idx]
        if block.shape[1] > WIDE_GROUP:
            stats = np.column_stack(
                [block.mean(1), block.std(1), block.min(1), block.max(1)]
            )
            blocks.append(stats)
            new_g[g] = np.arange(col, col + 4)
            new_names += [f"{g}_mean", f"{g}_std", f"{g}_min", f"{g}_max"]
            col += 4
        else:
            blocks.append(block)
            new_g[g] = np.arange(col, col + block.shape[1])
            new_names += gn
            col += block.shape[1]
    return np.hstack(blocks), new_g, new_names


def pick_loud_group(g_new: Mapping, logo: Optional[Mapping] = None) -> tuple[Optional[dict], list]:
    """Grain first: a 2× MMD peak is covariate; a 2× PO peak is concept.

    Localization, not a unique decomp. T is the batch label.
    """
    names = list(g_new)
    mmds = {g: max(float(g_new[g]["mmd"]), 0.0) for g in names}
    pos = {g: max(float(g_new[g]["po"] or 0.0), 0.0) for g in names}
    rows = [
        {
            "group": g,
            "cov": float(mmds[g]),
            "con": float(pos[g]),
            "kind": "mmd" if mmds[g] >= pos[g] else "po",
            "score": float(max(mmds[g], pos[g])),
        }
        for g in names
    ]
    g_m = max(mmds, key=mmds.get)
    rest_m = [mmds[g] for g in names if g != g_m]
    if mmds[g_m] >= QUIET_FLOOR_MMD and mmds[g_m] >= 2.0 * max(rest_m + [1e-12]):
        hit = {"group": g_m, "kind": "mmd", "cov": mmds[g_m], "con": pos[g_m], "score": mmds[g_m]}
        rows.sort(key=lambda r: r["score"], reverse=True)
        return hit, rows
    g_p = max(pos, key=pos.get)
    rest_p = [pos[g] for g in names if g != g_p]
    if pos[g_p] >= QUIET_FLOOR_PO and pos[g_p] >= 2.0 * max(rest_p + [1e-12]):
        hit = {"group": g_p, "kind": "po", "cov": mmds[g_p], "con": pos[g_p], "score": pos[g_p]}
        rows.sort(key=lambda r: r["score"], reverse=True)
        return hit, rows
    if logo is not None:
        pi_mmd = logo.get("pi_mmd") or {}
        pi_po = logo.get("pi_po") or {}
        if pi_mmd and max(pi_mmd.values()) >= max(list(pi_po.values()) + [0.0]) and max(pi_mmd.values()) > 0:
            g = max(pi_mmd, key=pi_mmd.get)
            hit = {"group": g, "kind": "mmd", "cov": mmds.get(g, 0.0), "con": pos.get(g, 0.0), "score": pi_mmd[g]}
            rows.sort(key=lambda r: r["score"], reverse=True)
            return hit, rows
        if pi_po and max(pi_po.values()) > 0:
            g = max(pi_po, key=pi_po.get)
            hit = {"group": g, "kind": "po", "cov": mmds.get(g, 0.0), "con": pos.get(g, 0.0), "score": pi_po[g]}
            rows.sort(key=lambda r: r["score"], reverse=True)
            return hit, rows
    rows.sort(key=lambda r: r["score"], reverse=True)
    if not rows or rows[0]["score"] <= 0:
        return None, rows
    return rows[0], rows


def rank_inside_group(X_e, X_n, Y_e, Y_n, names, groups, loud_name, seed: int):
    idx = _group_idx(groups, loud_name)
    sub = [names[int(j)] for j in idx]
    assert_not_outcome(sub)
    return fsds_rank_columns(X_e[:, idx], X_n[:, idx], Y_e, Y_n, sub, seed=seed)


def overlay_fsds_plan(
    logo: Mapping,
    mmd_broken: bool,
    po_broken: bool,
    loud: Optional[Mapping] = None,
) -> dict:
    """Overlay only when the global MMD or PO gate is actually broken.

    A 2× peak among groups on a quiet window is sampling noise, not a hop.
    """
    plan = dict(logo["plan"])
    raw = str(plan.get("global_action"))
    if not mmd_broken and not po_broken:
        return plan
    if loud and float(loud.get("score") or 0) > 0:
        g = str(loud["group"])
        if loud["kind"] == "mmd" and mmd_broken:
            plan = plan_next_batch(ACTION_XSHIFT, logo["ratios"])
            if g in plan["towers"]:
                plan["towers"][g] = {
                    "tower": TOWER_STEM,
                    "i_star": 0,
                    "reason": "FSDS grain MMD-loud: stem-adapt",
                }
            plan["overlay"] = "fsds_group_mmd"
            plan["global_action_raw"] = raw
            return plan
        if loud["kind"] == "po" and po_broken:
            plan = plan_next_batch(ACTION_FREEZE, logo["ratios"])
            if g in plan["towers"]:
                plan["towers"][g] = {
                    "tower": TOWER_TOP,
                    "i_star": "top",
                    "reason": "FSDS grain PO-loud: train top + prefix",
                }
            plan["overlay"] = "fsds_group_po"
            plan["global_action_raw"] = raw
            return plan
    if raw == ACTION_KEEP and mmd_broken and not po_broken:
        plan = plan_next_batch(ACTION_XSHIFT, logo["ratios"])
        plan["overlay"] = "fsds_mmd_keep"
        plan["global_action_raw"] = raw
    return plan


def _action_scale(action: str, cfg: Config) -> float:
    if action == TOWER_STEM:
        return cfg.lr_hi
    if action == TOWER_TOP:
        return cfg.lr_hi
    if action == TOWER_FULL:
        return 1.0
    if action in (TOWER_FREEZE, TOWER_INFER):
        return cfg.lr_lo
    return 1.0


def module_scales_from_plan(plan: Mapping, cfg: Config) -> dict[str, float]:
    """Map LOGO tower actions onto the style-transfer modules.

    content / noise → ContentTower (covariate lives in the token stem)
    style_src       → StyleTower + prefix (concept lives in P(Y|X))
    Shares are localization, not Shapley / CATE. T is the batch label.
    """
    scales = {
        "content": cfg.lr_lo,
        "style": 1.0,
        "fusion": 1.0,
        "prefix": 1.0,
        "decoder": 1.0,
        "style_cls": 1.0,
        "content_cls": cfg.lr_lo,
        "mine": 0.5,
    }
    mapping = {
        "content": ("content", "content_cls"),
        "style_src": ("style", "style_cls", "prefix"),
        "noise": ("content",),
    }
    for g, mods in mapping.items():
        tower = plan.get("towers", {}).get(g, {})
        action = tower.get("tower", TOWER_FULL)
        sc = _action_scale(str(action), cfg)
        for m in mods:
            scales[m] = max(float(scales[m]), float(sc))
    fusion = str(plan.get("fusion", FUSION_FULL))
    if fusion == FUSION_FULL:
        scales["fusion"] = max(scales["fusion"], 1.0)
        scales["decoder"] = max(scales["decoder"], 1.0)
    elif fusion == FUSION_HEAD:
        scales["fusion"] = cfg.lr_hi
        scales["decoder"] = cfg.lr_hi
        scales["prefix"] = max(scales["prefix"], cfg.lr_hi)
        scales["style_cls"] = max(scales["style_cls"], cfg.lr_hi)
    else:
        scales["fusion"] = min(scales["fusion"], cfg.lr_lo)
        scales["decoder"] = min(scales["decoder"], cfg.lr_lo)
        if str(plan.get("overlay", "")).endswith("_mmd"):
            scales["prefix"] = min(scales["prefix"], cfg.lr_lo)
    return scales


def fsds_route_xy(
    X_e,
    Y_e,
    X_n,
    Y_n,
    groups,
    names,
    cfg: Config,
    *,
    seed: int = 2026,
    with_logo: bool = True,
) -> dict:
    """RF-domain FSDS + MMD-LOCO + PO-risk. Y is not in X."""
    X_e = np.asarray(X_e, dtype=float)
    X_n = np.asarray(X_n, dtype=float)
    Y_e = np.asarray(Y_e, dtype=float).ravel()
    Y_n = np.asarray(Y_n, dtype=float).ravel()
    names = [str(n) for n in names]
    if X_e.shape[1] != len(names):
        raise ValueError(f"X cols ({X_e.shape[1]}) != names ({len(names)})")
    assert_not_outcome(names)
    raw_groups = as_groups(groups)
    raw_names = list(names)
    X_e, groups, names = sketch_groups(X_e, raw_groups, raw_names)
    X_n, _, _ = sketch_groups(X_n, raw_groups, raw_names)

    metrics = three_metrics(X_e, Y_e, X_n, Y_n, seed=seed, with_po=True)
    Xq1, Yq1, Xq2, Yq2 = _quiet_pair(X_e, Y_e, seed)
    quiet = three_metrics(Xq1, Yq1, Xq2, Yq2, seed=seed, with_po=True)
    g_new = group_three_metrics(X_e, Y_e, X_n, Y_n, groups, seed)

    po_base = quiet.get("po")
    if po_base is None:
        po_base = QUIET_FLOOR_PO
    else:
        po_base = max(float(po_base), QUIET_FLOOR_PO)
    mmd_base = max(float(quiet.get("mmd") or 0.0), QUIET_FLOOR_MMD)
    mse_base = quiet.get("mse")

    logo = None
    if with_logo:
        logo = logo_batch(
            X_e,
            Y_e,
            X_n,
            Y_n,
            groups,
            seed=seed,
            po_base=po_base,
            mse_base=mse_base,
            mmd_base=mmd_base,
        )
    loud, group_board = pick_loud_group(g_new, logo)
    mmd_broken = bool(large_deviation(metrics["mmd"], mmd_base))
    po_broken = bool(metrics["po"] is not None and large_deviation(metrics["po"], po_base))
    if not mmd_broken and not po_broken:
        loud = None
    if loud is not None:
        ranked = rank_inside_group(X_e, X_n, Y_e, Y_n, names, groups, loud["group"], seed)
    else:
        ranked = fsds_rank_columns(X_e, X_n, Y_e, Y_n, names, seed=seed)

    cov = _excess(metrics["mmd"], quiet["mmd"], QUIET_FLOOR_MMD)
    con = _excess(metrics["po"], quiet["po"], QUIET_FLOOR_PO)
    if metrics["po"] is None:
        con = _excess(abs(metrics["cmean_y"]), abs(quiet["cmean_y"]), QUIET_FLOOR_CMEAN_Y)
    if loud is not None:
        cov, con = float(loud["cov"]), float(loud["con"])
    alpha = attribution_to_alpha(cov, con, cfg)
    plan = None if logo is None else overlay_fsds_plan(logo, mmd_broken, po_broken, loud)
    return {
        "metrics": metrics,
        "quiet": quiet,
        "group_metrics": g_new,
        "loud": loud,
        "group_board": group_board,
        "cov": float(cov),
        "con": float(con),
        "alpha": float(alpha),
        "lr_scale": float(lr_scale(alpha, cfg)),
        "fsds_top": ranked[:8],
        "fsds_all": ranked,
        "feature_names": names,
        "groups": {k: _group_idx(groups, k).tolist() for k in groups},
        "logo": logo,
        "plan": plan,
        "y_in_x": False,
        "mmd_broken": mmd_broken,
        "po_broken": po_broken,
    }


def fsds_route(
    exist: Mapping,
    new: Mapping,
    cfg: Config,
    *,
    seed: int = 2026,
    with_logo: bool = True,
) -> dict:
    X_e, Y_e = serving_xy(exist)
    X_n, Y_n = serving_xy(new)
    return fsds_route_xy(
        X_e,
        Y_e,
        X_n,
        Y_n,
        style_groups(cfg),
        feature_names(cfg),
        cfg,
        seed=seed,
        with_logo=with_logo,
    )


def styleTransferFSDS(
    df,
    groups,
    ref_batch_size,
    batch_size,
    seed=2026,
    names: Optional[Sequence[str]] = None,
    cfg: Optional[Config] = None,
) -> list[dict]:
    """Last column is Y. Groups index X only. Each trail batch vs frozen D_ref."""
    cfg = cfg or Config()
    df = np.asarray(df, dtype=float)
    if df.ndim != 2 or df.shape[1] < 2:
        raise ValueError("df needs X columns and Y last")
    Y = df[:, -1]
    X = df[:, :-1]
    g = as_groups(groups)
    max_j = max(int(idx.max()) for idx in g.values())
    if max_j >= X.shape[1]:
        raise ValueError("group index hits Y or past X")
    names = [str(n) for n in names] if names is not None else [f"x{j}" for j in range(X.shape[1])]
    assert_not_outcome(names)
    n_ref = int(ref_batch_size)
    X_e, Y_e = X[:n_ref], Y[:n_ref]
    recs = []
    t = 0
    i = n_ref
    bs = int(batch_size)
    while i < len(df):
        X_n, Y_n = X[i : i + bs], Y[i : i + bs]
        if len(X_n) < 8:
            break
        route = fsds_route_xy(X_e, Y_e, X_n, Y_n, g, names, cfg, seed=int(seed) + t)
        plan = route.get("plan") or {}
        loud = route.get("loud")
        recs.append(
            {
                "t": t,
                "n": int(len(X_n)),
                "alpha": route["alpha"],
                "loud_group": None if loud is None else loud["group"],
                "loud_kind": None if loud is None else loud["kind"],
                "action": plan.get("global_action"),
                "update": plan.get("update"),
                "overlay": plan.get("overlay"),
                "towers": {k: v["tower"] for k, v in plan.get("towers", {}).items()},
                "mmd": route["metrics"]["mmd"],
                "po": route["metrics"]["po"],
                "y_in_x": False,
            }
        )
        t += 1
        i += bs
    return recs


def routing_effect(cfg: Config, n_rep: int = 6, n: int = 64) -> dict:
    """Hit rate: quiet stays quiet; covariate → noise/mmd; concept → style_src/po."""
    lex = make_lexicon(cfg, seed=0)
    scenes = [
        ("quiet", {}, None, None),
        ("covariate", {"drift": True}, "noise", "mmd"),
        ("concept", {"concept": True}, "style_src", "po"),
    ]
    out = {}
    for name, kw, want_g, want_k in scenes:
        ok = 0
        for r in range(int(n_rep)):
            exist = make_data(n, cfg, seed=20 + r, lexicon=lex)
            new = make_data(n, cfg, seed=200 + r, lexicon=lex, **kw)
            route = fsds_route(exist, new, cfg, seed=7 + r)
            loud = route.get("loud")
            overlay = (route.get("plan") or {}).get("overlay")
            if name == "quiet":
                hit = overlay is None and loud is None
            else:
                hit = (
                    loud is not None
                    and loud.get("group") == want_g
                    and loud.get("kind") == want_k
                )
            ok += int(hit)
        out[name] = {
            "hit": int(ok),
            "n": int(n_rep),
            "rate": float(ok) / max(int(n_rep), 1),
            "want": "quiet" if name == "quiet" else f"{want_g}/{want_k}",
        }
    return out


def effect_table_latex(effect: Mapping) -> str:
    lines = [
        r"\begin{tabular}{@{}lccc@{}}",
        r"\toprule",
        r"scene & want & hit & rate \\",
        r"\midrule",
    ]
    for name in ("quiet", "covariate", "concept"):
        row = effect.get(name) or {}
        lines.append(
            f"{name} & {row.get('want','')} & {row.get('hit',0)}/{row.get('n',0)} & {row.get('rate',0):.2f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def routing_board_latex(rows: Sequence[Mapping]) -> str:
    """Two-scene board: which grain rang, which tower to touch next."""
    lines = [
        r"\begin{tabular}{@{}llll@{}}",
        r"\toprule",
        r"scene & FSDS grain & next tower & note \\",
        r"\midrule",
    ]
    for r in rows:
        scene = str(r.get("scene", r.get("t", "")))
        g = r.get("loud_group")
        k = r.get("loud_kind")
        grain = "---" if not g else f"{g}/{k}"
        towers = r.get("towers") or {}
        nxt = ",".join(f"{k}:{v}" for k, v in towers.items() if v not in (TOWER_INFER, TOWER_FREEZE))
        if not nxt:
            nxt = r.get("update", "")
        note = str(r.get("overlay") or r.get("action") or "")
        lines.append(f"{scene} & {grain} & {nxt} & {note} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def build_optimizer_with_routing(
    model,
    cfg: Config,
    alpha: float,
    plan: Optional[Mapping] = None,
) -> torch.optim.Optimizer:
    """Soft adapter: alpha + LOGO plan → per-tower LR. No hard freeze."""
    if plan is None:
        scale_style = lr_scale(alpha, cfg)
        scale_content = 1.0 / max(scale_style, 0.5)
        scales = {
            "style": scale_style,
            "content": 0.3 * scale_content,
            "fusion": scale_style,
            "prefix": scale_style,
            "decoder": scale_style,
            "style_cls": scale_style,
            "content_cls": 0.3 * scale_content,
            "mine": 0.5,
        }
    else:
        scales = module_scales_from_plan(plan, cfg)
    groups = [
        {"params": model.style.parameters(), "lr": cfg.lr_base * scales["style"]},
        {"params": model.content.parameters(), "lr": cfg.lr_base * scales["content"]},
        {"params": model.fusion.parameters(), "lr": cfg.lr_base * scales["fusion"]},
        {"params": model.prefix.parameters(), "lr": cfg.lr_base * scales["prefix"]},
        {"params": model.decoder.parameters(), "lr": cfg.lr_base * scales["decoder"]},
        {"params": model.style_cls.parameters(), "lr": cfg.lr_base * scales["style_cls"]},
        {"params": model.content_cls.parameters(), "lr": cfg.lr_base * scales["content_cls"]},
        {"params": model.mine.parameters(), "lr": cfg.lr_base * scales["mine"]},
    ]
    return torch.optim.Adam(groups)


# ============================================================
# 5. Training
# ============================================================

@torch.no_grad()
def evaluate(model, loader, cfg: Config) -> dict:
    model.eval()
    total = {"style_acc": 0.0, "content_acc": 0.0, "n": 0}
    for batch in loader:
        batch = {k: v.to(cfg.device) for k, v in batch.items()}
        out = model(batch["src"], batch["s_tgt"])
        style_pred = out["style_logits"].argmax(-1)
        content_pred = out["content_logits"].argmax(-1)
        total["style_acc"] += (style_pred == batch["s_tgt"]).float().sum().item()
        total["content_acc"] += (content_pred == batch["topic"]).float().sum().item()
        total["n"] += batch["src"].size(0)
    for k in ("style_acc", "content_acc"):
        total[k] /= max(total["n"], 1)
    return total


def _print_route(route: Mapping) -> None:
    m = route["metrics"]
    po = m["po"]
    po_s = "na" if po is None else f"{po:.4f}"
    loud = route.get("loud")
    loud_s = "none" if not loud else f"{loud['group']}/{loud['kind']}={loud['score']:.2f}"
    print(
        f"[FSDS] mmd={m['mmd']:.4f}  po={po_s}  "
        f"cmean_x={m['cmean_x']:.4f}  cmean_y={m['cmean_y']:+.4f}  "
        f"cov={route['cov']:.3f}  con={route['con']:.3f}  "
        f"alpha={route['alpha']:.3f}  loud={loud_s}"
    )
    top = ", ".join(f"{r['feature']}:{r['score']:.2f}" for r in route["fsds_top"][:5])
    print(f"[FSDS] top = {top}")
    plan = route.get("plan")
    if not plan:
        return
    towers = {g: v["tower"] for g, v in plan.get("towers", {}).items()}
    extra = ""
    if plan.get("overlay"):
        extra = f"  overlay={plan.get('overlay')}  raw={plan.get('global_action_raw')}"
    print(
        f"[LOGO] action={plan.get('global_action')}  update={plan.get('update')}  "
        f"dominant={plan.get('dominant')}  fusion={plan.get('fusion')}{extra}"
    )
    print(f"[LOGO] towers = {towers}")


def _write_board(rows, cfg: Config, effect: Optional[Mapping] = None) -> Path:
    tex = routing_board_latex(rows)
    extra = ""
    if effect:
        extra = r"\vspace{1.2em}" + "\n" + r"\paragraph{Hit rate.}" + "\n" + effect_table_latex(effect)
    body = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs}",
            r"\begin{document}",
            r"\paragraph{Claim.}",
            r"FSDS ranks a modality grain first (content / style / noise),",
            r"then columns inside that grain. $Y$ is the target style, never a feature.",
            r"$T$ is the batch label. Shares are localization, not Shapley / CATE.",
            r"MMD-loud grain $\to$ stem-adapt the content tower.",
            r"PO-loud grain $\to$ train style top + prefix (GPT-2 slot).",
            r"A quiet window must not overlay freeze. Overlay only if global MMD or PO is broken.",
            r"\vspace{0.8em}",
            tex,
            extra,
            r"\end{document}",
        ]
    )
    out_dir = Path(__file__).resolve().parents[2] / "docs"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "style_transfer_fsds.tex"
    path.write_text(body, encoding="utf-8")
    res = Path(__file__).resolve().parents[2] / "results" / "style_transfer_fsds"
    res.mkdir(parents=True, exist_ok=True)
    (res / "board.tex").write_text(tex + ("\n" + extra if extra else ""), encoding="utf-8")
    if effect:
        import json

        (res / "effect.json").write_text(json.dumps(effect, indent=2), encoding="utf-8")
    return path


def train(cfg: Config):
    set_seed(cfg.seed)
    device = cfg.device
    lex = make_lexicon(cfg, seed=0)

    train_data = make_data(cfg.n_train, cfg, seed=1, lexicon=lex)
    val_data = make_data(cfg.n_val, cfg, seed=2, lexicon=lex)
    exist_data = make_data(cfg.n_attr, cfg, seed=3, lexicon=lex)
    cov_data = make_data(cfg.n_attr, cfg, seed=4, drift=True, lexicon=lex)
    con_data = make_data(cfg.n_attr, cfg, seed=5, concept=True, lexicon=lex)

    train_loader = DataLoader(PairDataset(train_data), batch_size=cfg.batch_size, shuffle=True)
    val_loader = DataLoader(PairDataset(val_data), batch_size=cfg.batch_size, shuffle=False)

    model = DisentangledStyleTransfer(cfg).to(device)
    loss_fn = MultiLoss(cfg)

    print("[FSDS] covariate walk vs D_ref (noise support) ...")
    route = fsds_route(exist_data, cov_data, cfg, seed=cfg.seed, with_logo=True)
    _print_route(route)
    print("[FSDS] concept remap of s_tgt vs D_ref ...")
    route_c = fsds_route(exist_data, con_data, cfg, seed=cfg.seed, with_logo=True)
    _print_route(route_c)
    quiet_data = make_data(cfg.n_attr, cfg, seed=6, lexicon=lex)
    print("[FSDS] quiet same-DGP vs D_ref ...")
    route_q = fsds_route(exist_data, quiet_data, cfg, seed=cfg.seed, with_logo=True)
    _print_route(route_q)

    print("[FSDS] routing hit rate (6 reps) ...")
    effect = routing_effect(cfg, n_rep=6, n=cfg.n_attr)
    for k, v in effect.items():
        print(f"  {k:10s}  want={v['want']:16s}  {v['hit']}/{v['n']}  ({v['rate']:.2f})")

    df = np.vstack([pack_style_df(exist_data), pack_style_df(cov_data), pack_style_df(con_data)])
    recs = styleTransferFSDS(
        df,
        style_groups(cfg),
        ref_batch_size=cfg.n_attr,
        batch_size=cfg.n_attr,
        seed=cfg.seed,
        names=feature_names(cfg),
        cfg=cfg,
    )
    print("[FSDS] styleTransferFSDS(df) batches:")
    for r in recs:
        print(
            f"  t={r['t']} loud={r['loud_group']}/{r['loud_kind']}  "
            f"update={r['update']}  overlay={r['overlay']}"
        )
    board_rows = [
        {
            "scene": "quiet",
            "loud_group": (route_q.get("loud") or {}).get("group"),
            "loud_kind": (route_q.get("loud") or {}).get("kind"),
            "towers": {k: v["tower"] for k, v in (route_q.get("plan") or {}).get("towers", {}).items()},
            "overlay": (route_q.get("plan") or {}).get("overlay"),
            "action": (route_q.get("plan") or {}).get("global_action"),
            "update": (route_q.get("plan") or {}).get("update"),
        },
        {
            "scene": "covariate noise",
            "loud_group": (route.get("loud") or {}).get("group"),
            "loud_kind": (route.get("loud") or {}).get("kind"),
            "towers": {k: v["tower"] for k, v in (route.get("plan") or {}).get("towers", {}).items()},
            "overlay": (route.get("plan") or {}).get("overlay"),
            "action": (route.get("plan") or {}).get("global_action"),
        },
        {
            "scene": "concept $s_{tgt}$",
            "loud_group": (route_c.get("loud") or {}).get("group"),
            "loud_kind": (route_c.get("loud") or {}).get("kind"),
            "towers": {k: v["tower"] for k, v in (route_c.get("plan") or {}).get("towers", {}).items()},
            "overlay": (route_c.get("plan") or {}).get("overlay"),
            "action": (route_c.get("plan") or {}).get("global_action"),
        },
    ]
    tex_path = _write_board(board_rows, cfg, effect=effect)
    print(f"[FSDS] board → {tex_path}")

    optimizer = build_optimizer_with_routing(model, cfg, route["alpha"], plan=route.get("plan"))

    print("[Train] start ...")
    step = 0
    while step < cfg.num_steps:
        for batch in train_loader:
            model.train()
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model(batch["src"], batch["s_tgt"])
            loss, parts = loss_fn(out, batch)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            if step % cfg.eval_every == 0:
                metrics = evaluate(model, val_loader, cfg)
                print(
                    f"step {step:4d}  "
                    f"loss={parts['total']:.4f}  "
                    f"style={parts['style']:.4f}  "
                    f"align={parts['align']:.4f}  "
                    f"fluency={parts['fluency']:.4f}  "
                    f"MI={parts['mi']:+.4f}  |  "
                    f"val style_acc={metrics['style_acc']:.3f}  "
                    f"content_acc={metrics['content_acc']:.3f}"
                )
            step += 1
            if step >= cfg.num_steps:
                break

    print("\n[Eval] final ...")
    final = evaluate(model, val_loader, cfg)
    print(f"  style_acc   = {final['style_acc']:.4f}")
    print(f"  content_acc = {final['content_acc']:.4f}")

    model.eval()
    with torch.no_grad():
        hs, ys = [], []
        for batch in val_loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            out = model(batch["src"], batch["s_tgt"])
            hs.append(out["h_c"].cpu())
            ys.append(batch["s_tgt"].cpu())
        h_c_np = torch.cat(hs).numpy()
        s_np = torch.cat(ys).numpy()
        mi_final = out["mi"].item()
        leak = float("nan")
        n = len(s_np)
        half = max(n // 2, 8)
        try:
            import warnings

            from sklearn.linear_model import LogisticRegression

            if len(set(s_np[:half].tolist())) > 1 and len(set(s_np[half:].tolist())) > 1:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    clf = LogisticRegression(max_iter=400).fit(h_c_np[:half], s_np[:half])
                leak = clf.score(h_c_np[half:], s_np[half:])
        except Exception:
            pass
        print(f"  MI lower bound    = {mi_final:+.4f}")
        print(f"  h_c -> style leak = {leak:.4f}  (held-out, lower is better)")

    return model, route, final


if __name__ == "__main__":
    cfg = Config()
    print("=" * 70)
    print("Disentangled Style Transfer + FSDS Attribution Routing")
    print("=" * 70)
    print(f"device     = {cfg.device}")
    print(f"num_styles = {cfg.num_styles}")
    print(f"vocab_size = {cfg.vocab_size}")
    print("=" * 70)
    train(cfg)
    print("\nDone.")
