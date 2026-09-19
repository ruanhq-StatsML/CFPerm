"""Disentangled Style Transfer with FSDS Attribution-Driven Routing.

Prototype:

  1. ContentTower + StyleTower
  2. Four losses: style / align / contrast / fluency + MINE
  3. FSDS (MMD / CMean / PO) + LOGO groups → alpha routing
  4. Attribution-driven LR modulation (soft adapter)

Serving table (Y is never a feature):
  X = source token ids as float columns
  Y = target style id (s_tgt)
  groups: content=topic words, style_src=source style markers, noise=pad

Real-deploy swap points:
  ContentTower.backbone → hfl/chinese-roberta-wwm-ext
  fsds_route            → RF-domain + MMD-LOCO + PO-risk (this file)
  decoder               → GPT-2 / Qwen prefix-tuning

Run:
    python Python/src/style_transfer_fsds.py
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional

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
    logo_batch,
)
from streaming_po_risk import large_deviation

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


# ============================================================
# 1. Synthetic data
# ============================================================

def make_data(n, cfg: Config, seed=0, drift=False, concept=False):
    """Synthetic style-transfer rows.

    content : topic id (0..n_topics-1)
    s_src   : source style
    s_tgt   : target style (≠ s_src unless concept=True redraws it)
    src     : topic words (4) + source style markers (4) + noise

    drift=True   : noise ~ Unif(0, vocab) instead of 400..vocab  (covariate)
    concept=True : s_tgt independent of src                     (P(Y|X) hop)

    Y for FSDS is s_tgt. It is not written into src.
    """
    rng = np.random.RandomState(seed)
    topics = rng.randint(0, cfg.n_topics, size=n)
    s_src = rng.randint(0, cfg.num_styles, size=n)
    offset = rng.randint(1, cfg.num_styles, size=n)
    s_tgt = (s_src + offset) % cfg.num_styles
    if concept:
        s_tgt = rng.randint(0, cfg.num_styles, size=n)

    topic_words = rng.randint(10, 400, size=(cfg.n_topics, N_TOPIC_TOKS))
    style_tokens = (
        500
        + np.arange(cfg.num_styles)[:, None] * 30
        + rng.randint(0, 30, size=(cfg.num_styles, cfg.n_style_markers))
    )

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

    Real deploy: self.backbone = AutoModel.from_pretrained(
        "hfl/chinese-roberta-wwm-ext")
    """

    def __init__(self, cfg: Config):
        super().__init__()
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
        x = self.encoder(x)
        x = x.mean(dim=1)
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
        h_f = self.fusion(torch.cat([h_c, h_s], dim=-1))
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
            + self.cfg.w_mi * l_mi
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


def attribution_to_alpha(cov: float, con: float, cfg: Config) -> float:
    """(cov, con) excess → alpha in (0, 1).

    cov high → P(X) walk → content / stem
    con high → P(Y|X) hop → style / top
    """
    s = cfg.cov_weight * float(cov) + cfg.con_weight * float(con)
    return float(1.0 / (1.0 + math.exp(-s / max(cfg.alpha_tau, 1e-6))))


def lr_scale(alpha: float, cfg: Config) -> float:
    return cfg.lr_lo + (cfg.lr_hi - cfg.lr_lo) * float(alpha)


def _quiet_split(X, Y, seed: int):
    rng = np.random.RandomState(int(seed))
    idx = rng.permutation(len(X))
    h = max(len(X) // 2, 1)
    return X[idx[:h]], Y[idx[:h]], X[idx[h:]], Y[idx[h:]]


def _excess(stream: Optional[float], quiet: Optional[float]) -> float:
    if stream is None:
        return 0.0
    q = 0.0 if quiet is None else float(quiet)
    return max(float(stream) / max(q, 1e-8) - 1.0, 0.0)


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
    style_src       → StyleTower   (concept lives in P(Y|X))
    Shares are localization, not Shapley / CATE. T is the batch label.
    """
    scales = {
        "content": cfg.lr_lo,
        "style": 1.0,
        "fusion": 1.0,
        "decoder": 1.0,
        "style_cls": 1.0,
        "content_cls": cfg.lr_lo,
        "mine": 0.5,
    }
    mapping = {
        "content": ("content", "content_cls"),
        "style_src": ("style", "style_cls"),
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
        scales["style_cls"] = max(scales["style_cls"], cfg.lr_hi)
    else:
        scales["fusion"] = min(scales["fusion"], cfg.lr_lo)
        scales["decoder"] = min(scales["decoder"], cfg.lr_lo)
    return scales


def fsds_route(
    exist: Mapping,
    new: Mapping,
    cfg: Config,
    *,
    seed: int = 2026,
    with_logo: bool = True,
) -> dict:
    """RF-domain FSDS + MMD-LOCO + PO-risk on the serving table.

    Replaces the toy neural FSDSAttribution head.
    """
    X_e, Y_e = serving_xy(exist)
    X_n, Y_n = serving_xy(new)
    names = feature_names(cfg)
    if X_e.shape[1] != len(names):
        raise ValueError(f"X cols ({X_e.shape[1]}) != names ({len(names)})")
    assert_not_outcome(names)

    metrics = three_metrics(X_e, Y_e, X_n, Y_n, seed=seed, with_po=True)
    Xq1, Yq1, Xq2, Yq2 = _quiet_split(X_e, Y_e, seed)
    quiet = three_metrics(Xq1, Yq1, Xq2, Yq2, seed=seed, with_po=True)

    ranked = fsds_rank_columns(X_e, X_n, Y_e, Y_n, names, seed=seed)
    groups = style_groups(cfg)
    logo = None
    if with_logo:
        logo = logo_batch(
            X_e,
            Y_e,
            X_n,
            Y_n,
            groups,
            seed=seed,
            po_base=quiet.get("po"),
            mse_base=quiet.get("mse"),
            mmd_base=quiet.get("mmd"),
        )

    cov = _excess(metrics["mmd"], quiet["mmd"])
    con = _excess(metrics["po"], quiet["po"])
    if metrics["po"] is None:
        con = _excess(abs(metrics["cmean_y"]), abs(quiet["cmean_y"]) + 1e-8)
    alpha = attribution_to_alpha(cov, con, cfg)
    plan = None if logo is None else logo["plan"]
    return {
        "metrics": metrics,
        "quiet": quiet,
        "cov": float(cov),
        "con": float(con),
        "alpha": float(alpha),
        "lr_scale": float(lr_scale(alpha, cfg)),
        "fsds_top": ranked[:8],
        "fsds_all": ranked,
        "feature_names": names,
        "groups": {k: [int(groups[k].start), int(groups[k].stop)] for k in groups},
        "logo": logo,
        "plan": plan,
        "y_in_x": False,
        "mmd_broken": bool(large_deviation(metrics["mmd"], quiet["mmd"])),
        "po_broken": bool(
            metrics["po"] is not None
            and quiet["po"] is not None
            and large_deviation(metrics["po"], quiet["po"])
        ),
    }


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
    print(
        f"[FSDS] mmd={m['mmd']:.4f}  po={po_s}  "
        f"cmean_x={m['cmean_x']:.4f}  cmean_y={m['cmean_y']:+.4f}  "
        f"cov={route['cov']:.3f}  con={route['con']:.3f}  "
        f"alpha={route['alpha']:.3f}  lr_scale={route['lr_scale']:.3f}"
    )
    top = ", ".join(f"{r['feature']}:{r['score']:.2f}" for r in route["fsds_top"][:5])
    print(f"[FSDS] top = {top}")
    plan = route.get("plan")
    if not plan:
        return
    towers = {g: v["tower"] for g, v in plan.get("towers", {}).items()}
    print(
        f"[LOGO] action={plan.get('global_action')}  update={plan.get('update')}  "
        f"dominant={plan.get('dominant')}  fusion={plan.get('fusion')}"
    )
    print(f"[LOGO] towers = {towers}")


def train(cfg: Config):
    set_seed(cfg.seed)
    device = cfg.device

    train_data = make_data(cfg.n_train, cfg, seed=1, drift=False)
    val_data = make_data(cfg.n_val, cfg, seed=2, drift=False)
    exist_data = make_data(cfg.n_attr, cfg, seed=3, drift=False)
    new_data = make_data(cfg.n_attr, cfg, seed=4, drift=True)

    train_loader = DataLoader(PairDataset(train_data), batch_size=cfg.batch_size, shuffle=True)
    val_loader = DataLoader(PairDataset(val_data), batch_size=cfg.batch_size, shuffle=False)

    model = DisentangledStyleTransfer(cfg).to(device)
    loss_fn = MultiLoss(cfg)

    print("[FSDS] covariate walk vs D_ref (noise support) ...")
    route = fsds_route(exist_data, new_data, cfg, seed=cfg.seed, with_logo=True)
    _print_route(route)

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
                    f"contrast={parts['contrast']:.4f}  "
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
        batch = next(iter(val_loader))
        batch = {k: v.to(device) for k, v in batch.items()}
        out = model(batch["src"], batch["s_tgt"])
        mi_final = out["mi"].item()
        leak = float("nan")
        try:
            from sklearn.linear_model import LogisticRegression

            h_c_np = out["h_c"].cpu().numpy()
            s_np = batch["s_tgt"].cpu().numpy()
            if len(set(s_np.tolist())) > 1:
                clf = LogisticRegression(max_iter=200).fit(h_c_np, s_np)
                leak = clf.score(h_c_np, s_np)
        except Exception:
            pass
        print(f"  MI lower bound    = {mi_final:+.4f}")
        print(f"  h_c -> style leak = {leak:.4f}  (lower is better)")

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
