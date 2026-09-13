#!/usr/bin/env python3
"""Gradient-cosine × per-modality LR on the online-window stream.

Characterize modality grads with cosine similarity, then optionally fold
alignment into the L2 actuator:

  equal        — dense equal LR
  soft         — continuous α→LR_m
  soft_gradcos — soft × ½(1+cos(g_m, g_shared))

Also logs temporal cos(g_t, g_{t-1}) per modality.

  PYTHONPATH=. python3 scripts/run_agod_gradcos_lr.py
  PYTHONPATH=. python3 scripts/run_agod_gradcos_lr.py --datasets msrvtt
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split
from torchvision import transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from agod.adapter import soft_lr_dispersion
from agod.lr_controller import (
    EMARouter,
    aligned_cos_sim,
    common_dim_grad_signature,
    cos_sim,
    modality_grad_cosine,
    schedule_modality_lr,
    softmax_scores,
)
from agod.shift import decompose_hybrid, decompose_mmd, residual_concept

OUT = ROOT / "results" / "agod_gradcos_lr"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_gradcos_lr")

SEED = 2026
EMA, TAU, BETA, HOLD = 0.40, 0.30, 0.10, 0.35
MSRVTT_MODS = ["video", "text", "audio"]
N_REF, T_WIN = 400, 6
N_CUR_MSRVTT, N_CUR_AMAZON = 240, 120
STEPS_M, BATCH, LR0, FUSE = 35, 64, 3e-3, 128
SCHEDULERS = ("equal", "soft", "soft_gradcos")
LAMBDA_ALIGN = 0.75


def load_msrvtt():
    d = ROOT / "data" / "msrvtt" / "packed"
    feats = {
        "video": np.load(d / "video_feat.npy").astype(np.float32),
        "text": np.load(d / "text_feat.npy").astype(np.float32),
        "audio": np.load(d / "audio_feat.npy").astype(np.float32),
    }
    y = np.load(d / "labelsmsr.npy").astype(np.int64)
    if y.ndim > 1:
        y = y.reshape(len(y), -1)[:, 0]
    if len(np.unique(y)) > 2:
        mode = int(np.bincount(y.astype(int)).argmax())
        y = (y == mode).astype(np.int64)
    return feats, y


def make_stream_y(y, n_cur: int):
    rng = np.random.default_rng(SEED)
    idx0, idx1 = np.where(y == 0)[0].copy(), np.where(y == 1)[0].copy()
    rng.shuffle(idx0)
    rng.shuffle(idx1)
    n0, n1 = min(N_REF // 2, len(idx0)), min(N_REF - N_REF // 2, len(idx1))
    ref = np.concatenate([idx0[:n0], idx1[:n1]])
    rng.shuffle(ref)
    rest0, rest1 = idx0[n0:], idx1[n1:]

    def take(pool, k, t):
        if len(pool) == 0 or k <= 0:
            return np.array([], dtype=int)
        s = (t * max(k, 1)) % len(pool)
        return np.concatenate([pool[s:], pool[:s]])[:k]

    windows = []
    for t in range(T_WIN):
        p1 = 0.25 + 0.50 * (t / max(T_WIN - 1, 1))
        k1 = int(n_cur * p1)
        cur = np.concatenate([take(rest0, n_cur - k1, t), take(rest1, k1, t)])
        rng.shuffle(cur)
        windows.append({"t": t, "p1": float(p1), "idx": cur})
    return {"ref_idx": ref, "windows": windows}


def blocks(feats, y, idx, mods):
    return {m: feats[m][idx] for m in mods}, y[idx].astype(float)


def domain_msg(b0, b1, y0, y1, mods, *, seed):
    auc, vimp, po = {}, {}, {}
    for i, m in enumerate(mods):
        X = np.vstack([b0[m], b1[m]])
        W = np.array([0] * len(b0[m]) + [1] * len(b1[m]))
        if len(np.unique(W)) < 2 or len(X) < 40:
            a, v = 0.5, 0.0
        else:
            Xtr, Xte, Wtr, Wte = train_test_split(
                X, W, test_size=0.3, random_state=seed + i, stratify=W
            )
            clf = RandomForestClassifier(
                n_estimators=40, max_depth=5, min_samples_leaf=3,
                random_state=seed + i, n_jobs=1,
            )
            clf.fit(Xtr, Wtr)
            try:
                a = float(roc_auc_score(Wte, clf.predict_proba(Xte)[:, 1]))
            except Exception:
                a = 0.5
            v = float(clf.feature_importances_.mean())
        auc[m], vimp[m] = a, v
        po[m] = residual_concept(b0[m], y0, b1[m], y1, seed=seed + 31 + i)
    return SimpleNamespace(auc=auc, vimp=vimp, po=po)


class Fusion(nn.Module):
    def __init__(self, dims, mods):
        super().__init__()
        self.mods = list(mods)
        self.projs = nn.ModuleDict(
            {
                m: nn.Sequential(nn.Linear(dims[m], FUSE), nn.ReLU(), nn.Dropout(0.1))
                for m in mods
            }
        )
        self.head = nn.Linear(FUSE, 2)

    def forward(self, batch):
        hs = [self.projs[m](batch[m]) for m in self.mods]
        return self.head(torch.stack(hs, 0).mean(0))

    def param_groups(self):
        return {m: list(self.projs[m].parameters()) for m in self.mods} | {
            "shared": list(self.head.parameters())
        }


def build_optim(model, mods):
    groups = model.param_groups()
    return torch.optim.Adam(
        [{"params": groups[m], "lr": LR0, "name": m} for m in list(mods) + ["shared"]]
    )


def set_lrs(opt, lr_mult):
    for g in opt.param_groups:
        g["lr"] = LR0 * float(lr_mult.get(g["name"], lr_mult.get("shared", 1.0)))


def to_batch(feats, y, idx, mods, device):
    batch = {
        m: torch.tensor(feats[m][idx], dtype=torch.float32, device=device) for m in mods
    }
    yy = torch.tensor(y[idx], dtype=torch.long, device=device)
    return batch, yy


@torch.no_grad()
def eval_acc(model, feats, y, idx, mods, device):
    model.eval()
    if len(idx) < 8:
        return {"acc": float("nan"), "n": len(idx)}
    logits, ys = [], []
    for s in range(0, len(idx), BATCH):
        sl = idx[s : s + BATCH]
        batch, yy = to_batch(feats, y, sl, mods, device)
        logits.append(model(batch).cpu().numpy())
        ys.append(yy.cpu().numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    return {"acc": float(accuracy_score(Y, L.argmax(1))), "n": int(len(Y))}


def split_hold(idx, seed):
    rng = np.random.default_rng(seed)
    idx = np.asarray(idx)
    perm = rng.permutation(len(idx))
    n_h = max(16, min(len(idx) // 2, int(len(idx) * HOLD)))
    return idx[perm[n_h:]], idx[perm[:n_h]]


def probe_grads(model, feats, y, idx, mods, device):
    """One forward/backward to snapshot modality / shared gradients."""
    model.train()
    model.zero_grad(set_to_none=True)
    rng = np.random.default_rng(SEED + 99)
    sel = rng.choice(idx, size=min(BATCH, len(idx)), replace=len(idx) < BATCH)
    batch, yy = to_batch(feats, y, sel, mods, device)
    loss = nn.CrossEntropyLoss()(model(batch), yy)
    loss.backward()
    groups = model.param_groups()
    # common-dim signatures so video/text/audio towers remain comparable
    grads = {m: common_dim_grad_signature(groups[m]) for m in mods}
    shared = common_dim_grad_signature(groups["shared"])
    return grads, shared, float(loss.item())


def train_window(model, opt, feats, y, idx, mods, device, *, lr_mult, steps):
    model.train()
    set_lrs(opt, lr_mult)
    crit = nn.CrossEntropyLoss()
    losses = []
    rng = np.random.default_rng(SEED + int(np.asarray(idx).sum()) % 100000)
    t0 = time.perf_counter()
    for _ in range(steps):
        sel = rng.choice(idx, size=min(BATCH, len(idx)), replace=len(idx) < BATCH)
        batch, yy = to_batch(feats, y, sel, mods, device)
        loss = crit(model(batch), yy)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        losses.append(float(loss.item()))
    return float(np.mean(losses)), (time.perf_counter() - t0) * 1000.0


def select_b5_raw(msg, b0, b1, y0, y1, mods, *, seed):
    mmd_d = decompose_mmd(b0, b1, y0, y1, mods, seed=seed)
    hy = decompose_hybrid(msg, mmd_d, mods)
    return softmax_scores(hy["score"], mods, TAU), hy


def run_msrvtt(device, sched: str):
    feats, y = load_msrvtt()
    mods = list(MSRVTT_MODS)
    stream = make_stream_y(y, N_CUR_MSRVTT)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = Fusion({m: feats[m].shape[1] for m in mods}, mods).to(device)
    opt = build_optim(model, mods)
    router = EMARouter(mods, ema=EMA)
    ref_idx = stream["ref_idx"].copy()
    ref_b, y_ref = blocks(feats, y, ref_idx, mods)
    warm, _ = split_hold(ref_idx, SEED)
    train_window(
        model, opt, feats, y, warm, mods, device,
        lr_mult={**{m: 1.0 for m in mods}, "shared": 1.0},
        steps=STEPS_M,
    )
    prev_grads = None
    traj = []
    for w in stream["windows"]:
        adapt, hold = split_hold(w["idx"], SEED + 13 * w["t"])
        cur_b, y_cur = blocks(feats, y, adapt, mods)
        msg = domain_msg(ref_b, cur_b, y_ref, y_cur, mods, seed=SEED + 10 * w["t"])
        raw, _ = select_b5_raw(
            msg, ref_b, cur_b, y_ref, y_cur, mods, seed=SEED + 10 * w["t"]
        )
        alpha = router.update(raw)

        grads, shared, _ = probe_grads(model, feats, y, adapt, mods, device)
        ginfo = modality_grad_cosine(grads, mods, shared=shared)
        temporal = {}
        if prev_grads is not None:
            for m in mods:
                temporal[m] = aligned_cos_sim(grads[m], prev_grads[m])
        prev_grads = {m: grads[m].copy() for m in mods}

        lr_mult = schedule_modality_lr(
            sched,
            alpha,
            mods,
            t=w["t"],
            t_max=T_WIN,
            beta=BETA,
            align_gain=ginfo["align_gain"],
            lambda_align=LAMBDA_ALIGN,
        )
        pre = eval_acc(model, feats, y, hold, mods, device)
        loss, wall = train_window(
            model, opt, feats, y, adapt, mods, device, lr_mult=lr_mult, steps=STEPS_M
        )
        post = eval_acc(model, feats, y, hold, mods, device)
        dacc = post["acc"] - pre["acc"]
        disp = soft_lr_dispersion(lr_mult, mods)
        keep = N_REF // 2
        ref_idx = np.concatenate([ref_idx[-keep:], adapt[: min(keep, len(adapt))]])
        ref_b, y_ref = blocks(feats, y, ref_idx, mods)
        row = {
            "t": w["t"],
            "n_adapt": int(len(adapt)),
            "alpha": {m: float(alpha[m]) for m in mods},
            "lr_mult": {k: float(v) for k, v in lr_mult.items()},
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": dacc,
            "train_loss": loss,
            "wall_ms": wall,
            "mean_pair_cos": ginfo["mean_pair_cos"],
            "frac_conflict": ginfo["frac_conflict"],
            "cos_to_shared": {m: float(ginfo["cos_to_shared"][m]) for m in mods},
            "align_gain": {m: float(ginfo["align_gain"][m]) for m in mods},
            "temporal_cos": temporal,
            **disp,
        }
        traj.append(row)
        print(
            f"[msrvtt/{sched}] t={w['t']} pair_cos={ginfo['mean_pair_cos']:+.2f} "
            f"conflict={ginfo['frac_conflict']:.2f} lr_ratio={disp['lr_ratio']:.2f} "
            f"acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj, mods


def _amazon_stream(samples, n_cur, base):
    by_cat: dict[str, list] = {}
    for s in samples:
        by_cat.setdefault(s["category"], []).append(s)
    cats = sorted(by_cat, key=lambda c: -len(by_cat[c]))
    ref = by_cat[cats[0]][: max(base.N_REF, min(250, len(by_cat[cats[0]])))]
    others = [c for c in cats[1:] if len(by_cat[c]) >= max(24, n_cur // 2)] or cats[1:4]
    windows = []
    for t in range(base.T_STEPS):
        cat = others[t % len(others)]
        pool = list(by_cat[cat])
        rng = np.random.default_rng(SEED + 7 * t)
        rng.shuffle(pool)
        n = min(int(n_cur), len(pool))
        mix_n = min(int((0.35 * (1 - t / max(base.T_STEPS - 1, 1))) * n), len(ref), n)
        cur = list(ref[:mix_n]) + pool[: n - mix_n]
        rng.shuffle(cur)
        windows.append({"t": t, "category": cat, "rows": cur})
    return {"ref": ref, "windows": windows}


@torch.no_grad()
def _amazon_eval(model, rows, device, image_tf, base):
    from torch.utils.data import DataLoader

    model.eval()
    if len(rows) < 8:
        return {"acc": float("nan"), "n": len(rows)}
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=base.BATCH, shuffle=False)
    logits, ys = [], []
    for imgs, txts, y in loader:
        logits.append(model(imgs.to(device), txts, return_mods=False).cpu().numpy())
        ys.append(y.numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    return {"acc": float(accuracy_score(Y, L.argmax(1))), "n": int(len(Y))}


def _amazon_probe_grads(model, rows, device, image_tf, base, mods):
    from torch.utils.data import DataLoader

    model.train()
    model.zero_grad(set_to_none=True)
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=base.BATCH, shuffle=True)
    imgs, txts, y = next(iter(loader))
    imgs, y = imgs.to(device), y.to(device)
    loss = nn.CrossEntropyLoss()(model(imgs, txts, return_mods=False), y)
    loss.backward()
    groups = model.param_groups()
    grads = {m: common_dim_grad_signature(groups[m]) for m in mods}
    shared = common_dim_grad_signature(groups["shared"])
    return grads, shared, float(loss.item())


def _amazon_train(model, opt, rows, device, image_tf, lr_mult, base, steps):
    from torch.utils.data import DataLoader

    model.train()
    if hasattr(model, "image_encoder"):
        model.image_encoder.eval()
    base.set_lrs(opt, lr_mult)
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=base.BATCH, shuffle=True)
    crit = nn.CrossEntropyLoss()
    losses, it = [], iter(loader)
    t0 = time.perf_counter()
    for _ in range(steps):
        try:
            imgs, txts, y = next(it)
        except StopIteration:
            it = iter(loader)
            imgs, txts, y = next(it)
        imgs, y = imgs.to(device), y.to(device)
        loss = crit(model(imgs, txts, return_mods=False), y)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        losses.append(float(loss.item()))
    return float(np.mean(losses)), (time.perf_counter() - t0) * 1000.0


def run_amazon(device, sched: str, *, samples=None, image_tf=None):
    import run_agod_amazon_modality_lr as base

    mods = list(base.MODS)
    if samples is None:
        shards = sorted(base.SHARD_DIR.glob("*.tar.gz"))
        if not shards:
            raise SystemExit(f"No Amazon shards in {base.SHARD_DIR}")
        print(f"loading {len(shards)} amazon shards…", flush=True)
        samples = base.load_shards(shards)
        print(
            f"samples={len(samples)} cats={Counter(s['category'] for s in samples).most_common(5)}",
            flush=True,
        )
    if image_tf is None:
        image_tf = transforms.Compose(
            [
                transforms.Resize((224, 224)),
                transforms.ToTensor(),
                transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
            ]
        )
    stream = _amazon_stream(samples, N_CUR_AMAZON, base)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = base.AmazonAGOD().to(device)
    opt = base.build_optim(model)
    router = EMARouter(mods, ema=EMA)
    ref = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref, device, image_tf)
    steps = getattr(base, "STEPS_PER_WIN", 8)
    prev_grads = None
    traj = []
    for w in stream["windows"]:
        rng = np.random.default_rng(SEED + 11 * w["t"])
        rows = list(w["rows"])
        rng.shuffle(rows)
        n_h = max(12, min(len(rows) // 2, int(len(rows) * 0.40)))
        hold, adapt = rows[:n_h], rows[n_h:]
        if len(adapt) < 8 or len(hold) < 8:
            continue
        cur_b, y_cur = base.extract_blocks(model, adapt, device, image_tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        raw = msg.alpha["B3"]
        alpha = router.update(raw)

        grads, shared, _ = _amazon_probe_grads(
            model, adapt, device, image_tf, base, mods
        )
        ginfo = modality_grad_cosine(grads, mods, shared=shared)
        temporal = {}
        if prev_grads is not None:
            for m in mods:
                temporal[m] = aligned_cos_sim(grads[m], prev_grads[m])
        prev_grads = {m: grads[m].copy() for m in mods}

        lr_mult = schedule_modality_lr(
            sched,
            alpha,
            mods,
            t=w["t"],
            t_max=base.T_STEPS,
            beta=BETA,
            align_gain=ginfo["align_gain"],
            lambda_align=LAMBDA_ALIGN,
        )
        pre = _amazon_eval(model, hold, device, image_tf, base)
        loss, wall = _amazon_train(
            model, opt, adapt, device, image_tf, lr_mult, base, steps
        )
        post = _amazon_eval(model, hold, device, image_tf, base)
        dacc = post["acc"] - pre["acc"]
        disp = soft_lr_dispersion(lr_mult, mods)
        ref_b = {
            m: np.vstack(
                [ref_b[m][-(base.N_REF // 2) :], cur_b[m][: base.N_REF // 2]]
            )
            for m in mods
        }
        y_ref = np.concatenate(
            [y_ref[-(base.N_REF // 2) :], y_cur[: base.N_REF // 2]]
        )
        row = {
            "t": w["t"],
            "category": w["category"],
            "n_adapt": int(len(adapt)),
            "alpha": {m: float(alpha[m]) for m in mods},
            "lr_mult": {k: float(v) for k, v in lr_mult.items()},
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": dacc,
            "train_loss": loss,
            "wall_ms": wall,
            "mean_pair_cos": ginfo["mean_pair_cos"],
            "frac_conflict": ginfo["frac_conflict"],
            "cos_to_shared": {m: float(ginfo["cos_to_shared"][m]) for m in mods},
            "align_gain": {m: float(ginfo["align_gain"][m]) for m in mods},
            "temporal_cos": temporal,
            **disp,
        }
        traj.append(row)
        print(
            f"[amazon/{sched}] t={w['t']} cat={str(w['category'])[:18]} "
            f"pair_cos={ginfo['mean_pair_cos']:+.2f} conflict={ginfo['frac_conflict']:.2f} "
            f"lr_ratio={disp['lr_ratio']:.2f} "
            f"acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj, mods, samples, image_tf


def summarize(traj, mods, *, dataset, sched):
    lifts = np.array([r["acc_lift"] for r in traj], float)
    pair = np.array([r["mean_pair_cos"] for r in traj], float)
    conf = np.array([r["frac_conflict"] for r in traj], float)
    temp = []
    for r in traj:
        if r.get("temporal_cos"):
            temp.extend(list(r["temporal_cos"].values()))
    return {
        "dataset": dataset,
        "scheduler": sched,
        "n_windows": len(traj),
        "mean_acc_lift": float(np.nanmean(lifts)),
        "frac_pos": float(np.mean(lifts > 0)),
        "mean_lr_ratio": float(np.nanmean([r["lr_ratio"] for r in traj])),
        "mean_pair_cos": float(np.nanmean(pair)),
        "mean_frac_conflict": float(np.nanmean(conf)),
        "mean_temporal_cos": float(np.nanmean(temp)) if temp else float("nan"),
        "mean_acc_post": float(np.nanmean([r["acc_post"] for r in traj])),
    }


def plot_board(cells, path: Path):
    datasets = sorted({c["dataset"] for c in cells})
    scheds = [s for s in SCHEDULERS if any(c["scheduler"] == s for c in cells)]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), facecolor="#f7f5f1")
    fig.suptitle(
        "Gradient-cosine × modality LR (Amazon / MSR-VTT)",
        fontsize=13,
        fontweight="bold",
    )
    x = np.arange(len(scheds))
    w = 0.35
    metrics = [
        ("mean_acc_lift", "mean Acc lift", True),
        ("mean_pair_cos", "mean pairwise grad cos", False),
        ("mean_frac_conflict", "frac conflict (cos<0)", False),
    ]
    for ax, (key, title, hline) in zip(axes, metrics):
        for i, ds in enumerate(datasets):
            vals = []
            for s in scheds:
                row = next(
                    (c for c in cells if c["dataset"] == ds and c["scheduler"] == s),
                    None,
                )
                vals.append(row[key] if row else np.nan)
            ax.bar(x + (i - 0.5) * w, vals, w, label=ds)
        if hline:
            ax.axhline(0, color="#999", ls="--", lw=0.8)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(scheds, rotation=20, ha="right")
        ax.legend(frameon=False, fontsize=8)
    fig.text(
        0.5,
        0.02,
        "align_gain=½(1+cos(g_m,g_shared)) · soft_gradcos = soft × [(1-λ)+λ·align] · FWD always on",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(cells, path: Path):
    rows = [
        f"| {c['dataset']} | {c['scheduler']} | {c['mean_acc_lift']:+.3f} | "
        f"{c['mean_pair_cos']:+.3f} | {c['mean_frac_conflict']:.2f} | "
        f"{c['mean_temporal_cos']:+.3f} | {c['mean_lr_ratio']:.2f} |"
        for c in sorted(cells, key=lambda r: (r["dataset"], r["scheduler"]))
    ]
    notes = []
    for ds in sorted({c["dataset"] for c in cells}):
        eq = next(c for c in cells if c["dataset"] == ds and c["scheduler"] == "equal")
        others = [c for c in cells if c["dataset"] == ds and c["scheduler"] != "equal"]
        best = max(others, key=lambda c: c["mean_acc_lift"])
        notes.append(
            f"- **{ds}**: best `{best['scheduler']}` lift {best['mean_acc_lift']:+.3f} "
            f"(vs equal {eq['mean_acc_lift']:+.3f}, "
            f"Δ={best['mean_acc_lift'] - eq['mean_acc_lift']:+.3f}); "
            f"pair_cos={best['mean_pair_cos']:+.3f}, conflict={best['mean_frac_conflict']:.2f}"
        )
    md = f"""# Gradient-cosine × modality LR

## Point

在 online-window 流上，用 **梯度余弦相似度** 刻画模态更新几何，并可选地折进 L2 步长：

| 量 | 含义 |
|---|---|
| `cos(g_m, g_shared)` | 模态投影梯度与 head 的对齐 |
| `pair_cos` | 模态两两梯度夹角（<0 ⇒ 冲突） |
| `temporal_cos` | 跨窗 `cos(g_t, g_{{t-1}})`（步长是否抖） |
| `align_gain` | `½(1+cos_to_shared)` ∈ [0,1] |
| `soft_gradcos` | `LR_m = soft(α_m) · [(1-λ)+λ·align_gain_m]` |

不必 dense equal：α 负责 **谁该大步**；grad-cos 负责 **谁的更新方向不打架**。

## Smoke

| Dataset | scheduler | Acc lift | pair_cos | conflict | temporal_cos | lr_ratio |
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

### Readout
{chr(10).join(notes)}

### Other opportunities
- **Conflict damp**: if `cos(g_m,g_shared)<0`, force damp schedule on that mod
- **Temporal gate**: low `temporal_cos` → raise β (flatten) until grads stabilize
- **Budget × align**: redistribute fixed excess-LR by `α_m · align_gain_m`
- **Sample-grain grads**: cosine on per-example grads → replay weights (different actuator)

```bash
PYTHONPATH=. python3 scripts/run_agod_gradcos_lr.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(cells, path: Path):
    lines = [
        f"{c['dataset']} & {c['scheduler'].replace('_', '\\\\_')} & "
        f"{c['mean_acc_lift']:+.3f} & {c['mean_pair_cos']:+.3f} & "
        f"{c['mean_frac_conflict']:.2f} & {c['mean_temporal_cos']:+.3f} & "
        f"{c['mean_lr_ratio']:.2f} \\\\"
        for c in sorted(cells, key=lambda r: (r["dataset"], r["scheduler"]))
    ]
    tex = (
        "% Gradient-cosine × modality LR\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Per-modality LR with gradient cosine alignment on the online-window stream.}\n"
        "\\label{tab:agod-gradcos-lr}\n"
        "\\begin{tabular}{llccccc}\\toprule\n"
        "Dataset & scheduler & Acc lift & pair\\_cos & conflict & temporal\\_cos & lr\\_ratio \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["msrvtt", "amazon"])
    ap.add_argument("--schedulers", nargs="+", default=list(SCHEDULERS))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    cells, traj_dump = [], {}
    amazon_cache = {"samples": None, "image_tf": None}

    if "msrvtt" in args.datasets:
        for sched in args.schedulers:
            traj, mods = run_msrvtt(device, sched)
            cells.append(summarize(traj, mods, dataset="msrvtt", sched=sched))
            traj_dump[f"msrvtt:{sched}"] = traj

    if "amazon" in args.datasets:
        for sched in args.schedulers:
            traj, mods, samples, image_tf = run_amazon(
                device,
                sched,
                samples=amazon_cache["samples"],
                image_tf=amazon_cache["image_tf"],
            )
            amazon_cache["samples"] = samples
            amazon_cache["image_tf"] = image_tf
            cells.append(summarize(traj, mods, dataset="amazon", sched=sched))
            traj_dump[f"amazon:{sched}"] = traj

    payload = {
        "agod_version": "0.1.0",
        "focus": "gradient cosine × modality LR on online-window stream",
        "schedulers": list(args.schedulers),
        "design": {
            "equal": "dense equal LR",
            "soft": "continuous α→LR_m",
            "soft_gradcos": "soft × [(1-λ)+λ·½(1+cos(g_m,g_shared))]",
            "pair_cos": "mean pairwise modality-grad cosine",
            "frac_conflict": "fraction of pairs with cos<0",
            "temporal_cos": "cos(g_t, g_{t-1}) per modality",
        },
        "cells": cells,
        "trajectory": traj_dump,
    }
    (OUT / "agod_gradcos_lr.json").write_text(json.dumps(payload, indent=2))
    plot_board(cells, OUT / "AGOD_GradCos_LR_Board.png")
    write_docs(cells, OUT / "README.md")
    write_latex(cells, OUT / "AGOD_gradcos_lr_tables_only.tex")
    write_docs(cells, DOCS / "AGOD_gradcos_lr.md")
    write_latex(cells, DOCS / "AGOD_gradcos_lr_tables_only.tex")
    shutil.copy2(OUT / "AGOD_GradCos_LR_Board.png", ART / "AGOD_GradCos_LR_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    print("\n=== grad-cos board ===", flush=True)
    for c in cells:
        print(
            f"  {c['dataset']:8s} {c['scheduler']:12s} "
            f"lift={c['mean_acc_lift']:+.3f} pair_cos={c['mean_pair_cos']:+.3f} "
            f"conflict={c['mean_frac_conflict']:.2f}",
            flush=True,
        )
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
