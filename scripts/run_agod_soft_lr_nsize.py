#!/usr/bin/env python3
"""Continuous soft-LR adapter + sample-size curves (Amazon / MSR-VTT).

Focus (not a gate zoo):
  1) How many adapt samples N are needed before Acc lift > 0 vs equal-LR?
  2) How to characterize / evaluate continuous next-stage step sizes (L2)?
  3) Soft continuous LR vs one sparse quantile gate (not dense equal).

Policies:
  equal     — dense equal LR (B1口径), flops_rel=1
  soft_lr   — continuous α→LR_m, all mods active (no hard gate)
  soft_sparse — soft_lr + quantile hard-gate (sparse adapt-BWD)

  PYTHONPATH=. python3 scripts/run_agod_soft_lr_nsize.py
  PYTHONPATH=. python3 scripts/run_agod_soft_lr_nsize.py --datasets msrvtt
"""
from __future__ import annotations

import argparse
import json
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

from agod.adapter import (
    AdapterConfig,
    LayeredAdapter,
    alpha_entropy,
    characterize_soft_lr_traj,
    flops_rel_proj,
    n_star_for_lift,
    soft_lr_dispersion,
)
from agod.lr_controller import EMARouter, alpha_to_lr, softmax_scores
from agod.shift import decompose_hybrid, decompose_mmd, residual_concept

OUT = ROOT / "results" / "agod_soft_lr_nsize"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_soft_lr_nsize")

SEED = 2026
EMA = 0.40
TAU = 0.30
BETA = 0.10
HOLD = 0.35
MSRVTT_MODS = ["video", "text", "audio"]
N_REF, T_WIN = 400, 6
STEPS_BASE, BATCH, LR0, FUSE = 35, 64, 3e-3, 128

# Sample-size grids (adapt pool size before holdout split)
MSRVTT_N_GRID = (80, 160, 240, 320, 480)
AMAZON_N_GRID = (40, 80, 120, 160, 200)
POLICIES = ("equal", "soft_lr", "soft_sparse")


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
    idx0 = np.where(y == 0)[0].copy()
    idx1 = np.where(y == 1)[0].copy()
    rng.shuffle(idx0)
    rng.shuffle(idx1)
    n0, n1 = N_REF // 2, N_REF - N_REF // 2
    n0, n1 = min(n0, len(idx0)), min(n1, len(idx1))
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
        windows.append({"t": t, "p1": float(p1), "idx": cur, "n_cur": int(n_cur)})
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
                n_estimators=40,
                max_depth=5,
                min_samples_leaf=3,
                random_state=seed + i,
                n_jobs=1,
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
        [
            {"params": groups[m], "lr": LR0, "name": m}
            for m in list(mods) + ["shared"]
        ]
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


def train_window(model, opt, feats, y, idx, mods, device, *, lr_mult, active, steps):
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
        groups = model.param_groups()
        for m in mods:
            if not active[m]:
                for p in groups[m]:
                    p.grad = None
        opt.step()
        losses.append(float(loss.item()))
    return float(np.mean(losses)), (time.perf_counter() - t0) * 1000.0


def select_b5_raw(msg, b0, b1, y0, y1, mods, *, seed):
    mmd_d = decompose_mmd(b0, b1, y0, y1, mods, seed=seed)
    hy = decompose_hybrid(msg, mmd_d, mods)
    return softmax_scores(hy["score"], mods, TAU), hy


def decide_policy(policy, alpha, mods, adapter: LayeredAdapter | None, *, t: int):
    if policy == "equal":
        lr = {**{m: 1.0 for m in mods}, "shared": 1.0}
        active = {m: True for m in mods}
        return lr, active, 0.0
    assert adapter is not None
    dec = adapter.step(alpha, seed=SEED, t=t)
    return dec.lr_mult, dec.active, dec.theta_used


def run_msrvtt_once(device, *, n_cur: int, policy: str, steps: int | None = None):
    steps = STEPS_BASE if steps is None else int(steps)
    feats, y = load_msrvtt()
    mods = list(MSRVTT_MODS)
    stream = make_stream_y(y, n_cur)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = Fusion({m: feats[m].shape[1] for m in mods}, mods).to(device)
    opt = build_optim(model, mods)
    router = EMARouter(mods, ema=EMA)
    gate = "none" if policy == "soft_lr" else ("quantile" if policy == "soft_sparse" else "none")
    adapter = None
    if policy != "equal":
        adapter = LayeredAdapter(
            mods,
            AdapterConfig(beta=BETA, q_drop=0.33),
            gate_policy=gate,
        )
    ref_idx = stream["ref_idx"].copy()
    ref_b, y_ref = blocks(feats, y, ref_idx, mods)
    warm, _ = split_hold(ref_idx, SEED)
    train_window(
        model,
        opt,
        feats,
        y,
        warm,
        mods,
        device,
        lr_mult={**{m: 1.0 for m in mods}, "shared": 1.0},
        active={m: True for m in mods},
        steps=steps,
    )
    traj = []
    for w in stream["windows"]:
        adapt, hold = split_hold(w["idx"], SEED + 13 * w["t"])
        cur_b, y_cur = blocks(feats, y, adapt, mods)
        msg = domain_msg(ref_b, cur_b, y_ref, y_cur, mods, seed=SEED + 10 * w["t"])
        raw, _ = select_b5_raw(
            msg, ref_b, cur_b, y_ref, y_cur, mods, seed=SEED + 10 * w["t"]
        )
        alpha = router.update(raw)
        lr_mult, active, th = decide_policy(policy, alpha, mods, adapter, t=w["t"])
        pre = eval_acc(model, feats, y, hold, mods, device)
        loss, wall = train_window(
            model,
            opt,
            feats,
            y,
            adapt,
            mods,
            device,
            lr_mult=lr_mult,
            active=active,
            steps=steps,
        )
        post = eval_acc(model, feats, y, hold, mods, device)
        dacc = post["acc"] - pre["acc"]
        fr = flops_rel_proj(active, mods)
        disp = soft_lr_dispersion(lr_mult, mods)
        keep = N_REF // 2
        ref_idx = np.concatenate([ref_idx[-keep:], adapt[: min(keep, len(adapt))]])
        ref_b, y_ref = blocks(feats, y, ref_idx, mods)
        row = {
            "t": w["t"],
            "n_cur": int(n_cur),
            "n_adapt": int(len(adapt)),
            "n_hold": int(len(hold)),
            "alpha": {m: float(alpha[m]) for m in mods},
            "lr_mult": {k: float(v) for k, v in lr_mult.items()},
            "active": active,
            "theta_used": th,
            "flops_rel": fr,
            "wall_ms": wall,
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": dacc,
            "train_loss": loss,
            "alpha_entropy": alpha_entropy(alpha, mods),
            **disp,
        }
        traj.append(row)
        print(
            f"[msrvtt/{policy}/N={n_cur}] t={w['t']} n_adapt={len(adapt)} "
            f"FLOPs={fr:.2f} lr_ratio={disp['lr_ratio']:.2f} "
            f"acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj, mods


def _amazon_make_stream(samples, n_cur: int, base):
    """Category-shift stream with controllable current-window size."""
    by_cat: dict[str, list] = {}
    for s in samples:
        by_cat.setdefault(s["category"], []).append(s)
    cats = sorted(by_cat, key=lambda c: -len(by_cat[c]))
    ref_cat = cats[0]
    ref = by_cat[ref_cat][: max(base.N_REF, min(250, len(by_cat[ref_cat])))]
    others = [c for c in cats[1:] if len(by_cat[c]) >= max(24, n_cur // 2)]
    if not others:
        others = cats[1:4] or cats[:1]
    windows = []
    for t in range(base.T_STEPS):
        cat = others[t % len(others)]
        pool = list(by_cat[cat])
        rng = np.random.default_rng(SEED + 7 * t)
        rng.shuffle(pool)
        n = min(int(n_cur), len(pool))
        # mix a shrinking fraction of ref to make gradual shift
        mix_n = int((0.35 * (1 - t / max(base.T_STEPS - 1, 1))) * n)
        mix_n = min(mix_n, len(ref), n)
        cur = list(ref[:mix_n]) + pool[: n - mix_n]
        rng.shuffle(cur)
        windows.append({"t": t, "category": cat, "rows": cur, "n_cur": n})
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


def _amazon_train(model, opt, rows, device, image_tf, lr_mult, active, mods, base, steps):
    from torch.utils.data import DataLoader

    model.train()
    if hasattr(model, "image_encoder"):
        model.image_encoder.eval()
    base.set_lrs(opt, lr_mult)
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=base.BATCH, shuffle=True)
    crit = nn.CrossEntropyLoss()
    losses, it = [], iter(loader)
    t0 = time.perf_counter()
    groups = model.param_groups()
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
        for m in mods:
            if not active.get(m, True):
                for p in groups[m]:
                    p.grad = None
        opt.step()
        losses.append(float(loss.item()))
    return float(np.mean(losses)), (time.perf_counter() - t0) * 1000.0


def run_amazon_once(device, *, n_cur: int, policy: str, samples=None, image_tf=None):
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
    stream = _amazon_make_stream(samples, n_cur, base)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = base.AmazonAGOD().to(device)
    opt = base.build_optim(model)
    router = EMARouter(mods, ema=EMA)
    gate = "none" if policy == "soft_lr" else ("quantile" if policy == "soft_sparse" else "none")
    adapter = None
    if policy != "equal":
        adapter = LayeredAdapter(
            mods,
            AdapterConfig(beta=BETA, q_drop=0.33, theta=0.40),
            gate_policy=gate,
        )
    ref = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref, device, image_tf)
    steps = getattr(base, "STEPS_PER_WIN", 8)
    traj = []
    for w in stream["windows"]:
        rng = np.random.default_rng(SEED + 11 * w["t"])
        rows = list(w["rows"])
        rng.shuffle(rows)
        n_h = max(12, min(len(rows) // 2, int(len(rows) * 0.40)))
        hold, adapt = rows[:n_h], rows[n_h:]
        if len(adapt) < 8 or len(hold) < 8:
            print(f"[amazon/{policy}/N={n_cur}] t={w['t']} skip tiny split", flush=True)
            continue
        cur_b, y_cur = base.extract_blocks(model, adapt, device, image_tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        raw = msg.alpha["B3"]
        alpha = router.update(raw)
        lr_mult, active, th = decide_policy(policy, alpha, mods, adapter, t=w["t"])
        pre = _amazon_eval(model, hold, device, image_tf, base)
        loss, wall = _amazon_train(
            model, opt, adapt, device, image_tf, lr_mult, active, mods, base, steps
        )
        post = _amazon_eval(model, hold, device, image_tf, base)
        dacc = post["acc"] - pre["acc"]
        fr = flops_rel_proj(active, mods)
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
            "n_cur": int(w["n_cur"]),
            "n_adapt": int(len(adapt)),
            "n_hold": int(len(hold)),
            "alpha": {m: float(alpha[m]) for m in mods},
            "lr_mult": {k: float(v) for k, v in lr_mult.items()},
            "active": active,
            "theta_used": th,
            "flops_rel": fr,
            "wall_ms": wall,
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": dacc,
            "train_loss": loss,
            "alpha_entropy": alpha_entropy(alpha, mods),
            **disp,
        }
        traj.append(row)
        print(
            f"[amazon/{policy}/N={n_cur}] t={w['t']} cat={str(w['category'])[:18]} "
            f"n_adapt={len(adapt)} FLOPs={fr:.2f} lr_ratio={disp['lr_ratio']:.2f} "
            f"acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj, mods, samples, image_tf


def summarize_cell(traj, mods, *, dataset, policy, n_cur):
    char = characterize_soft_lr_traj(traj, mods)
    n_adapt_mean = float(np.mean([r["n_adapt"] for r in traj])) if traj else float("nan")
    return {
        "dataset": dataset,
        "policy": policy,
        "n_cur": int(n_cur),
        "n_adapt_mean": n_adapt_mean,
        **char,
    }


def find_n_stars(cells, *, dataset, policy):
    rows = [
        c
        for c in cells
        if c["dataset"] == dataset and c["policy"] == policy
    ]
    # use mean adapt samples as x-axis
    curve = [
        {
            "n_adapt": c["n_adapt_mean"],
            "n_cur": c["n_cur"],
            "mean_acc_lift": c["mean_acc_lift"],
        }
        for c in rows
    ]
    return n_star_for_lift(curve, n_key="n_adapt", lift_key="mean_acc_lift", eps=0.0)


def plot_board(cells, nstars, path: Path):
    datasets = sorted({c["dataset"] for c in cells})
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), facecolor="#f7f5f1")
    fig.suptitle(
        "Continuous soft-LR · Acc lift vs adapt sample size",
        fontsize=13,
        fontweight="bold",
    )
    for ax, ds in zip(axes, datasets):
        for pol, ls in zip(POLICIES, ["-", "--", ":"]):
            rows = sorted(
                [c for c in cells if c["dataset"] == ds and c["policy"] == pol],
                key=lambda r: r["n_adapt_mean"],
            )
            if not rows:
                continue
            xs = [r["n_adapt_mean"] for r in rows]
            ys = [r["mean_acc_lift"] for r in rows]
            ax.plot(xs, ys, ls, marker="o", label=pol)
            st = nstars.get(f"{ds}:{pol}")
            if st and st.get("n_star") is not None:
                ax.axvline(st["n_star"], color="#999", ls=":", lw=0.8, alpha=0.7)
        ax.axhline(0, color="#bbb", ls="--", lw=0.8)
        ax.set_title(ds)
        ax.set_xlabel("mean adapt samples / window")
        ax.set_ylabel("mean Acc lift")
        ax.legend(frameon=False, fontsize=8)
    fig.text(
        0.5,
        0.02,
        "equal = dense equal-LR · soft_lr = continuous α→LR · soft_sparse = soft + quantile gate",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_latex(cells, nstars, path: Path):
    lines = []
    for c in sorted(cells, key=lambda r: (r["dataset"], r["policy"], r["n_cur"])):
        lines.append(
            f"{c['dataset']} & {c['policy'].replace('_', '\\_')} & {c['n_cur']} & "
            f"{c['n_adapt_mean']:.0f} & {c['mean_acc_lift']:+.3f} & "
            f"{c['mean_flops_rel']:.3f} & {c['mean_lr_ratio']:.2f} & "
            f"{c['mean_alpha_entropy']:.3f} \\\\"
        )
    nstar_lines = []
    for k, st in sorted(nstars.items()):
        ds, pol = k.split(":", 1)
        ns = "---" if st.get("n_star") is None else f"{st['n_star']:.0f}"
        nstar_lines.append(
            f"{ds} & {pol.replace('_', '\\_')} & {ns} & "
            f"{st.get('max_lift', float('nan')):+.3f} & "
            f"{st.get('best_n') if st.get('best_n') is not None else '---'} \\\\"
        )
    tex = (
        "% Continuous soft-LR sample-size curves (Amazon / MSR-VTT)\n"
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{Acc lift vs adapt sample size. Continuous next-stage step sizes "
        "characterized by lr\\_ratio $=\\max_m\\mathrm{LR}_m/\\min_m\\mathrm{LR}_m$ "
        "and $\\alpha$-entropy.}\n"
        "\\label{tab:agod-soft-lr-nsize}\n"
        "\\begin{tabular}{llrrrrrr}\n\\toprule\n"
        "Dataset & policy & $N_{\\mathrm{cur}}$ & $N_{\\mathrm{adapt}}$ & Acc lift & "
        "flops\\_rel & lr\\_ratio & $H(\\alpha)$ \\\\\n\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n\n"
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{$N^\\star$: smallest mean adapt-sample size with mean Acc lift $>0$.}\n"
        "\\label{tab:agod-soft-lr-nstar}\n"
        "\\begin{tabular}{llrrr}\n\\toprule\n"
        "Dataset & policy & $N^\\star$ & max lift & best $N$ \\\\\n\\midrule\n"
        + "\n".join(nstar_lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def write_docs(cells, nstars, path: Path):
    rows = [
        f"| {c['dataset']} | {c['policy']} | {c['n_cur']} | {c['n_adapt_mean']:.0f} | "
        f"{c['mean_acc_lift']:+.3f} | {c['mean_flops_rel']:.3f} | "
        f"{c['mean_lr_ratio']:.2f} | {c['mean_alpha_entropy']:.3f} |"
        for c in sorted(cells, key=lambda r: (r["dataset"], r["policy"], r["n_cur"]))
    ]
    ns = [
        f"| {k.split(':')[0]} | {k.split(':')[1]} | "
        f"{'—' if st.get('n_star') is None else f'{st[\"n_star\"]:.0f}'} | "
        f"{st.get('max_lift', float('nan')):+.3f} |"
        for k, st in sorted(nstars.items())
    ]
    md = f"""# Continuous soft-LR · sample size · evaluation

## Point (fewer switches)

不需要一堆 L3 gate 变体。核心是：

1. **连续 next-stage 步长**（L2 soft LR）：`LR_m = lr0 · (β + (1-β)·α_m·|M|)`
2. **要多少 adapt samples** 才能相对 equal-LR 出现 Acc lift（`N*`)
3. **不必 dense equal adapter**：用 attribution 加权的连续步长；可选一个 sparse gate 做对照

## How to characterize continuous step sizes

| 量 | 含义 |
|---|---|
| `lr_ratio` = max LR / min LR | 模态步长拉开程度（=1 ⇒ 退化成 equal） |
| `lr_std` / `lr_cv` | 步长离散度 |
| `H(α)` | routing 熵（equal ⇒ log\\|M\\|；越小越尖） |
| Acc lift(N) | 效用随 sample size 的曲线 |
| `flops_rel` | soft_lr=1；soft_sparse<1（只省 adapt BWD） |

## How to evaluate

- **相对 equal**：同一 N、同一 window 流，比 mean Acc lift
- **N\\***：升序 N 曲线上，**第一个** mean Acc lift > 0 的 mean `N_adapt`
- **不是越尖越好**：看 lift 与 `lr_ratio` / `H(α)` 是否同向；尖但 lift≤0 ⇒ 过冲
- **sparse 不是必须**：soft_lr 已是非 equal；soft_sparse 只在要省 adapt FLOPs 时加

## Smoke board

| Dataset | policy | N_cur | N_adapt | Acc lift | flops_rel | lr_ratio | H(α) |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

### N* (first mean Acc lift > 0)

| Dataset | policy | N* | max lift |
|---|---|---:|---:|
{chr(10).join(ns)}

```bash
PYTHONPATH=. python3 scripts/run_agod_soft_lr_nsize.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["msrvtt", "amazon"])
    ap.add_argument("--msrvtt-n", nargs="+", type=int, default=list(MSRVTT_N_GRID))
    ap.add_argument("--amazon-n", nargs="+", type=int, default=list(AMAZON_N_GRID))
    ap.add_argument("--policies", nargs="+", default=list(POLICIES))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    cells = []
    traj_dump = {}
    amazon_cache = {"samples": None, "image_tf": None}

    if "msrvtt" in args.datasets:
        for n in args.msrvtt_n:
            for pol in args.policies:
                traj, mods = run_msrvtt_once(device, n_cur=n, policy=pol)
                cell = summarize_cell(traj, mods, dataset="msrvtt", policy=pol, n_cur=n)
                cells.append(cell)
                traj_dump[f"msrvtt:{pol}:N{n}"] = traj

    if "amazon" in args.datasets:
        for n in args.amazon_n:
            for pol in args.policies:
                traj, mods, samples, image_tf = run_amazon_once(
                    device,
                    n_cur=n,
                    policy=pol,
                    samples=amazon_cache["samples"],
                    image_tf=amazon_cache["image_tf"],
                )
                amazon_cache["samples"] = samples
                amazon_cache["image_tf"] = image_tf
                cell = summarize_cell(traj, mods, dataset="amazon", policy=pol, n_cur=n)
                cells.append(cell)
                traj_dump[f"amazon:{pol}:N{n}"] = traj

    nstars = {}
    for ds in sorted({c["dataset"] for c in cells}):
        for pol in args.policies:
            nstars[f"{ds}:{pol}"] = find_n_stars(cells, dataset=ds, policy=pol)

    payload = {
        "agod_version": "0.1.0",
        "focus": [
            "continuous_soft_lr_step_sizes",
            "sample_size_to_positive_lift",
            "sparse_vs_dense_equal",
        ],
        "policies": list(args.policies),
        "design": {
            "equal": "dense equal LR (B1); flops_rel=1",
            "soft_lr": "continuous α→LR_m; all mods active; characterize via lr_ratio/H(α)",
            "soft_sparse": "soft_lr + quantile hard-gate; sparse adapt-BWD",
            "n_star": "smallest mean N_adapt with mean Acc lift > 0 on ascending curve",
        },
        "cells": cells,
        "n_stars": nstars,
        "trajectory": traj_dump,
    }
    (OUT / "agod_soft_lr_nsize.json").write_text(json.dumps(payload, indent=2))
    plot_board(cells, nstars, OUT / "AGOD_Soft_LR_NSize_Board.png")
    write_latex(cells, nstars, OUT / "AGOD_soft_lr_nsize_tables_only.tex")
    write_docs(cells, nstars, OUT / "README.md")
    write_latex(cells, nstars, DOCS / "AGOD_soft_lr_nsize_tables_only.tex")
    write_docs(cells, nstars, DOCS / "AGOD_soft_lr_nsize.md")
    # mirror board to artifacts
    import shutil

    shutil.copy2(OUT / "AGOD_Soft_LR_NSize_Board.png", ART / "AGOD_Soft_LR_NSize_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    print("\n=== N* summary ===", flush=True)
    for k, st in sorted(nstars.items()):
        print(f"  {k}: n_star={st.get('n_star')} max_lift={st.get('max_lift'):+.3f}", flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
