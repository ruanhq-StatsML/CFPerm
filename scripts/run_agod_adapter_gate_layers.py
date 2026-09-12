#!/usr/bin/env python3
"""Adapter layers + dynamic L3 hard-gate on Amazon + MSR-VTT.

Fixed routing prior; vary only L3:
  none | fixed | quantile | ema_theta | hysteresis | random

  L0 sensor → L1 EMA α → L2 soft LR(α) → L3 hard adapt-gate (proj BWD)

  PYTHONPATH=. python3 scripts/run_agod_adapter_gate_layers.py
  PYTHONPATH=. python3 scripts/run_agod_adapter_gate_layers.py --datasets msrvtt
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

import agod
from agod.adapter import (
    ADAPTER_LAYERS,
    GATE_POLICIES,
    AdapterConfig,
    LayeredAdapter,
    characterize_gate_traj,
    flops_rel_proj,
)
from agod.lr_controller import EMARouter, softmax_scores
from agod.shift import decompose_hybrid, decompose_mmd, residual_concept

OUT = ROOT / "results" / "agod_adapter_layers"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_adapter_layers")

SEED = 2026
EMA = 0.40
TAU = 0.30
BETA = 0.10
HOLD = 0.35
GATE_VARIANTS = ("none", "fixed", "quantile", "ema_theta", "hysteresis", "random")

# MSR-VTT knobs
MSRVTT_MODS = ["video", "text", "audio"]
N_REF, N_CUR, T_WIN = 400, 320, 6
STEPS, BATCH, LR0, FUSE = 35, 64, 3e-3, 128


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


def make_stream_y(y):
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
        k1 = int(N_CUR * p1)
        cur = np.concatenate([take(rest0, N_CUR - k1, t), take(rest1, k1, t)])
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
                n_estimators=50,
                max_depth=6,
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
        self.gate = nn.Linear(FUSE * len(mods), len(mods))
        self.head = nn.Sequential(nn.Linear(FUSE, 64), nn.ReLU(), nn.Linear(64, 2))

    def forward(self, batch):
        zs = [self.projs[m](batch[m]) for m in self.mods]
        g = torch.softmax(self.gate(torch.cat(zs, 1)), 1)
        fused = sum(g[:, i : i + 1] * zs[i] for i in range(len(self.mods)))
        return self.head(fused)

    def param_groups(self):
        g = {m: list(self.projs[m].parameters()) for m in self.mods}
        g["shared"] = list(self.gate.parameters()) + list(self.head.parameters())
        return g


def build_optim(model, mods):
    g = model.param_groups()
    return torch.optim.AdamW(
        [{"params": g[m], "lr": LR0, "name": m} for m in mods]
        + [{"params": g["shared"], "lr": LR0, "name": "shared"}]
    )


def set_lrs(opt, mult):
    for g in opt.param_groups:
        g["lr"] = LR0 * float(mult.get(g.get("name", "shared"), 1.0))


def to_batch(feats, y, idx, mods, device):
    return (
        {m: torch.from_numpy(feats[m][idx]).to(device) for m in mods},
        torch.from_numpy(y[idx].astype(np.int64)).to(device),
    )


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
    n_h = max(40, min(len(idx) // 2, int(len(idx) * HOLD)))
    return idx[perm[n_h:]], idx[perm[:n_h]]


def train_window(model, opt, feats, y, idx, mods, device, *, lr_mult, active):
    model.train()
    set_lrs(opt, lr_mult)
    crit = nn.CrossEntropyLoss()
    losses = []
    rng = np.random.default_rng(SEED + int(np.asarray(idx).sum()) % 100000)
    t0 = time.perf_counter()
    for _ in range(STEPS):
        sel = rng.choice(idx, size=min(BATCH, len(idx)), replace=False)
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


def run_msrvtt_gate(device, gate_policy: str):
    feats, y = load_msrvtt()
    mods = list(MSRVTT_MODS)
    stream = make_stream_y(y)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = Fusion({m: feats[m].shape[1] for m in mods}, mods).to(device)
    opt = build_optim(model, mods)
    router = EMARouter(mods, ema=EMA)
    adapter = LayeredAdapter(
        mods,
        AdapterConfig(beta=BETA, theta=0.28, p_keep_rand=0.67),
        gate_policy=gate_policy,
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
        dec = adapter.step(alpha, seed=SEED, t=w["t"])
        pre = eval_acc(model, feats, y, hold, mods, device)
        loss, wall = train_window(
            model,
            opt,
            feats,
            y,
            adapt,
            mods,
            device,
            lr_mult=dec.lr_mult,
            active=dec.active,
        )
        post = eval_acc(model, feats, y, hold, mods, device)
        dacc = post["acc"] - pre["acc"]
        fr = flops_rel_proj(dec.active, mods)
        keep = N_REF // 2
        ref_idx = np.concatenate([ref_idx[-keep:], adapt[: min(keep, len(adapt))]])
        ref_b, y_ref = blocks(feats, y, ref_idx, mods)
        row = {
            "t": w["t"],
            "p1": w["p1"],
            "alpha": dec.alpha,
            "lr_mult": dec.lr_mult,
            "active": dec.active,
            "theta_used": dec.theta_used,
            "flops_rel": fr,
            "wall_ms": wall,
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": dacc,
            "train_loss": loss,
        }
        traj.append(row)
        print(
            f"[msrvtt/{gate_policy}] t={w['t']} θ={dec.theta_used:.3f} "
            f"sel={ {m: int(dec.active[m]) for m in mods} } "
            f"FLOPs={fr:.2f} | acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj, mods


def run_amazon_gate(device, gate_policy: str):
    import run_agod_amazon_modality_lr as base

    mods = list(base.MODS)
    shards = sorted(base.SHARD_DIR.glob("*.tar.gz"))
    if not shards:
        raise SystemExit(f"No Amazon shards in {base.SHARD_DIR}")
    print(f"loading {len(shards)} amazon shards…", flush=True)
    samples = base.load_shards(shards)
    print(
        f"samples={len(samples)} cats={Counter(s['category'] for s in samples).most_common(5)}",
        flush=True,
    )
    image_tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
        ]
    )
    stream = base.make_stream(samples)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = base.AmazonAGOD().to(device)
    opt = base.build_optim(model)
    router = EMARouter(mods, ema=EMA)
    adapter = LayeredAdapter(
        mods,
        AdapterConfig(
            beta=BETA,
            theta=0.40,
            theta_on=0.45,
            theta_off=0.35,
            p_keep_rand=0.67,
        ),
        gate_policy=gate_policy,
    )
    ref = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref, device, image_tf)
    traj = []
    for w in stream["windows"]:
        rng = np.random.default_rng(SEED + 11 * w["t"])
        rows = list(w["rows"])
        rng.shuffle(rows)
        n_h = max(24, min(len(rows) // 2, int(len(rows) * 0.40)))
        hold, adapt = rows[:n_h], rows[n_h:]
        cur_b, y_cur = base.extract_blocks(model, adapt, device, image_tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        raw = msg.alpha["B3"]
        alpha = router.update(raw)
        dec = adapter.step(alpha, seed=SEED, t=w["t"])
        pre = _amazon_eval(model, hold, device, image_tf, base)
        loss, wall = _amazon_train(model, opt, adapt, device, image_tf, dec, mods, base)
        post = _amazon_eval(model, hold, device, image_tf, base)
        dacc = post["acc"] - pre["acc"]
        fr = flops_rel_proj(dec.active, mods)
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
            "alpha": dec.alpha,
            "lr_mult": dec.lr_mult,
            "active": dec.active,
            "theta_used": dec.theta_used,
            "flops_rel": fr,
            "wall_ms": wall,
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": dacc,
            "train_loss": loss,
        }
        traj.append(row)
        print(
            f"[amazon/{gate_policy}] t={w['t']} cat={str(w['category'])[:18]} "
            f"θ={dec.theta_used:.3f} "
            f"sel={ {m: int(dec.active[m]) for m in mods} } "
            f"FLOPs={fr:.2f} | acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj, mods


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


def _amazon_train(model, opt, rows, device, image_tf, dec, mods, base):
    from torch.utils.data import DataLoader

    model.train()
    # freeze encoder if present
    if hasattr(model, "image_encoder"):
        model.image_encoder.eval()
    base.set_lrs(opt, dec.lr_mult)
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=base.BATCH, shuffle=True)
    crit = nn.CrossEntropyLoss()
    losses, it = [], iter(loader)
    t0 = time.perf_counter()
    steps = getattr(base, "STEPS_PER_WIN", 8)
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
            if not dec.active.get(m, True):
                for p in groups[m]:
                    p.grad = None
        opt.step()
        losses.append(float(loss.item()))
    return float(np.mean(losses)), (time.perf_counter() - t0) * 1000.0


def summarize_ds(traj, mods):
    char = characterize_gate_traj(traj, mods)
    return {
        "mean_acc_lift": char.get("mean_acc_lift", float("nan")),
        "mean_acc_post": char.get("mean_acc_post", float("nan")),
        "mean_flops_rel": char.get("mean_flops_rel", float("nan")),
        "mean_cost_utility": char.get("mean_cost_utility", float("nan")),
        "switches_per_window": char.get("switches_per_window", float("nan")),
        "mean_theta_used": char.get("mean_theta_used", float("nan")),
        "frac_active": char.get("frac_active", {}),
        "mean_active_mods": char.get("mean_active_mods", float("nan")),
        "n_windows": char.get("n_windows", len(traj)),
    }


def plot_board(all_summary, path: Path):
    datasets = list(all_summary.keys())
    gates = [g for g in GATE_VARIANTS if g in next(iter(all_summary.values()))]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.8), facecolor="#f7f5f1")
    fig.suptitle(
        "AGOD adapter layers · dynamic L3 gate (Amazon / MSR-VTT)",
        fontsize=13,
        fontweight="bold",
    )
    x = np.arange(len(gates))
    w = 0.35
    for i, ds in enumerate(datasets):
        lifts = [all_summary[ds][g]["mean_acc_lift"] for g in gates]
        flops = [all_summary[ds][g]["mean_flops_rel"] for g in gates]
        axes[0].bar(x + (i - 0.5) * w, lifts, w, label=ds)
        axes[1].bar(x + (i - 0.5) * w, flops, w, label=ds)
    axes[0].axhline(0, color="#999", ls="--", lw=0.8)
    axes[0].set_title("mean Acc lift")
    axes[1].set_title("mean flops_rel")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(gates, rotation=20, ha="right")
        ax.legend(frameon=False, fontsize=8)
    fig.text(
        0.5,
        0.02,
        "L2 soft LR(α) always on · L3 hard adapt-gate variants · FWD always paid",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_latex(all_summary, path: Path):
    gates = [g for g in GATE_VARIANTS if g in next(iter(all_summary.values()))]
    rows = []
    for ds, sm in all_summary.items():
        for g in gates:
            s = sm[g]
            rows.append(
                f"{ds} & {g.replace('_', '\\_')} & {s['mean_acc_lift']:+.3f} & "
                f"{s['mean_flops_rel']:.3f} & {s['mean_cost_utility']:+.3f} & "
                f"{s['switches_per_window']:.2f} & {s['mean_theta_used']:.3f} \\\\"
            )
    tex = (
        "% AGOD adapter layers + dynamic L3 gate (Amazon / MSR-VTT)\n"
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{Adapter-layer smoke: fixed routing prior, vary only L3 hard-gate. "
        "L2 soft LR$(\\alpha)$ always on. FWD inference always paid.}\n"
        "\\label{tab:agod-adapter-gate-layers}\n"
        "\\begin{tabular}{llccccc}\n\\toprule\n"
        "Dataset & L3 gate & Acc lift & flops\\_rel & util & switches/win & "
        "$\\theta_t$ \\\\\n\\midrule\n"
        + "\n".join(rows)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.write_text(tex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["msrvtt", "amazon"],
        choices=["msrvtt", "amazon"],
    )
    ap.add_argument("--gates", nargs="+", default=list(GATE_VARIANTS), choices=list(GATE_POLICIES))
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device} agod={agod.__version__}", flush=True)
    print(f"layers={ADAPTER_LAYERS}", flush=True)
    print(f"gates={args.gates}", flush=True)

    all_summary, all_traj = {}, {}
    for ds in args.datasets:
        all_summary[ds], all_traj[ds] = {}, {}
        for g in args.gates:
            print(f"\n===== {ds} / L3={g} =====", flush=True)
            if ds == "msrvtt":
                traj, mods = run_msrvtt_gate(device, g)
            else:
                traj, mods = run_amazon_gate(device, g)
            summary = summarize_ds(traj, mods)
            all_summary[ds][g] = summary
            all_traj[ds][g] = traj
            print(
                f"→ lift={summary['mean_acc_lift']:+.3f} "
                f"flops={summary['mean_flops_rel']:.3f} "
                f"util={summary['mean_cost_utility']:+.3f} "
                f"switches/win={summary['switches_per_window']:.2f}",
                flush=True,
            )

    payload = {
        "agod_version": agod.__version__,
        "adapter_layers": list(ADAPTER_LAYERS),
        "gate_policies": list(args.gates),
        "design": {
            "L0": "shift sensors (MMD/PO on MSR-VTT; MSG on Amazon)",
            "L1": "EMA Softmax alpha state",
            "L2": "soft LR adapter LR_m=lr0*(beta+(1-beta)*alpha_m*|M|)",
            "L3": "hard adapt-gate on proj BWD; FWD always on",
        },
        "summary": all_summary,
        "trajectory": all_traj,
    }
    jp = OUT / "agod_adapter_gate_layers.json"
    jp.write_text(json.dumps(payload, indent=2, default=float))
    board = OUT / "AGOD_Adapter_Gate_Layers_Board.png"
    plot_board(all_summary, board)
    tex = OUT / "AGOD_adapter_gate_layers_tables_only.tex"
    write_latex(all_summary, tex)
    (DOCS / tex.name).write_text(tex.read_text())

    md_rows = []
    for ds, sm in all_summary.items():
        for g in args.gates:
            s = sm[g]
            md_rows.append(
                f"| {ds} | {g} | {s['mean_acc_lift']:+.3f} | "
                f"{s['mean_flops_rel']:.3f} | {s['mean_cost_utility']:+.3f} | "
                f"{s['switches_per_window']:.2f} |"
            )
    note = DOCS / "AGOD_adapter_gate_layers.md"
    note.write_text(
        "# AGOD adapter layers + dynamic L3 gate\n\n"
        "## How many layers / how to characterize\n\n"
        "| Layer | Role | Adjustable knobs |\n"
        "|---|---|---|\n"
        "| L0 sensor | cov/concept scores | MMD / PO / MSG |\n"
        "| L1 state | EMA Softmax α | ema, τ |\n"
        "| L2 soft LR | continuous next-stage step sizes | β, lr0 |\n"
        "| L3 hard gate | sparse adapt (proj BWD only) | θ / quantile / hyst / rand |\n\n"
        "FWD inference always paid. Efficiency claim is **adapt FLOPs**, not latency.\n\n"
        "## Dynamic L3 policies\n\n"
        "- `none` — soft LR only\n"
        "- `fixed` — α ≥ θ\n"
        "- `quantile` — drop bottom-α mass\n"
        "- `ema_theta` — θ tracks EMA(mean α)\n"
        "- `hysteresis` — on/off thresholds (less chatter)\n"
        "- `random` — non-attribution sparsity control\n\n"
        "## Smoke board (Amazon + MSR-VTT)\n\n"
        "| Dataset | L3 gate | Acc lift | flops_rel | util | switches/win |\n"
        "|---|---|---:|---:|---:|---:|\n"
        + "\n".join(md_rows)
        + "\n\n```bash\nPYTHONPATH=. python3 scripts/run_agod_adapter_gate_layers.py\n```\n"
    )
    (OUT / "README.md").write_text(note.read_text())
    for p in [jp, board, tex, note]:
        (ART / Path(p).name).write_bytes(Path(p).read_bytes())
    print("\n=== BOARD ===", flush=True)
    print(note.read_text(), flush=True)


if __name__ == "__main__":
    main()
