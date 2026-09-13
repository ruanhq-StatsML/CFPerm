#!/usr/bin/env python3
"""MSR-VTT 3-mod: FSDS + instance-disc + proto-drift + grad-memory landscape.

  L0 sensors : PO/MMD + NN instance discrimination + prototype drift
  L1 routing : compose_proto_fsds -> Softmax -> EMA alpha (proportions)
  L2 actuator: next_step_lr(alpha or alpha_hat) * landscape_adapt_gain

Policies: equal | fsds | fsds_disc_proto | fused_landscape

  PYTHONPATH=. python3 scripts/run_agod_fsds_landscape_smoke.py
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.grad_memory import GradProjMemory, landscape_adapt_gain
from agod.instance_disc import disc_as_concept, nn_discrimination_scores
from agod.lr_controller import (
    EMARouter,
    alpha_to_lr,
    common_dim_grad_signature,
    softmax_scores,
)
from agod.mmd import rbf_mmd2
from agod.next_step import next_step_lr, predict_next_alpha, proportion_report
from agod.proto_drift import ModalityPrototypeBank, compose_proto_fsds
from agod.shift import residual_concept

OUT = ROOT / "results" / "agod_fsds_landscape"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_fsds_landscape")

MODS = ["video", "text", "audio"]
SEED = 2026
N_REF, N_CUR, T_WIN = 320, 280, 6
HOLD, STEPS, BATCH = 0.35, 28, 64
LR0, EMA, TAU, BETA = 3e-3, 0.40, 0.30, 0.10
FUSE = 128
POLICIES = ("equal", "fsds", "fsds_disc_proto", "fused_landscape")


def load_pack(pack: str = "packed"):
    d = ROOT / "data" / "msrvtt" / pack
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
    return feats, y.astype(np.int64), d


def make_stream(y):
    rng = np.random.default_rng(SEED)
    idx = rng.permutation(len(y))
    ref = idx[:N_REF]
    rest = idx[N_REF:]
    wins = []
    for t in range(T_WIN):
        sl = rest[t * N_CUR : (t + 1) * N_CUR]
        if len(sl) < 80:
            break
        wins.append({"t": t, "idx": sl})
    return {"ref": ref, "windows": wins}


def blocks(feats, y, idx):
    return {m: feats[m][idx] for m in MODS}, y[idx].astype(float)


def mmd_scores(b0, b1):
    out = {}
    for m in MODS:
        try:
            out[m] = float(max(rbf_mmd2(b0[m], b1[m], seed=SEED), 0.0))
        except TypeError:
            out[m] = float(max(rbf_mmd2(b0[m], b1[m]), 0.0))
        except Exception:
            out[m] = 0.0
    return out


def po_scores(b0, y0, b1, y1):
    return {
        m: float(max(residual_concept(b0[m], y0, b1[m], y1, seed=SEED + i), 0.0))
        for i, m in enumerate(MODS)
    }


class MSRFusion(nn.Module):
    def __init__(self, dims):
        super().__init__()
        self.video_proj = nn.Sequential(nn.Linear(dims["video"], FUSE), nn.ReLU(), nn.Dropout(0.1))
        self.text_proj = nn.Sequential(nn.Linear(dims["text"], FUSE), nn.ReLU(), nn.Dropout(0.1))
        self.audio_proj = nn.Sequential(nn.Linear(dims["audio"], FUSE), nn.ReLU(), nn.Dropout(0.1))
        self.gate = nn.Linear(FUSE * 3, 3)
        self.head = nn.Sequential(nn.Linear(FUSE, 64), nn.ReLU(), nn.Linear(64, 2))

    def forward(self, batch):
        v = self.video_proj(batch["video"])
        t = self.text_proj(batch["text"])
        a = self.audio_proj(batch["audio"])
        g = torch.softmax(self.gate(torch.cat([v, t, a], 1)), 1)
        fused = g[:, 0:1] * v + g[:, 1:2] * t + g[:, 2:3] * a
        return self.head(fused)

    def param_groups(self):
        return {
            "video": list(self.video_proj.parameters()),
            "text": list(self.text_proj.parameters()),
            "audio": list(self.audio_proj.parameters()),
            "shared": list(self.gate.parameters()) + list(self.head.parameters()),
        }


def build_optim(model):
    g = model.param_groups()
    return torch.optim.AdamW(
        [
            {"params": g["video"], "lr": LR0, "name": "video"},
            {"params": g["text"], "lr": LR0, "name": "text"},
            {"params": g["audio"], "lr": LR0, "name": "audio"},
            {"params": g["shared"], "lr": LR0, "name": "shared"},
        ]
    )


def set_lrs(opt, mult):
    for g in opt.param_groups:
        g["lr"] = LR0 * float(mult.get(g.get("name", "shared"), 1.0))


def to_batch(feats, y, idx, device):
    return (
        {m: torch.from_numpy(feats[m][idx]).to(device) for m in MODS},
        torch.from_numpy(y[idx].astype(np.int64)).to(device),
    )


@torch.no_grad()
def eval_acc(model, feats, y, idx, device):
    model.eval()
    if len(idx) < 8:
        return float("nan")
    logits, ys = [], []
    for s in range(0, len(idx), BATCH):
        batch, yy = to_batch(feats, y, idx[s : s + BATCH], device)
        logits.append(model(batch).cpu().numpy())
        ys.append(yy.cpu().numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    return float((L.argmax(1) == Y).mean())


def split_hold(idx, seed):
    rng = np.random.default_rng(seed)
    p = rng.permutation(len(idx))
    n_h = max(24, int(len(idx) * HOLD))
    return idx[p[n_h:]], idx[p[:n_h]]


def grad_sigs(model, feats, y, idx, device):
    model.train()
    model.zero_grad(set_to_none=True)
    sel = idx[: min(BATCH, len(idx))]
    batch, yy = to_batch(feats, y, sel, device)
    loss = nn.CrossEntropyLoss()(model(batch), yy)
    loss.backward()
    groups = model.param_groups()
    sigs = {m: common_dim_grad_signature(groups[m]) for m in MODS}
    model.zero_grad(set_to_none=True)
    return sigs


def train_window(model, opt, feats, y, idx, device, *, lr_mult):
    model.train()
    set_lrs(opt, lr_mult)
    crit = nn.CrossEntropyLoss()
    losses = []
    t0 = time.perf_counter()
    for _ in range(STEPS):
        sel = idx[np.random.randint(0, len(idx), size=min(BATCH, len(idx)))]
        batch, yy = to_batch(feats, y, sel, device)
        loss = crit(model(batch), yy)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        losses.append(float(loss.item()))
    return float(np.mean(losses)), (time.perf_counter() - t0) * 1000.0


def run_policy(feats, y, stream, device, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    dims = {m: feats[m].shape[1] for m in MODS}
    model = MSRFusion(dims).to(device)
    opt = build_optim(model)
    ref_idx = stream["ref"]
    ref_b, y_ref = blocks(feats, y, ref_idx)

    router = EMARouter(MODS, ema=EMA)
    proto = ModalityPrototypeBank(MODS, ema=0.85)
    proto.update(ref_b, y_ref)
    gmem = GradProjMemory(MODS, capacity=6)
    alpha_hist = []
    sensors_prev = None
    traj = []

    for w in stream["windows"]:
        adapt, hold = split_hold(w["idx"], SEED + 11 * w["t"])
        cur_b, y_cur = blocks(feats, y, adapt)
        po = po_scores(ref_b, y_ref, cur_b, y_cur)
        mmd = mmd_scores(ref_b, cur_b)
        disc_raw = nn_discrimination_scores(ref_b, cur_b, MODS, k=5, seed=SEED + w["t"])
        disc = disc_as_concept(disc_raw, MODS)
        proto_d = proto.drift(cur_b, y_cur)

        if policy == "equal":
            raw = {m: 1.0 / len(MODS) for m in MODS}
            alpha = dict(raw)
        elif policy == "fsds":
            raw = compose_proto_fsds(
                po, mmd, {m: 0.0 for m in MODS}, {m: 0.0 for m in MODS}, MODS
            )
            alpha = router.update(softmax_scores(raw, MODS, TAU))
        else:
            # fsds_disc_proto and fused_landscape share the same alpha sensor
            raw = compose_proto_fsds(po, mmd, proto_d, disc, MODS)
            alpha = router.update(softmax_scores(raw, MODS, TAU))

        sigs = grad_sigs(model, feats, y, adapt, device)
        resid = gmem.residual_energy(sigs)
        land = gmem.landscape()
        gain = landscape_adapt_gain(resid, land, MODS)
        gmem.push(sigs)

        sensors_now = {"g": raw, "po": po, "mmd": mmd, "disc": disc, "proto_d": proto_d}
        alpha_hat = predict_next_alpha(
            alpha_hist + [alpha],
            MODS,
            sensors_now=sensors_now,
            sensors_prev=sensors_prev,
        )
        prop = proportion_report(alpha, alpha_hat, MODS)

        if policy == "fused_landscape":
            lr_mult = next_step_lr(alpha_hat, MODS, landscape_gain=gain, beta=BETA)
        else:
            lr_mult = alpha_to_lr(alpha, MODS, beta=BETA)

        pre = eval_acc(model, feats, y, hold, device)
        loss, wall = train_window(model, opt, feats, y, adapt, device, lr_mult=lr_mult)
        post = eval_acc(model, feats, y, hold, device)

        proto.update(cur_b, y_cur)
        keep = N_REF // 2
        ref_idx = np.concatenate([ref_idx[-keep:], adapt[: min(keep, len(adapt))]])
        ref_b, y_ref = blocks(feats, y, ref_idx)

        alpha_hist.append(alpha)
        sensors_prev = sensors_now
        prop_mae = float(prop.get("prop_mae", prop.get("mae", 0.0)))
        top_mod = prop.get("top_mod", max(MODS, key=lambda m: alpha[m]))
        erank = float(land.get("effective_rank", land.get("erank", len(MODS))))
        row = {
            "t": int(w["t"]),
            "alpha": {m: float(alpha[m]) for m in MODS},
            "alpha_hat_next": {m: float(alpha_hat[m]) for m in MODS},
            "prop_mae": prop_mae,
            "top_mod": top_mod,
            "po": {m: float(po[m]) for m in MODS},
            "mmd": {m: float(mmd[m]) for m in MODS},
            "disc": {m: float(disc_raw[m]) for m in MODS},
            "proto_d": {m: float(proto_d[m]) for m in MODS},
            "residual": {m: float(resid[m]) for m in MODS},
            "erank": erank,
            "lr_mult": {k: float(v) for k, v in lr_mult.items()},
            "acc_pre": pre,
            "acc_post": post,
            "acc_lift": post - pre,
            "train_loss": loss,
            "wall_ms": wall,
        }
        traj.append(row)
        print(
            f"[msrvtt/{policy}] t={w['t']} "
            f"Acc {pre:.3f}->{post:.3f} (d={row['acc_lift']:+.3f}) "
            f"a={{{', '.join(f'{m}:{alpha[m]:.2f}' for m in MODS)}}} "
            f"hat_mae={prop_mae:.3f} erank={erank:.2f}",
            flush=True,
        )
    return traj


def summarize(traj, policy):
    return {
        "dataset": "msrvtt",
        "policy": policy,
        "n_windows": len(traj),
        "mean_acc_lift": float(np.mean([r["acc_lift"] for r in traj])),
        "mean_acc_post": float(np.mean([r["acc_post"] for r in traj])),
        "mean_prop_mae": float(np.mean([r["prop_mae"] for r in traj])),
        "mean_erank": float(np.mean([r["erank"] for r in traj])),
        "mean_alpha": {m: float(np.mean([r["alpha"][m] for r in traj])) for m in MODS},
    }


def plot_board(cells, trajs, path: Path):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), facecolor="#f7f5f1")
    names = [c["policy"] for c in cells]
    ax = axes[0]
    x = np.arange(len(names))
    ax.bar(x, [c["mean_acc_lift"] for c in cells], color="#2f5d50")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=18, ha="right")
    ax.axhline(0, color="#999", ls=":", lw=0.9)
    ax.set_ylabel("Acc lift")
    ax.set_title("MSR-VTT: Acc up by policy")

    ax = axes[1]
    for pol, traj in trajs.items():
        ax.plot([r["t"] for r in traj], [r["acc_lift"] for r in traj], marker="o", label=pol, lw=1.5)
    ax.axhline(0, color="#999", ls=":", lw=0.9)
    ax.set_xlabel("window t")
    ax.set_ylabel("Acc lift")
    ax.set_title("Acc trajectory")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[2]
    target = "fused_landscape" if "fused_landscape" in trajs else names[-1]
    mean_a = cells[[c["policy"] for c in cells].index(target)]["mean_alpha"]
    ax.bar(MODS, [mean_a[m] for m in MODS], color=["#b85c38", "#3d5a80", "#6a994e"])
    ax.set_ylim(0, 1)
    ax.set_ylabel("mean alpha")
    ax.set_title(f"3-mod proportions ({target})")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(cells, path: Path):
    rows = []
    for c in cells:
        a = " / ".join(f"{m}:{c['mean_alpha'][m]:.2f}" for m in MODS)
        rows.append(
            f"| `{c['policy']}` | {c['mean_acc_lift']:+.3f} | {c['mean_acc_post']:.3f} | "
            f"{c['mean_prop_mae']:.3f} | {c['mean_erank']:.2f} | {a} |"
        )
    best = max(cells, key=lambda c: c["mean_acc_lift"])
    md = f"""# FSDS + instance-disc + proto-drift + grad-memory (MSR-VTT)

## Logic

1. **FSDS hybrid**: PO (concept) - MMD (covariate) as the base attribution score.
2. **Nonparametric instance discrimination**: kNN domain purity on each modality block
   -> feature-level separability sensor (no deep NCE).
3. **Prototype drift**: EMA modality centroids; d=1-cos(proto, batch_mean);
   high PO + high proto-drift => true concept move; high MMD + low proto-drift => mush.
4. **Grad-projection memory + landscape**: bank of common-dim grad signatures;
   residual energy + Gram erank -> LR gains (new direction up, collinear basin down).
5. **Proportion + next-step**: alpha_video/alpha_text/alpha_audio now; a-hat from EMA+sensor delta;
   fused policy uses a-hat * landscape_gain as the adapt step.

## Smoke (MSR-VTT packed, 3-mod)

| policy | Acc up | Acc post | prop MAE | erank | alpha (v/t/a) |
|---|---:|---:|---:|---:|---|
{chr(10).join(rows)}

Best Acc up: **`{best['policy']}`** ({best['mean_acc_lift']:+.3f}).

```bash
PYTHONPATH=. python3 scripts/run_agod_fsds_landscape_smoke.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--pack", default="packed")
    ap.add_argument("--policies", nargs="+", default=list(POLICIES), choices=list(POLICIES))
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    print("loading MSR-VTT pack...", flush=True)
    feats, y, d = load_pack(args.pack)
    print(f"n={len(y)} dims={[feats[m].shape[1] for m in MODS]} @ {d}", flush=True)
    stream = make_stream(y)
    device = args.device

    cells, trajs = [], {}
    for pol in args.policies:
        traj = run_policy(feats, y, stream, device, pol)
        trajs[pol] = traj
        cells.append(summarize(traj, pol))

    payload = {
        "agod_version": "0.1.0",
        "focus": "FSDS + instance-disc + proto-drift + grad-memory landscape -> next adapt",
        "mods": MODS,
        "rule": {
            "sensors": "PO, MMD, NN-disc, proto-drift, grad-residual, erank",
            "proportions": "alpha_m = EMA(Softmax(compose_proto_fsds(...)))",
            "next_step": "a_hat=predict_next_alpha; LR=next_step_lr(a_hat)*landscape_gain",
        },
        "cells": cells,
        "trajectory": trajs,
    }
    (OUT / "agod_fsds_landscape.json").write_text(json.dumps(payload, indent=2))
    plot_board(cells, trajs, OUT / "AGOD_FSDS_Landscape_Board.png")
    write_docs(cells, OUT / "README.md")
    write_docs(cells, DOCS / "AGOD_fsds_landscape.md")
    shutil.copy2(OUT / "AGOD_FSDS_Landscape_Board.png", ART / "AGOD_FSDS_Landscape_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")

    print("=== summary ===")
    for c in cells:
        print(
            f"  {c['policy']}: Acc+={c['mean_acc_lift']:+.3f} Acc_post={c['mean_acc_post']:.3f} "
            f"prop_mae={c['mean_prop_mae']:.3f} erank={c['mean_erank']:.2f} "
            f"a=" + "/".join(f"{c['mean_alpha'][m]:.2f}" for m in MODS)
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
