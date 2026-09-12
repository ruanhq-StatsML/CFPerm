#!/usr/bin/env python3
"""Amazon online: drift-vs-noise gate + smooth EMA routing.

Policies
  equal               — uniform α / LR
  msg_softmax         — Softmax(AUC·VIMP + γ·PO) + EMA  (old MSG)
  concept_minus_cov   — Softmax(concept − covariate) + EMA
  smooth_drift_noise  — PO gated by uni-Acc + VIMP/Fisher, then
                        g=norm(w1·PO_g + w2·MMD + w3·VIMP),
                        α=EMA(softmax(g/τ)),
                        w=w0(β+(1-β)α|M|)

Reports holdout Acc↑ and Brier/MSE.

  PYTHONPATH=. python3 scripts/run_agod_smooth_drift_noise.py
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from torch.utils.data import DataLoader
from torchvision import transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import run_agod_amazon_modality_lr as base
from agod.smooth_router import (
    SmoothDriftNoiseRouter,
    SmoothRouterConfig,
    normalize_nonneg,
    softmax_tau,
    weights_from_alpha,
)

OUT = ROOT / "results" / "agod_smooth_drift_noise"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_smooth_drift_noise")

SEED = base.SEED
MODS = list(base.MODS)
POLICIES = ("equal", "msg_softmax", "concept_minus_cov", "smooth_drift_noise")
EMA, TAU, BETA = 0.40, 0.30, 0.10
HOLD = 0.40
STEPS = 10  # CPU-friendly; base uses STEPS_PER_WIN=8
BATCH = base.BATCH


def split_rows(rows, seed):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    n_h = max(24, min(len(rows) // 2, int(len(rows) * HOLD)))
    return [rows[i] for i in idx[n_h:]], [rows[i] for i in idx[:n_h]]


def probs_mse(logits: np.ndarray, y: np.ndarray) -> float:
    z = logits - logits.max(axis=1, keepdims=True)
    e = np.exp(z)
    p = e / np.clip(e.sum(axis=1, keepdims=True), 1e-12, None)
    oh = np.zeros_like(p)
    oh[np.arange(len(y)), y.astype(int)] = 1.0
    return float(np.mean((p - oh) ** 2))


@torch.no_grad()
def eval_hold(model, rows, device, image_tf):
    model.eval()
    if len(rows) < 8:
        return {"acc": float("nan"), "mse": float("nan"), "auc": float("nan"), "n": len(rows)}
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=BATCH, shuffle=False)
    logits, ys = [], []
    for imgs, txts, y in loader:
        logits.append(model(imgs.to(device), txts, return_mods=False).cpu().numpy())
        ys.append(y.numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    out = {
        "acc": float(accuracy_score(Y, L.argmax(1))),
        "mse": probs_mse(L, Y),
        "auc": float("nan"),
        "n": int(len(Y)),
    }
    if len(np.unique(Y)) > 1:
        try:
            p = torch.softmax(torch.from_numpy(L), 1).numpy()[:, 1]
            out["auc"] = float(roc_auc_score(Y, p))
        except Exception:
            pass
    return out


def unimodal_acc(Xtr, ytr, Xte, yte, *, seed: int) -> float:
    if len(Xtr) < 16 or len(Xte) < 8 or len(np.unique(ytr)) < 2:
        return 0.5
    try:
        clf = RandomForestClassifier(
            n_estimators=40,
            max_depth=6,
            min_samples_leaf=2,
            n_jobs=1,
            random_state=seed,
        )
        clf.fit(Xtr, ytr.astype(int))
        return float(accuracy_score(yte.astype(int), clf.predict(Xte)))
    except Exception:
        return 0.5


def measure_uni_acc(b_adapt, y_adapt, b_hold, y_hold, mods, *, seed: int) -> dict[str, float]:
    return {
        m: unimodal_acc(
            b_adapt[m], y_adapt, b_hold[m], y_hold, seed=seed + 11 * i
        )
        for i, m in enumerate(mods)
    }


def mmd_proxy_from_msg(msg) -> dict[str, float]:
    return {m: max(float(msg.auc[m]) - 0.5, 0.0) for m in MODS}


def train_window(model, opt, rows, device, image_tf, *, lr_mult):
    model.train()
    model.image_encoder.eval()
    base.set_lrs(opt, lr_mult)
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=BATCH, shuffle=True)
    crit = nn.CrossEntropyLoss()
    losses, it = [], iter(loader)
    t0 = time.perf_counter()
    for _ in range(STEPS):
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


def alpha_equal():
    return {m: 1.0 / len(MODS) for m in MODS}


def run_policy(stream, device, image_tf, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = base.AmazonAGOD().to(device)
    opt = base.build_optim(model)
    ref = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref, device, image_tf)

    ref_tr, ref_te = split_rows(ref, SEED - 1)
    b_tr, y_tr = base.extract_blocks(model, ref_tr, device, image_tf)
    b_te, y_te = base.extract_blocks(model, ref_te, device, image_tf)
    uni_ref = measure_uni_acc(b_tr, y_tr, b_te, y_te, MODS, seed=SEED)

    ema = alpha_equal()
    router = SmoothDriftNoiseRouter(MODS, SmoothRouterConfig(tau=TAU, ema=EMA, beta=BETA))
    traj = []

    for w in stream["windows"]:
        adapt, hold = split_rows(w["rows"], SEED + 13 * w["t"])
        cur_b, y_cur = base.extract_blocks(model, adapt, device, image_tf)
        hold_b, y_hold = base.extract_blocks(model, hold, device, image_tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        uni = measure_uni_acc(cur_b, y_cur, hold_b, y_hold, MODS, seed=SEED + 7 * w["t"])
        mmd = mmd_proxy_from_msg(msg)

        gate_info: dict = {}
        if policy == "equal":
            alpha = alpha_equal()
            lr_mult = weights_from_alpha(alpha, MODS, beta=BETA)
        elif policy == "msg_softmax":
            raw = msg.alpha["B3"]
            for m in MODS:
                ema[m] = EMA * ema[m] + (1.0 - EMA) * float(raw[m])
            s = sum(ema.values())
            alpha = {m: ema[m] / s for m in MODS}
            lr_mult = weights_from_alpha(alpha, MODS, beta=BETA)
        elif policy == "concept_minus_cov":
            cov = {
                m: max(float(msg.auc[m]) - 0.5, 0.0) * (1.0 + float(msg.vimp[m]))
                for m in MODS
            }
            con = {m: max(float(msg.po[m]), 0.0) for m in MODS}
            cov_n, con_n = normalize_nonneg(cov, MODS), normalize_nonneg(con, MODS)
            score = {m: con_n[m] - cov_n[m] for m in MODS}
            raw = softmax_tau(score, MODS, TAU)
            for m in MODS:
                ema[m] = EMA * ema[m] + (1.0 - EMA) * raw[m]
            s = sum(ema.values())
            alpha = {m: ema[m] / s for m in MODS}
            lr_mult = weights_from_alpha(alpha, MODS, beta=BETA)
        else:
            packed = router.update(
                po=msg.po,
                mmd=mmd,
                vimp=msg.vimp,
                uni_acc=uni,
                uni_acc_ref=uni_ref,
            )
            alpha = packed["alpha"]
            lr_mult = packed["lr_mult"]
            gate_info = {
                "frac_signal": packed["gate"]["frac_signal"],
                "reasons": packed["gate"]["reasons"],
                "is_signal": packed["gate"]["is_signal"],
                "adapter_only": packed["gate"]["adapter_only"],
                "g": packed["g"],
            }

        pre = eval_hold(model, hold, device, image_tf)
        loss, wall = train_window(model, opt, adapt, device, image_tf, lr_mult=lr_mult)
        post = eval_hold(model, hold, device, image_tf)

        keep = base.N_REF // 2
        ref_b = {
            m: np.vstack([ref_b[m][-keep:], cur_b[m][: min(keep, len(cur_b[m]))]])
            for m in MODS
        }
        y_ref = np.concatenate([y_ref[-keep:], y_cur[: min(keep, len(y_cur))]])

        row = {
            "t": int(w["t"]),
            "category": w.get("category"),
            "alpha": {m: float(alpha[m]) for m in MODS},
            "lr_mult": {k: float(v) for k, v in lr_mult.items()},
            "po": {m: float(msg.po[m]) for m in MODS},
            "vimp": {m: float(msg.vimp[m]) for m in MODS},
            "mmd_proxy": {m: float(mmd[m]) for m in MODS},
            "uni_acc": {m: float(uni[m]) for m in MODS},
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": post["acc"] - pre["acc"],
            "mse_pre": pre["mse"],
            "mse_post": post["mse"],
            "mse_drop": pre["mse"] - post["mse"],
            "train_loss": loss,
            "wall_ms": wall,
            "gate": gate_info,
        }
        traj.append(row)
        sig = gate_info.get("frac_signal", float("nan"))
        print(
            f"[amazon/{policy}] t={w['t']} "
            f"MSE {pre['mse']:.4f}→{post['mse']:.4f} (drop={row['mse_drop']:+.4f}) "
            f"Acc {pre['acc']:.3f}→{post['acc']:.3f} (d={row['acc_lift']:+.3f}) "
            f"frac_signal={sig}",
            flush=True,
        )
    return traj


def summarize(traj, policy):
    sigs = [
        r["gate"].get("frac_signal", np.nan) for r in traj if r.get("gate")
    ]
    return {
        "dataset": "amazon",
        "policy": policy,
        "n_windows": len(traj),
        "mean_mse_pre": float(np.mean([r["mse_pre"] for r in traj])),
        "mean_mse_post": float(np.mean([r["mse_post"] for r in traj])),
        "mean_mse_drop": float(np.mean([r["mse_drop"] for r in traj])),
        "frac_mse_drop_pos": float(np.mean([r["mse_drop"] > 0 for r in traj])),
        "mean_acc_lift": float(np.mean([r["acc_lift"] for r in traj])),
        "mean_acc_post": float(np.mean([r["acc_post"] for r in traj])),
        "mean_frac_signal": float(np.nanmean(sigs)) if sigs else float("nan"),
    }


def plot_board(cells, trajs, path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0), facecolor="#f7f5f1")
    names = [c["policy"] for c in cells]
    ax = axes[0]
    x = np.arange(len(names))
    ax.bar(x - 0.15, [c["mean_mse_pre"] for c in cells], 0.3, label="MSE pre", color="#9aa0a6")
    ax.bar(x + 0.15, [c["mean_mse_post"] for c in cells], 0.3, label="MSE post", color="#b85c38")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=18, ha="right")
    ax.set_ylabel("holdout MSE (Brier)")
    ax.set_title("Amazon: MSE pre → post")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    for pol, traj in trajs.items():
        ax.plot(
            [r["t"] for r in traj],
            [r["acc_lift"] for r in traj],
            marker="o",
            label=pol,
            lw=1.6,
        )
    ax.axhline(0.0, color="#999", ls=":", lw=0.9)
    ax.set_xlabel("window t")
    ax.set_ylabel("Acc lift")
    ax.set_title("Acc↑ by policy")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(cells, path: Path):
    rows = [
        f"| `{c['policy']}` | {c['mean_mse_pre']:.4f} | {c['mean_mse_post']:.4f} | "
        f"{c['mean_mse_drop']:+.4f} | {c['mean_acc_lift']:+.3f} | "
        f"{c['mean_acc_post']:.3f} | {c['mean_frac_signal']:.2f} |"
        for c in cells
    ]
    best_mse = max(cells, key=lambda c: c["mean_mse_drop"])
    best_acc = max(cells, key=lambda c: c["mean_acc_lift"])
    md = f"""# Smooth drift-vs-noise router (Amazon)

## Rule incorporated

1. **Drift vs noise gate** (per modality):
   - boost only if `PO` high **and** unimodal Acc not bad **and** VIMP/Fisher not low
   - else damp PO (noise path) / near adapter-only LR floor
2. **Smooth control law**:
   ```
   g = normalize(w1·PO_gated + w2·MMD + w3·VIMP)
   α = EMA(Softmax(g / τ))
   w = w0 · (β + (1-β)·α·|M|)
   ```

## Smoke

| policy | MSE pre | MSE post | MSE drop↑ | Acc↑ | Acc post | frac signal |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

Best MSE-drop: **`{best_mse['policy']}`** ({best_mse['mean_mse_drop']:+.4f}).
Best Acc↑: **`{best_acc['policy']}`** ({best_acc['mean_acc_lift']:+.3f}).

```bash
PYTHONPATH=. python3 scripts/run_agod_smooth_drift_noise.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(cells, path: Path):
    lines = [
        f"{c['policy'].replace('_', '\\_')} & {c['mean_mse_pre']:.4f} & "
        f"{c['mean_mse_post']:.4f} & {c['mean_mse_drop']:+.4f} & "
        f"{c['mean_acc_lift']:+.3f} \\\\"
        for c in cells
    ]
    tex = (
        "% Smooth drift-vs-noise router\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Amazon online: smooth drift-vs-noise gate + EMA routing "
        "vs equal / MSG / concept$-$cov (holdout Brier/MSE and Acc).}\n"
        "\\label{tab:agod-smooth-drift-noise}\n"
        "\\begin{tabular}{lrrr}\\toprule\n"
        "policy & MSE pre & MSE post & MSE drop & Acc$\\uparrow$ \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--policies", nargs="+", default=list(POLICIES), choices=list(POLICIES))
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    print("loading Amazon samples…", flush=True)
    shards = sorted(base.SHARD_DIR.glob("*.tar.gz"))[:3]
    samples = base.load_shards(shards)
    stream = base.make_stream(samples)
    image_tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize(
                mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]
            ),
        ]
    )
    device = args.device

    cells, trajs = [], {}
    for pol in args.policies:
        traj = run_policy(stream, device, image_tf, pol)
        trajs[pol] = traj
        cells.append(summarize(traj, pol))

    payload = {
        "agod_version": "0.1.0",
        "focus": "drift-vs-noise gate + smooth EMA routing on Amazon",
        "rule": {
            "gate": "PO high + uni Acc ok + VIMP/Fisher ok → boost; else damp/adapter-only",
            "control": "g=norm(w1 POg+w2 MMD+w3 VIMP); α=EMA(softmax(g/τ)); w=w0(β+(1-β)α|M|)",
        },
        "cells": cells,
        "trajectory": trajs,
    }
    (OUT / "agod_smooth_drift_noise.json").write_text(json.dumps(payload, indent=2))
    plot_board(cells, trajs, OUT / "AGOD_Smooth_Drift_Noise_Board.png")
    write_docs(cells, OUT / "README.md")
    write_latex(cells, OUT / "AGOD_smooth_drift_noise_tables_only.tex")
    write_docs(cells, DOCS / "AGOD_smooth_drift_noise.md")
    write_latex(cells, DOCS / "AGOD_smooth_drift_noise_tables_only.tex")
    shutil.copy2(
        OUT / "AGOD_Smooth_Drift_Noise_Board.png",
        ART / "AGOD_Smooth_Drift_Noise_Board.png",
    )
    shutil.copy2(OUT / "README.md", ART / "README.md")

    print("=== summary ===")
    for c in cells:
        print(
            f"  {c['policy']}: MSE drop={c['mean_mse_drop']:+.4f} "
            f"Acc↑={c['mean_acc_lift']:+.3f} Acc_post={c['mean_acc_post']:.3f} "
            f"frac_signal={c['mean_frac_signal']:.2f}"
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
