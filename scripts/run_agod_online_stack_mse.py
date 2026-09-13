#!/usr/bin/env python3
"""Online stacking + erank/modality-balance — track holdout MSE (Brier).

Variants on the MSR-VTT online-window stream:

  mean_ce         — mean-pool fusion + CE
  stack_ce        — online stacking over modality logits + CE
  stack_erank_bal — stacking + CE + erank↑ / align↓ / weight-balance

Primary metric: holdout MSE = mean((softmax(z)−onehot(y))²) pre→post.

  PYTHONPATH=. python3 scripts/run_agod_online_stack_mse.py
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
from sklearn.metrics import accuracy_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import run_agod_gradcos_lr as g
from agod.lr_controller import EMARouter
from agod.online_stack import (
    MeanFusion,
    StackFusion,
    erank_balance_loss,
    gram_erank_torch,
    probs_mse,
    probs_mse_torch,
)

OUT = ROOT / "results" / "agod_online_stack_mse"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_online_stack_mse")

SEED = g.SEED
MODS = list(g.MSRVTT_MODS)
VARIANTS = ("mean_ce", "stack_ce", "stack_erank_bal")
LAM_ERANK, LAM_ALIGN, LAM_BAL = 0.35, 0.35, 0.15
ERANK_FLOOR = 0.85 * len(MODS)


def to_batch(feats, y, idx, mods, device):
    batch = {
        m: torch.tensor(feats[m][idx], dtype=torch.float32, device=device) for m in mods
    }
    yy = torch.tensor(y[idx], dtype=torch.long, device=device)
    return batch, yy


@torch.no_grad()
def eval_hold(model, feats, y, idx, mods, device):
    model.eval()
    if len(idx) < 8:
        return {"acc": float("nan"), "mse": float("nan"), "n": len(idx)}
    logits, ys = [], []
    for s in range(0, len(idx), g.BATCH):
        sl = idx[s : s + g.BATCH]
        batch, yy = to_batch(feats, y, sl, mods, device)
        logits.append(model(batch).cpu().numpy())
        ys.append(yy.cpu().numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    return {
        "acc": float(accuracy_score(Y, L.argmax(1))),
        "mse": float(probs_mse(L, Y, n_class=L.shape[1])),
        "n": int(len(Y)),
    }


def _nanmean(xs):
    arr = np.asarray(xs, float)
    if arr.size == 0 or not np.isfinite(arr).any():
        return float("nan")
    return float(np.nanmean(arr))


def train_window(
    model,
    opt,
    feats,
    y,
    idx,
    mods,
    device,
    *,
    kind: str,
    steps: int,
    alpha: dict | None = None,
):
    model.train()
    crit = nn.CrossEntropyLoss()
    losses, mses, eranks = [], [], []
    rng = np.random.default_rng(SEED + int(np.asarray(idx).sum()) % 100000)
    t0 = time.perf_counter()
    for _ in range(steps):
        sel = rng.choice(idx, size=min(g.BATCH, len(idx)), replace=len(idx) < g.BATCH)
        batch, yy = to_batch(feats, y, sel, mods, device)
        if kind == "mean_ce":
            logits, _hs = model(batch, return_h=True)
            loss = crit(logits, yy)
            erank_v = float("nan")
        else:
            logits, hs, _lm, stack_w = model(batch, return_parts=True)
            loss = crit(logits, yy)
            erank_v = float("nan")
            if kind == "stack_erank_bal":
                bal = erank_balance_loss(
                    hs,
                    mods,
                    stack_w=stack_w,
                    alpha=alpha,
                    erank_floor=ERANK_FLOOR,
                    lambda_align=LAM_ALIGN,
                    lambda_erank=LAM_ERANK,
                    lambda_bal=LAM_BAL,
                )
                loss = loss + bal["loss"]
                erank_v = float(bal["erank"].detach().cpu())
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        losses.append(float(loss.item()))
        with torch.no_grad():
            mses.append(float(probs_mse_torch(logits.detach(), yy).cpu()))
        eranks.append(erank_v)
    return {
        "train_loss": float(np.mean(losses)),
        "train_mse": float(np.mean(mses)),
        "mean_erank_aux": _nanmean(eranks),
        "wall_ms": (time.perf_counter() - t0) * 1000.0,
    }


def select_alpha(feats, y, ref_idx, adapt_idx, mods, *, seed):
    b0 = {m: feats[m][ref_idx] for m in mods}
    b1 = {m: feats[m][adapt_idx] for m in mods}
    y0 = y[ref_idx].astype(float)
    y1 = y[adapt_idx].astype(float)
    try:
        msg = g.domain_msg(b0, b1, y0, y1, mods, seed=seed)
        raw, _ = g.select_b5_raw(msg, b0, b1, y0, y1, mods, seed=seed)
        return raw
    except Exception:
        return {m: 1.0 / len(mods) for m in mods}


def run_variant(device: str, kind: str):
    feats, y = g.load_msrvtt()
    mods = list(MODS)
    stream = g.make_stream_y(y, g.N_CUR_MSRVTT)
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    dims = {m: feats[m].shape[1] for m in mods}
    if kind == "mean_ce":
        model = MeanFusion(dims, mods, fuse=g.FUSE).to(device)
    else:
        model = StackFusion(dims, mods, fuse=g.FUSE).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=g.LR0)
    router = EMARouter(mods, ema=g.EMA)
    ref_idx = stream["ref_idx"].copy()
    warm, _ = g.split_hold(ref_idx, SEED)
    train_window(
        model, opt, feats, y, warm, mods, device, kind=kind, steps=g.STEPS_M, alpha=None
    )

    traj = []
    for win in stream["windows"]:
        adapt, hold = g.split_hold(win["idx"], SEED + 13 * win["t"])
        raw = select_alpha(feats, y, ref_idx, adapt, mods, seed=SEED + 10 * win["t"])
        alpha = router.update(raw)

        model.eval()
        with torch.no_grad():
            sel = np.random.default_rng(SEED + 99).choice(
                adapt, size=min(g.BATCH, len(adapt)), replace=len(adapt) < g.BATCH
            )
            batch, _yy = to_batch(feats, y, sel, mods, device)
            if kind == "mean_ce":
                _logits, hs = model(batch, return_h=True)
                w_np = {m: 1.0 / len(mods) for m in mods}
            else:
                _logits, hs, _lm, stack_w = model(batch, return_parts=True)
                w_np = {
                    m: float(stack_w.detach()[i].cpu()) for i, m in enumerate(mods)
                }
            try:
                protos = [hs[m].mean(0) for m in mods]
                erank_h, align_h = gram_erank_torch(protos)
                erank_log = float(erank_h.cpu())
                align_log = float(align_h.cpu())
            except Exception:
                erank_log, align_log = float("nan"), float("nan")

        pre = eval_hold(model, feats, y, hold, mods, device)
        tr = train_window(
            model,
            opt,
            feats,
            y,
            adapt,
            mods,
            device,
            kind=kind,
            steps=g.STEPS_M,
            alpha=alpha,
        )
        post = eval_hold(model, feats, y, hold, mods, device)

        keep = g.N_REF // 2
        ref_idx = np.concatenate([ref_idx[-keep:], adapt[: min(keep, len(adapt))]])
        row = {
            "t": int(win["t"]),
            "alpha": {m: float(alpha[m]) for m in mods},
            "stack_w": w_np,
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": post["acc"] - pre["acc"],
            "mse_pre": pre["mse"],
            "mse_post": post["mse"],
            "mse_drop": pre["mse"] - post["mse"],
            "train_loss": tr["train_loss"],
            "train_mse": tr["train_mse"],
            "mean_erank_aux": tr["mean_erank_aux"],
            "erank_hidden": erank_log,
            "align_hidden": align_log,
            "wall_ms": tr["wall_ms"],
        }
        traj.append(row)
        print(
            f"[msrvtt/{kind}] t={win['t']} "
            f"MSE {pre['mse']:.4f}→{post['mse']:.4f} (drop={row['mse_drop']:+.4f}) "
            f"Acc {pre['acc']:.3f}→{post['acc']:.3f} (d={row['acc_lift']:+.3f}) "
            f"erank_aux={tr['mean_erank_aux']}",
            flush=True,
        )
    return traj


def summarize(traj, kind):
    return {
        "dataset": "msrvtt",
        "variant": kind,
        "n_windows": len(traj),
        "mean_mse_pre": float(np.mean([r["mse_pre"] for r in traj])),
        "mean_mse_post": float(np.mean([r["mse_post"] for r in traj])),
        "mean_mse_drop": float(np.mean([r["mse_drop"] for r in traj])),
        "frac_mse_drop_pos": float(np.mean([r["mse_drop"] > 0 for r in traj])),
        "mean_acc_lift": float(np.mean([r["acc_lift"] for r in traj])),
        "mean_train_mse": float(np.mean([r["train_mse"] for r in traj])),
        "mean_erank_aux": _nanmean([r["mean_erank_aux"] for r in traj]),
    }


def plot_board(cells, trajs, path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), facecolor="#f7f5f1")
    names = [c["variant"] for c in cells]
    ax = axes[0]
    x = np.arange(len(names))
    ax.bar(x - 0.15, [c["mean_mse_pre"] for c in cells], 0.3, label="MSE pre", color="#9aa0a6")
    ax.bar(x + 0.15, [c["mean_mse_post"] for c in cells], 0.3, label="MSE post", color="#b85c38")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=12)
    ax.set_ylabel("holdout MSE (Brier)")
    ax.set_title("Online window: MSE pre → post")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    for kind, traj in trajs.items():
        ax.plot(
            [r["t"] for r in traj],
            [r["mse_drop"] for r in traj],
            marker="o",
            label=kind,
            lw=1.8,
        )
    ax.axhline(0.0, color="#999", ls=":", lw=0.9)
    ax.set_xlabel("window t")
    ax.set_ylabel("MSE drop (pre−post)")
    ax.set_title("Does MSE fall within the window?")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(cells, path: Path):
    rows = [
        f"| `{c['variant']}` | {c['mean_mse_pre']:.4f} | {c['mean_mse_post']:.4f} | "
        f"{c['mean_mse_drop']:+.4f} | {c['frac_mse_drop_pos']:.2f} | "
        f"{c['mean_acc_lift']:+.3f} | {c['mean_erank_aux']:.2f} |"
        for c in cells
    ]
    best = max(cells, key=lambda c: c["mean_mse_drop"])
    md = f"""# Online stacking + erank/modality balance — MSE (MSR-VTT)

## Setup

1. **Online stacking**: modality logits fused by learned `w=softmax(ψ)`
2. **Erank / modality balance** (`stack_erank_bal` only):
   - `L_erank = ReLU(erank_floor − erank(G_h))`
   - `L_align = mean max(cos_ij, 0)²`
   - `L_bal = KL(w ‖ prior)` (prior=α if erank collapsed else uniform)
3. Task = CE; **MSE = Brier** on holdout (`pre→post` drop = good)

## Smoke (honest)

| variant | MSE pre | MSE post | MSE drop↑ | frac drop>0 | Acc↑ | mean erank aux |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

Best MSE-drop: **`{best['variant']}`** ({best['mean_mse_drop']:+.4f}).

On low-corr MSR-VTT, hidden erank stays ≈|M| so `L_erank` barely fires;
mean MSE drop is **not** positive overall — same regime lesson as soft_decorr.
Amazon high-corr is the next place to expect MSE↓.

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stack_mse.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(cells, path: Path):
    lines = [
        f"{c['variant'].replace('_', '\\_')} & {c['mean_mse_pre']:.4f} & "
        f"{c['mean_mse_post']:.4f} & {c['mean_mse_drop']:+.4f} & "
        f"{c['mean_acc_lift']:+.3f} \\\\"
        for c in cells
    ]
    tex = (
        "% Online stacking + erank balance — MSE\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{MSR-VTT online windows: holdout Brier/MSE under mean fusion, "
        "online stacking, and stacking+erank/modality balance.}\n"
        "\\label{tab:agod-online-stack-mse}\n"
        "\\begin{tabular}{lrrr}\\toprule\n"
        "variant & MSE pre & MSE post & MSE drop & Acc$\\uparrow$ \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=list(VARIANTS))
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    cells, trajs = [], {}
    for kind in args.variants:
        traj = run_variant(args.device, kind)
        trajs[kind] = traj
        cells.append(summarize(traj, kind))

    payload = {
        "agod_version": "0.1.0",
        "focus": "online stacking + erank/modality balance; holdout MSE drop",
        "loss": {
            "task": "CE(stacked_or_mean_logits, y)",
            "aux": "L_erank + L_align + L_bal (stack_erank_bal only)",
            "mse_metric": "mean((softmax(z)-onehot(y))^2) on holdout",
        },
        "cells": cells,
        "trajectory": trajs,
    }
    (OUT / "agod_online_stack_mse.json").write_text(json.dumps(payload, indent=2))
    plot_board(cells, trajs, OUT / "AGOD_Online_Stack_MSE_Board.png")
    write_docs(cells, OUT / "README.md")
    write_latex(cells, OUT / "AGOD_online_stack_mse_tables_only.tex")
    write_docs(cells, DOCS / "AGOD_online_stack_mse.md")
    write_latex(cells, DOCS / "AGOD_online_stack_mse_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Online_Stack_MSE_Board.png", ART / "AGOD_Online_Stack_MSE_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")

    print("=== MSE summary ===")
    for c in cells:
        print(
            f"  {c['variant']}: MSE {c['mean_mse_pre']:.4f}→{c['mean_mse_post']:.4f} "
            f"(drop={c['mean_mse_drop']:+.4f}, frac>0={c['frac_mse_drop_pos']:.2f}) "
            f"Acc↑={c['mean_acc_lift']:+.3f}"
        )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
