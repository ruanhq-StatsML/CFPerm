#!/usr/bin/env python3
"""Amazon smoke harness for AGOD MMD concept↑ / covariate↓ LR.

Thin CLI over ``agod`` controllers + ``run_agod_amazon_modality_lr`` IO/model.
Default Acc policy is B5 (MMD cov + PO concept). See
``docs/agod/AGOD_ml_infra_justification.md``.

  PYTHONPATH=. python3 scripts/run_agod_amazon_mmd_lr.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, roc_auc_score
from torch.utils.data import DataLoader
from torchvision import transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import agod  # noqa: E402
import run_agod_amazon_modality_lr as base  # noqa: E402
from agod.lr_controller import EMARouter, alpha_to_lr  # noqa: E402
from agod.policies import POLICIES, select_policy  # noqa: E402

OUT = ROOT / "results" / "agod_amazon"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_amazon")

SEED = base.SEED
MODS = base.MODS
BATCH = base.BATCH
STEPS = 32
HOLD = 0.40
EMA = 0.40
TAU = 0.30
BETA = 0.10
LR0_SCALE = 2.0
KAPPA = 1.25


def split_rows(rows, seed):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    n_h = max(24, min(len(rows) // 2, int(len(rows) * HOLD)))
    return [rows[i] for i in idx[n_h:]], [rows[i] for i in idx[:n_h]]


@torch.no_grad()
def eval_acc(model, rows, device, tf):
    model.eval()
    if len(rows) < 8:
        return {"acc": float("nan"), "auc": float("nan"), "n": len(rows)}
    loader = DataLoader(base.AmazonDS(rows, tf), batch_size=BATCH, shuffle=False)
    logits, ys = [], []
    for imgs, txts, y in loader:
        logits.append(
            model(imgs.to(device), txts, return_mods=False).cpu().numpy()
        )
        ys.append(y.numpy())
    L, Y = np.concatenate(logits), np.concatenate(ys)
    acc = float(accuracy_score(Y, L.argmax(1)))
    auc = float("nan")
    if len(np.unique(Y)) > 1:
        try:
            p = torch.softmax(torch.from_numpy(L), 1).numpy()[:, 1]
            auc = float(roc_auc_score(Y, p))
        except Exception:
            pass
    return {"acc": acc, "auc": auc, "n": int(len(Y))}


def train_window(model, opt, rows, device, tf, *, lr_mult):
    model.train()
    model.image_encoder.eval()
    base.set_lrs(opt, lr_mult)
    for g in opt.param_groups:
        g["lr"] *= LR0_SCALE
    loader = DataLoader(base.AmazonDS(rows, tf), batch_size=BATCH, shuffle=True)
    crit = nn.CrossEntropyLoss()
    losses, it = [], iter(loader)
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
    return float(np.mean(losses))


def run_policy(ctor, stream, device, tf, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = ctor().to(device)
    opt = base.build_optim(model)
    ref_b, y_ref = base.extract_blocks(model, stream["ref"], device, tf)
    router = EMARouter(MODS, ema=EMA)
    traj = []

    for w in stream["windows"]:
        adapt, hold = split_rows(w["rows"], SEED + 13 * w["t"])
        cur_b, y_cur = base.extract_blocks(model, adapt, device, tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        raw, decomp, gain = select_policy(
            policy,
            msg=msg,
            blocks0=ref_b,
            blocks1=cur_b,
            y0=y_ref,
            y1=y_cur,
            mods=MODS,
            tau=TAU,
            kappa=KAPPA,
            seed=SEED + 10 * w["t"],
        )
        alpha = router.update(raw)
        lr_mult = alpha_to_lr(alpha, MODS, beta=BETA, gain=gain)

        pre = eval_acc(model, hold, device, tf)
        loss = train_window(model, opt, adapt, device, tf, lr_mult=lr_mult)
        post = eval_acc(model, hold, device, tf)
        dacc = post["acc"] - pre["acc"]
        dauc = (
            post["auc"] - pre["auc"]
            if np.isfinite(pre["auc"]) and np.isfinite(post["auc"])
            else float("nan")
        )

        ref_b = {
            m: np.vstack(
                [ref_b[m][-(base.N_REF // 2) :], cur_b[m][: base.N_REF // 2]]
            )
            for m in MODS
        }
        y_ref = np.concatenate(
            [y_ref[-(base.N_REF // 2) :], y_cur[: base.N_REF // 2]]
        )

        traj.append(
            {
                "t": w["t"],
                "category": w["category"],
                "cov": decomp["cov"],
                "concept": decomp["con"],
                "score": decomp["score"],
                "alpha": alpha,
                "gain": gain,
                "lr_mult": lr_mult,
                "acc_pre": pre["acc"],
                "acc_post": post["acc"],
                "acc_lift": dacc,
                "auc_lift": dauc,
                "train_loss": loss,
            }
        )
        print(
            f"[{policy}] t={w['t']} {w['category'][:22]:<22} "
            f"con={ {m: round(decomp['con'][m], 2) for m in MODS} } "
            f"cov={ {m: round(decomp['cov'][m], 2) for m in MODS} } "
            f"LR×={ {m: round(lr_mult[m], 2) for m in MODS} } | "
            f"acc {pre['acc']:.3f}→{post['acc']:.3f} (Δ={dacc:+.3f})",
            flush=True,
        )
    return traj


def summarize(results):
    out = {}
    for pol, traj in results.items():
        lifts = [r["acc_lift"] for r in traj]
        posts = [r["acc_post"] for r in traj]
        out[pol] = {
            "mean_acc_lift": float(np.nanmean(lifts)),
            "mean_acc_post": float(np.nanmean(posts)),
            "mean_auc_lift": float(np.nanmean([r["auc_lift"] for r in traj])),
            "wins_vs_zero": int(sum(1 for x in lifts if x > 0)),
            "n_windows": len(traj),
        }
    for a, b in [
        ("B3", "B1"),
        ("B3", "B2"),
        ("B5", "B1"),
        ("B5", "B2"),
        ("B5", "B3"),
        ("B5", "B4"),
    ]:
        out[f"{a}_minus_{b}_lift"] = out[a]["mean_acc_lift"] - out[b]["mean_acc_lift"]
        out[f"{a}_minus_{b}_post"] = out[a]["mean_acc_post"] - out[b]["mean_acc_post"]
    return out


def plot_dash(results, summary, path: Path):
    fig = plt.figure(figsize=(14.0, 9.0), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.28)
    fig.suptitle(
        "Amazon AGOD (system): MMD cov↓ + PO concept↑ — Acc",
        fontsize=14,
        fontweight="bold",
    )
    b5 = results["B5"]
    ts = [r["t"] + 1 for r in b5]

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(ts, [r["concept"]["image"] for r in b5], "o-", color="#C53030", label="concept_image")
    ax.plot(ts, [r["cov"]["image"] for r in b5], "s--", color="#2B6CB0", label="cov_image")
    ax.plot(ts, [r["concept"]["text"] for r in b5], "o-", color="#DD6B20", label="concept_text")
    ax.plot(ts, [r["cov"]["text"] for r in b5], "s--", color="#38A169", label="cov_text")
    ax.set_title("B5 sensors (PO concept / MMD cov)")
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(ts, [r["lr_mult"]["image"] for r in b5], "o-", color="#C05621", label="LR× image")
    ax.plot(ts, [r["lr_mult"]["text"] for r in b5], "o-", color="#2B6CB0", label="LR× text")
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("B5 actuator (per-modality LR)")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 0])
    colors = ["#718096", "#DD6B20", "#2B6CB0", "#C53030", "#805AD5"]
    marks = ["o", "s", "^", "D", "P"]
    for pol, c, mk in zip(POLICIES, colors, marks):
        ax.plot(
            [r["t"] + 1 for r in results[pol]],
            [r["acc_lift"] for r in results[pol]],
            marker=mk,
            color=c,
            label=pol,
        )
    ax.axhline(0.0, color="#999", ls="--", lw=0.8)
    ax.set_title("Held-out Acc lift")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    xs = np.arange(len(POLICIES))
    lifts = [summary[p]["mean_acc_lift"] for p in POLICIES]
    posts = [summary[p]["mean_acc_post"] for p in POLICIES]
    ax.bar(xs - 0.15, lifts, 0.3, color="#C53030", label="mean Acc lift")
    ax.bar(xs + 0.15, posts, 0.3, color="#4A5568", label="mean Acc post")
    ax.set_xticks(xs)
    ax.set_xticklabels(list(POLICIES))
    ax.axhline(0.0, color="#999", ls="--", lw=0.7)
    ax.set_title(
        f"B5−B1={summary['B5_minus_B1_lift']:+.3f}  "
        f"B5−B2={summary['B5_minus_B2_lift']:+.3f}"
    )
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.012,
        f"agod {agod.__version__}: sensors in agod.mmd/shift, actuator in agod.lr_controller; "
        "Amazon harness owns IO/model only",
        ha="center",
        fontsize=8.5,
        color="#333",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device} agod={agod.__version__}", flush=True)
    shards = sorted(base.SHARD_DIR.glob("*.tar.gz"))
    samples = base.load_shards(shards)
    print(f"samples={len(samples)}", flush=True)
    tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
        ]
    )
    stream = base.make_stream(samples)
    print(
        f"ref={stream['ref_cat']} wins={[w['category'] for w in stream['windows']]}",
        flush=True,
    )

    results = {}
    for pol in POLICIES:
        print(f"\n===== {pol} =====", flush=True)
        results[pol] = run_policy(base.AmazonAGOD, stream, device, tf, pol)

    summary = summarize(results)
    print("\n=== ACCURACY SUMMARY ===", flush=True)
    print(json.dumps(summary, indent=2), flush=True)

    payload = {
        "agod_version": agod.__version__,
        "default_policy": "B5",
        "rule": "B5: MMD² cov + PO concept → Softmax → LR",
        "ml_infra": "docs/agod/AGOD_ml_infra_justification.md",
        "summary": summary,
        "trajectory": results,
    }
    jp = OUT / "agod_amazon_mmd_lr.json"
    jp.write_text(json.dumps(payload, indent=2, default=float))
    dash = OUT / "AGOD_Amazon_MMD_LR_Acc_Dashboard.png"
    plot_dash(results, summary, dash)

    infra = DOCS / "AGOD_ml_infra_justification.md"
    assert infra.exists(), infra
    for p in [dash, jp, infra]:
        (ART / Path(p).name).write_bytes(Path(p).read_bytes())
    print("dashboard", dash, flush=True)
    print(
        "B5-B1",
        summary["B5_minus_B1_lift"],
        "B5-B2",
        summary["B5_minus_B2_lift"],
        flush=True,
    )


if __name__ == "__main__":
    main()
