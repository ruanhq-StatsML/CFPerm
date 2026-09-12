#!/usr/bin/env python3
"""Image/Text AGOD portable prototype (same control plane as Amazon / MSR-VTT).

Datasets under data/img_txt/<name>/ with packed img_feats + txt_feats (+ optional bbox).
Policies: B1 equal | B2 RF con-cov | B5 MMD-cov+PO | B5g α-gate | B5r random-shutdown.

Hard-gate semantics (important):
  α_m < θ → skip modality proj *BWD / adapt* (structured adapt-dropout).
  FWD inference still runs. This is NOT "OOD large ⇒ skip inference".

  PYTHONPATH=. python3 scripts/run_agod_imgtxt_mmd_lr.py
  PYTHONPATH=. python3 scripts/run_agod_imgtxt_mmd_lr.py --datasets coco_outdoor_indoor
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import agod
from agod.lr_controller import EMARouter, alpha_to_lr, softmax_scores
from agod.shift import decompose_hybrid, decompose_mmd, decompose_rf, residual_concept

DATA_ROOT = ROOT / "data" / "img_txt"
OUT = ROOT / "results" / "agod_imgtxt"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_imgtxt")

SEED = 2026
N_REF = 400
N_CUR = 320
T_WIN = 6
HOLD = 0.35
STEPS = 35
BATCH = 64
LR0 = 3e-3
EMA = 0.40
TAU = 0.30
BETA = 0.10
GATE_TH = 0.30
FUSE = 128
C_PF, C_PB, C_SF, C_SB = 1.0, 2.0, 0.5, 1.0

DEFAULT_DATASETS = (
    "coco_outdoor_indoor",
    "coco_time_order",
    "coco_center_split",
    "mm_train_test",
    "fashion_iq",
    "indiana_cxr",
    "microscopy_clip",
)
POLICIES = ("B1", "B2", "B5", "B5g", "B5r")
P_KEEP_RAND = 0.67  # random adapt-shutdown keep prob; FWD always on


def load_dataset(name: str):
    d = DATA_ROOT / name
    img = np.load(d / "img_feats.npy").astype(np.float32)
    txt = np.load(d / "txt_feats.npy").astype(np.float32)
    if (d / "labels.npy").exists():
        y_raw = np.load(d / "labels.npy").astype(np.int64)
    else:
        import pandas as pd
        meta = pd.read_csv(d / "df_metadata.csv")
        if "label" not in meta.columns:
            raise FileNotFoundError(f"{name}: need labels.npy or df_metadata.label")
        y_raw = meta["label"].to_numpy().astype(np.int64)
    mode = int(np.bincount(y_raw).argmax())
    y = (y_raw == mode).astype(np.int64)
    feats = {"image": img, "text": txt}
    if (d / "bbox_named_feats.npy").exists():
        feats["bbox"] = np.load(d / "bbox_named_feats.npy").astype(np.float32)
    return feats, y, {"name": name, "n": int(len(y)), "pos_rate": float(y.mean()), "mode_label": mode}


def make_stream(y):
    rng = np.random.default_rng(SEED)
    idx0 = np.where(y == 0)[0].copy()
    idx1 = np.where(y == 1)[0].copy()
    rng.shuffle(idx0)
    rng.shuffle(idx1)
    n0, n1 = N_REF // 2, N_REF - N_REF // 2
    n0, n1 = min(n0, len(idx0)), min(n1, len(idx1))
    ref_idx = np.concatenate([idx0[:n0], idx1[:n1]])
    rng.shuffle(ref_idx)
    rest0, rest1 = idx0[n0:], idx1[n1:]

    def take(pool, k, t):
        if len(pool) == 0 or k <= 0:
            return np.array([], dtype=int)
        start = (t * max(k, 1)) % len(pool)
        return np.concatenate([pool[start:], pool[:start]])[:k]

    windows = []
    for t in range(T_WIN):
        p1 = 0.25 + 0.50 * (t / max(T_WIN - 1, 1))
        k1 = int(N_CUR * p1)
        k0 = N_CUR - k1
        cur = np.concatenate([take(rest0, k0, t), take(rest1, k1, t)])
        rng.shuffle(cur)
        windows.append(
            {
                "t": t,
                "p1": float(p1),
                "idx": cur,
                "n0": int((y[cur] == 0).sum()),
                "n1": int((y[cur] == 1).sum()),
            }
        )
    return {"ref_idx": ref_idx, "windows": windows}


def blocks_from_idx(feats, y, idx, mods):
    return {m: feats[m][idx] for m in mods}, y[idx].astype(float)


def domain_msg(b0, b1, y0, y1, mods, *, seed: int):
    auc, vimp, po = {}, {}, {}
    for i, m in enumerate(mods):
        X0, X1 = b0[m], b1[m]
        X = np.vstack([X0, X1])
        W = np.array([0] * len(X0) + [1] * len(X1))
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
        po[m] = residual_concept(X0, y0, X1, y1, seed=seed + 31 + i)
    return SimpleNamespace(auc=auc, vimp=vimp, po=po)


def flops_rel(active: dict, mods) -> float:
    fwd = len(mods) * C_PF + C_SF
    bwd = sum(C_PB for m in mods if active[m]) + C_SB
    full = len(mods) * C_PF + C_SF + len(mods) * C_PB + C_SB
    return float((fwd + bwd) / full)


class ImgTxtFusion(nn.Module):
    def __init__(self, dims: dict[str, int], mods: list[str]):
        super().__init__()
        self.mods = list(mods)
        self.projs = nn.ModuleDict(
            {m: nn.Sequential(nn.Linear(dims[m], FUSE), nn.ReLU(), nn.Dropout(0.1)) for m in mods}
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


def build_optim(model: ImgTxtFusion, mods):
    g = model.param_groups()
    return torch.optim.AdamW(
        [{"params": g[m], "lr": LR0, "name": m} for m in mods]
        + [{"params": g["shared"], "lr": LR0, "name": "shared"}]
    )


def set_lrs(opt, mult: dict):
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
        return {"acc": float("nan"), "auc": float("nan"), "n": len(idx)}
    logits, ys = [], []
    for s in range(0, len(idx), BATCH):
        sl = idx[s : s + BATCH]
        batch, yy = to_batch(feats, y, sl, mods, device)
        logits.append(model(batch).cpu().numpy())
        ys.append(yy.cpu().numpy())
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


def select_raw(policy, msg, b0, b1, y0, y1, mods, *, seed):
    unit = {m: 1.0 for m in mods}
    if policy == "B1":
        decomp = {"cov": unit, "con": unit, "score": {m: 0.0 for m in mods}}
        return {m: 1.0 / len(mods) for m in mods}, decomp
    if policy == "B2":
        decomp = decompose_rf(msg, mods)
        return softmax_scores(decomp["score"], mods, TAU), decomp
    mmd_d = decompose_mmd(b0, b1, y0, y1, mods, seed=seed)
    hy = decompose_hybrid(msg, mmd_d, mods)
    return softmax_scores(hy["score"], mods, TAU), hy


def run_policy(feats, y, stream, mods, device, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    dims = {m: feats[m].shape[1] for m in mods}
    model = ImgTxtFusion(dims, mods).to(device)
    opt = build_optim(model, mods)
    router = EMARouter(mods, ema=EMA)
    ref_idx = stream["ref_idx"].copy()
    ref_b, y_ref = blocks_from_idx(feats, y, ref_idx, mods)

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
        cur_b, y_cur = blocks_from_idx(feats, y, adapt, mods)
        msg = domain_msg(ref_b, cur_b, y_ref, y_cur, mods, seed=SEED + 10 * w["t"])
        raw, decomp = select_raw(
            policy, msg, ref_b, cur_b, y_ref, y_cur, mods, seed=SEED + 10 * w["t"]
        )
        alpha = router.update(raw)
        lr_mult = alpha_to_lr(alpha, mods, beta=BETA)

        if policy == "B5g":
            active = {m: alpha[m] >= GATE_TH for m in mods}
            if not any(active.values()):
                mstar = max(mods, key=lambda m: alpha[m])
                active = {m: m == mstar for m in mods}
        elif policy == "B5r":
            rng_g = np.random.default_rng(SEED + 101 * int(w["t"]) + 17)
            active = {m: bool(rng_g.random() < P_KEEP_RAND) for m in mods}
            if not any(active.values()):
                active[mods[int(rng_g.integers(0, len(mods)))]] = True
        else:
            active = {m: True for m in mods}

        pre = eval_acc(model, feats, y, hold, mods, device)
        loss, wall = train_window(
            model, opt, feats, y, adapt, mods, device, lr_mult=lr_mult, active=active
        )
        post = eval_acc(model, feats, y, hold, mods, device)
        dacc = post["acc"] - pre["acc"]
        fr = flops_rel(active, mods)
        util = dacc / fr if fr > 0 else float("nan")

        keep = N_REF // 2
        ref_idx = np.concatenate([ref_idx[-keep:], adapt[: min(keep, len(adapt))]])
        ref_b, y_ref = blocks_from_idx(feats, y, ref_idx, mods)

        row = {
            "t": w["t"],
            "p1": w["p1"],
            "n0": w["n0"],
            "n1": w["n1"],
            "cov": decomp["cov"],
            "concept": decomp["con"],
            "score": decomp["score"],
            "alpha": alpha,
            "lr_mult": lr_mult,
            "selected": active,
            "flops_rel": fr,
            "wall_ms": wall,
            "acc_pre": pre["acc"],
            "acc_post": post["acc"],
            "acc_lift": dacc,
            "cost_utility": util,
            "train_loss": loss,
        }
        traj.append(row)
        print(
            f"[{policy}] t={w['t']} p1={w['p1']:.2f} "
            f"a={ {m: round(alpha[m], 2) for m in mods} } "
            f"sel={ {m: int(active[m]) for m in mods} } "
            f"FLOPs={fr:.2f} | acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj


def summarize(results, mods):
    out = {}
    for pol, traj in results.items():
        lifts = [r["acc_lift"] for r in traj]
        out[pol] = {
            "mean_acc_lift": float(np.nanmean(lifts)),
            "mean_acc_post": float(np.nanmean([r["acc_post"] for r in traj])),
            "mean_flops_rel": float(np.mean([r["flops_rel"] for r in traj])),
            "mean_wall_ms": float(np.mean([r["wall_ms"] for r in traj])),
            "mean_cost_utility": float(np.nanmean([r["cost_utility"] for r in traj])),
            "wins_vs_zero": int(sum(1 for x in lifts if x > 0)),
            "n_windows": len(traj),
            **{f"mean_lr_{m}": float(np.mean([r["lr_mult"][m] for r in traj])) for m in mods},
            **{f"frac_{m}_sel": float(np.mean([r["selected"][m] for r in traj])) for m in mods},
        }
    for a, b in [
        ("B5", "B1"),
        ("B5", "B2"),
        ("B5g", "B1"),
        ("B5g", "B5"),
        ("B5r", "B1"),
        ("B5r", "B5"),
        ("B5g", "B5r"),
    ]:
        out[f"{a}_minus_{b}_lift"] = out[a]["mean_acc_lift"] - out[b]["mean_acc_lift"]
        out[f"{a}_minus_{b}_flops"] = out[a]["mean_flops_rel"] - out[b]["mean_flops_rel"]
    return out


def plot_dataset(name, results, summary, mods, path: Path):
    pols = list(POLICIES)
    colors = ["#C53030", "#2B6CB0", "#38A169", "#805AD5", "#DD6B20"]
    fig = plt.figure(figsize=(14.0, 9.0), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.28)
    fig.suptitle(
        f"AGOD img/txt prototype · {name} (gate = adapt-dropout, not skip-infer)",
        fontsize=13,
        fontweight="bold",
    )
    b5 = results["B5"]
    ts = [r["t"] + 1 for r in b5]

    ax = fig.add_subplot(gs[0, 0])
    for m, c in zip(mods, colors):
        ax.plot(ts, [r["concept"][m] for r in b5], "o-", color=c, label=f"con_{m}")
        ax.plot(ts, [r["cov"][m] for r in b5], "s--", color=c, alpha=0.45, label=f"cov_{m}")
    ax.set_title("B5 sensors")
    ax.legend(frameon=False, fontsize=6, ncol=2)

    ax = fig.add_subplot(gs[0, 1])
    for m, c in zip(mods, colors):
        ax.plot(ts, [r["lr_mult"][m] for r in b5], "o-", color=c, label=f"LR* {m}")
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("B5 actuator LR")
    ax.legend(frameon=False, fontsize=7)

    ax = fig.add_subplot(gs[1, 0])
    for pol, c, mk in zip(pols, ["#718096", "#DD6B20", "#805AD5", "#C53030", "#2B6CB0"], ["o", "s", "D", "P", "^"]):
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
    xs = np.arange(len(pols))
    lifts = [summary[p]["mean_acc_lift"] for p in pols]
    flops = [summary[p]["mean_flops_rel"] for p in pols]
    ax.bar(xs - 0.15, lifts, 0.3, color="#C53030", label="mean Acc lift")
    ax.bar(xs + 0.15, flops, 0.3, color="#4A5568", label="mean flops_rel")
    ax.set_xticks(xs)
    ax.set_xticklabels(pols)
    ax.axhline(0.0, color="#999", ls="--", lw=0.7)
    ax.set_title(
        f"B5-B1={summary['B5_minus_B1_lift']:+.3f}  "
        f"B5g flops={summary['B5g']['mean_flops_rel']:.2f}"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.text(
        0.5,
        0.012,
        "hard-gate skips proj BWD when α low (often high cov); FWD still paid",
        ha="center",
        fontsize=8.5,
        color="#333",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_board(all_summary: dict, path: Path):
    names = list(all_summary.keys())
    pols = list(POLICIES)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), facecolor="#f7f5f1")
    fig.suptitle(
        "AGOD img/txt multi-dataset smoke · Acc lift vs adapt FLOPs",
        fontsize=13,
        fontweight="bold",
    )
    x = np.arange(len(names))
    w = 0.15
    for i, pol in enumerate(pols):
        axes[0].bar(
            x + (i - 2.0) * w,
            [all_summary[n][pol]["mean_acc_lift"] for n in names],
            w,
            label=pol,
        )
        axes[1].bar(
            x + (i - 2.0) * w,
            [all_summary[n][pol]["mean_flops_rel"] for n in names],
            w,
            label=pol,
        )
    axes[0].axhline(0, color="#999", ls="--", lw=0.8)
    axes[0].set_title("mean Acc lift")
    axes[1].set_title("mean flops_rel (adapt units)")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=15, ha="right")
        ax.legend(frameon=False, fontsize=7)
    fig.text(
        0.5,
        0.02,
        "B5g=α-gate · B5r=random-shutdown · FWD always on",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.92])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def run_one(name: str, device):
    feats, y, meta = load_dataset(name)
    mods = list(feats.keys())
    print(
        f"\n######## {name} n={meta['n']} mods={mods} pos={meta['pos_rate']:.3f} ########",
        flush=True,
    )
    results = {}
    for pol in POLICIES:
        print(f"\n===== {name}/{pol} =====", flush=True)
        results[pol] = run_policy(feats, y, make_stream(y), mods, device, pol)

    summary = summarize(results, mods)
    ds_out = OUT / name
    ds_out.mkdir(parents=True, exist_ok=True)
    payload = {
        "dataset": name,
        "meta": meta,
        "modalities": mods,
        "agod_version": agod.__version__,
        "gate_semantics": {
            "means": "skip modality proj BWD when alpha_m < theta",
            "not": "skip inference / FWD",
            "analogy": "structured adapt-dropout / sparse update",
            "theta": GATE_TH,
        },
        "summary": summary,
        "trajectory": results,
    }
    jp = ds_out / "agod_imgtxt_mmd_lr.json"
    jp.write_text(json.dumps(payload, indent=2, default=float))
    dash = ds_out / "AGOD_ImgTxt_MMD_LR_Acc_Dashboard.png"
    plot_dataset(name, results, summary, mods, dash)
    return summary, mods, dash, jp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device} agod={agod.__version__}", flush=True)

    all_summary = {}
    board_rows = []
    for name in args.datasets:
        summary, mods, dash, jp = run_one(name, device)
        all_summary[name] = summary
        board_rows.append(
            f"| {name} | {','.join(mods)} | "
            f"{summary['B1']['mean_acc_lift']:+.3f} | "
            f"{summary['B5']['mean_acc_lift']:+.3f} | "
            f"{summary['B5g']['mean_acc_lift']:+.3f} | "
            f"{summary['B5r']['mean_acc_lift']:+.3f} | "
            f"{summary['B5g']['mean_flops_rel']:.3f} | "
            f"{summary['B5r']['mean_flops_rel']:.3f} | "
            f"{summary['B5g_minus_B5r_lift']:+.3f} |"
        )
        for p in [dash, jp]:
            (ART / f"{name}_{Path(p).name}").write_bytes(Path(p).read_bytes())

    board = OUT / "AGOD_ImgTxt_MultiDataset_Board.png"
    plot_board(all_summary, board)
    (ART / board.name).write_bytes(board.read_bytes())

    payload = {
        "agod_version": agod.__version__,
        "gate_semantics": "adapt-dropout (skip proj BWD); FWD inference always on",
        "datasets": all_summary,
    }
    jp_all = OUT / "agod_imgtxt_multidataset.json"
    jp_all.write_text(json.dumps(payload, indent=2, default=float))

    note = DOCS / "AGOD_imgtxt_mmd_lr_prototype.md"
    note.write_text(
        "# Image/Text AGOD portable prototype\n\n"
        "## Gate semantics\n\n"
        "- **Yes, dropout-like**: hard-gate zeros modality **proj gradients** when "
        "`α_m < θ` → structured **adapt-dropout** / sparse update.\n"
        "- **Not** \"OOD large ⇒ skip inference\": **FWD still runs**; only adapt "
        "BWD for low-α modalities is skipped. High covariate often lowers α "
        "(Acc rule: cov↑ → LR↓), so gate fires more often under strong OOD.\n\n"
        "## Smoke board\n\n"
        "| Dataset | mods | B1 | B5 | B5g | B5r | B5g flops | B5r flops | B5g−B5r |\n"
        "|---|---|---:|---:|---:|---:|---:|---:|---:|\n"
        + "\n".join(board_rows)
        + "\n\n"
        "```bash\nPYTHONPATH=. python3 scripts/run_agod_imgtxt_mmd_lr.py\n```\n"
    )
    (OUT / "README.md").write_text(note.read_text())
    (ART / note.name).write_bytes(note.read_bytes())
    print("\n=== BOARD ===", flush=True)
    print(note.read_text(), flush=True)


if __name__ == "__main__":
    main()
