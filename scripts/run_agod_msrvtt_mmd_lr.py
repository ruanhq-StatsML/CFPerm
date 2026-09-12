#!/usr/bin/env python3
"""MSR-VTT portable AGOD prototype (same control plane as Amazon).

Packed feats: video/text/audio + binary labels.
  sensor  : MMD2 covariate + PO/residual concept
  state   : EMA Softmax alpha
  actuator: AdamW param-group LR_m(alpha)
  adaptor : hard-gate modality proj BWD if alpha_m < theta (B5g)

Policies: B1 equal | B2 RF con-cov | B5 MMD-cov+PO | B5g α-gate | B5r random-shutdown

  PYTHONPATH=. python3 scripts/run_agod_msrvtt_mmd_lr.py
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
from agod.lr_controller import EMARouter, alpha_to_lr, softmax_scores, z_norm
from agod.shift import decompose_hybrid, decompose_mmd, decompose_rf, residual_concept

DATA_ROOT = ROOT / "data" / "msrvtt"
OUT = ROOT / "results" / "agod_msrvtt"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_msrvtt")

MODS = ["video", "text", "audio"]
SEED = 2026
N_REF = 400
N_CUR = 350
T_WIN = 6
HOLD = 0.35
STEPS = 40
BATCH = 64
LR0 = 3e-3
EMA = 0.40
TAU = 0.30
BETA = 0.10
GATE_TH = 0.28
P_KEEP_RAND = 0.67
POLICIES = ("B1", "B2", "B5", "B5g", "B5r")
FUSE = 128
C_PF, C_PB, C_SF, C_SB = 1.0, 2.0, 0.5, 1.0


def load_pack(pack: str = "packed"):
    d = DATA_ROOT / pack
    if pack == "features_window":
        feats = {
            "video": np.load(d / "video_feat.npy").astype(np.float32),
            "text": np.load(d / "text_feat.npy").astype(np.float32),
            "audio": np.load(d / "audio_feat.npy").astype(np.float32),
        }
        y = np.load(d / "video_labels.npy").astype(np.int64)
    else:
        feats = {
            "video": np.load(d / "video_feat.npy").astype(np.float32),
            "text": np.load(d / "text_feat.npy").astype(np.float32),
            "audio": np.load(d / "audio_feat.npy").astype(np.float32),
        }
        y = np.load(d / "labelsmsr.npy").astype(np.int64)
    # binary if needed
    if y.ndim > 1:
        y = y.reshape(len(y), -1)[:, 0]
    if len(np.unique(y)) > 2:
        mode = int(np.bincount(y.astype(int)).argmax())
        y = (y == mode).astype(np.int64)
    else:
        y = y.astype(np.int64)
    return feats, y, d


def make_stream(y):
    rng = np.random.default_rng(SEED)
    idx0 = np.where(y == 0)[0].copy()
    idx1 = np.where(y == 1)[0].copy()
    rng.shuffle(idx0)
    rng.shuffle(idx1)
    n0, n1 = N_REF // 2, N_REF - N_REF // 2
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


def blocks_from_idx(feats, y, idx):
    return {m: feats[m][idx] for m in MODS}, y[idx].astype(float)


def domain_msg(b0, b1, y0, y1, *, seed: int):
    auc, vimp, po = {}, {}, {}
    for i, m in enumerate(MODS):
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
                n_estimators=60,
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


def flops_rel(active: dict) -> float:
    fwd = len(MODS) * C_PF + C_SF
    bwd = sum(C_PB for m in MODS if active[m]) + C_SB
    full = len(MODS) * C_PF + C_SF + len(MODS) * C_PB + C_SB
    return float((fwd + bwd) / full)


class MSRFusion(nn.Module):
    def __init__(self):
        super().__init__()
        self.video_proj = nn.Sequential(nn.Linear(768, FUSE), nn.ReLU(), nn.Dropout(0.1))
        self.text_proj = nn.Sequential(nn.Linear(768, FUSE), nn.ReLU(), nn.Dropout(0.1))
        self.audio_proj = nn.Sequential(nn.Linear(512, FUSE), nn.ReLU(), nn.Dropout(0.1))
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


def build_optim(model: MSRFusion):
    g = model.param_groups()
    return torch.optim.AdamW(
        [
            {"params": g["video"], "lr": LR0, "name": "video"},
            {"params": g["text"], "lr": LR0, "name": "text"},
            {"params": g["audio"], "lr": LR0, "name": "audio"},
            {"params": g["shared"], "lr": LR0, "name": "shared"},
        ]
    )


def set_lrs(opt, mult: dict):
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
        return {"acc": float("nan"), "auc": float("nan"), "n": len(idx)}
    logits, ys = [], []
    for s in range(0, len(idx), BATCH):
        sl = idx[s : s + BATCH]
        batch, yy = to_batch(feats, y, sl, device)
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


def train_window(model, opt, feats, y, idx, device, *, lr_mult, active):
    model.train()
    set_lrs(opt, lr_mult)
    crit = nn.CrossEntropyLoss()
    losses = []
    rng = np.random.default_rng(SEED + int(np.asarray(idx).sum()) % 100000)
    t0 = time.perf_counter()
    for _ in range(STEPS):
        sel = rng.choice(idx, size=min(BATCH, len(idx)), replace=False)
        batch, yy = to_batch(feats, y, sel, device)
        loss = crit(model(batch), yy)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        groups = model.param_groups()
        for m in MODS:
            if not active[m]:
                for p in groups[m]:
                    p.grad = None
        opt.step()
        losses.append(float(loss.item()))
    return float(np.mean(losses)), (time.perf_counter() - t0) * 1000.0


def select_raw(policy, msg, b0, b1, y0, y1, *, seed):
    unit = {m: 1.0 for m in MODS}
    if policy == "B1":
        decomp = {"cov": unit, "con": unit, "score": {m: 0.0 for m in MODS}}
        return {m: 1.0 / len(MODS) for m in MODS}, decomp
    if policy == "B2":
        decomp = decompose_rf(msg, MODS)
        return softmax_scores(decomp["score"], MODS, TAU), decomp
    mmd_d = decompose_mmd(b0, b1, y0, y1, MODS, seed=seed)
    hy = decompose_hybrid(msg, mmd_d, MODS)
    return softmax_scores(hy["score"], MODS, TAU), hy


def run_policy(feats, y, stream, device, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = MSRFusion().to(device)
    opt = build_optim(model)
    router = EMARouter(MODS, ema=EMA)
    ref_idx = stream["ref_idx"].copy()
    ref_b, y_ref = blocks_from_idx(feats, y, ref_idx)

    warm, _ = split_hold(ref_idx, SEED)
    train_window(
        model,
        opt,
        feats,
        y,
        warm,
        device,
        lr_mult={**{m: 1.0 for m in MODS}, "shared": 1.0},
        active={m: True for m in MODS},
    )

    traj = []
    for w in stream["windows"]:
        adapt, hold = split_hold(w["idx"], SEED + 13 * w["t"])
        cur_b, y_cur = blocks_from_idx(feats, y, adapt)
        msg = domain_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        raw, decomp = select_raw(
            policy, msg, ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"]
        )
        alpha = router.update(raw)
        lr_mult = alpha_to_lr(alpha, MODS, beta=BETA)

        if policy == "B5g":
            active = {m: alpha[m] >= GATE_TH for m in MODS}
            if not any(active.values()):
                mstar = max(MODS, key=lambda m: alpha[m])
                active = {m: m == mstar for m in MODS}
        elif policy == "B5r":
            rng_g = np.random.default_rng(SEED + 101 * int(w["t"]) + 17)
            active = {m: bool(rng_g.random() < P_KEEP_RAND) for m in MODS}
            if not any(active.values()):
                active[MODS[int(rng_g.integers(0, len(MODS)))]] = True
        else:
            active = {m: True for m in MODS}

        pre = eval_acc(model, feats, y, hold, device)
        loss, wall = train_window(
            model, opt, feats, y, adapt, device, lr_mult=lr_mult, active=active
        )
        post = eval_acc(model, feats, y, hold, device)
        dacc = post["acc"] - pre["acc"]
        dauc = (
            post["auc"] - pre["auc"]
            if np.isfinite(pre["auc"]) and np.isfinite(post["auc"])
            else float("nan")
        )
        fr = flops_rel(active)
        util = dacc / fr if fr > 0 else float("nan")

        keep = N_REF // 2
        ref_idx = np.concatenate([ref_idx[-keep:], adapt[: min(keep, len(adapt))]])
        ref_b, y_ref = blocks_from_idx(feats, y, ref_idx)

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
            "auc_lift": dauc,
            "cost_utility": util,
            "train_loss": loss,
        }
        traj.append(row)
        print(
            f"[{policy}] t={w['t']} p1={w['p1']:.2f} "
            f"a={ {m: round(alpha[m], 2) for m in MODS} } "
            f"LR*={ {m: round(lr_mult[m], 2) for m in MODS} } "
            f"sel={ {m: int(active[m]) for m in MODS} } "
            f"FLOPs={fr:.2f} | acc {pre['acc']:.3f}->{post['acc']:.3f} (d={dacc:+.3f})",
            flush=True,
        )
    return traj


def summarize(results):
    out = {}
    for pol, traj in results.items():
        lifts = [r["acc_lift"] for r in traj]
        out[pol] = {
            "mean_acc_lift": float(np.nanmean(lifts)),
            "mean_acc_post": float(np.nanmean([r["acc_post"] for r in traj])),
            "mean_flops_rel": float(np.mean([r["flops_rel"] for r in traj])),
            "mean_wall_ms": float(np.mean([r["wall_ms"] for r in traj])),
            "mean_cost_utility": float(np.nanmean([r["cost_utility"] for r in traj])),
            "mean_lr_video": float(np.mean([r["lr_mult"]["video"] for r in traj])),
            "mean_lr_text": float(np.mean([r["lr_mult"]["text"] for r in traj])),
            "mean_lr_audio": float(np.mean([r["lr_mult"]["audio"] for r in traj])),
            "frac_video_sel": float(np.mean([r["selected"]["video"] for r in traj])),
            "frac_text_sel": float(np.mean([r["selected"]["text"] for r in traj])),
            "frac_audio_sel": float(np.mean([r["selected"]["audio"] for r in traj])),
            "wins_vs_zero": int(sum(1 for x in lifts if x > 0)),
            "n_windows": len(traj),
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


def plot_dash(results, summary, path: Path):
    pols = list(POLICIES)
    fig = plt.figure(figsize=(14.0, 9.0), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.28)
    fig.suptitle(
        "MSR-VTT AGOD prototype: MMD cov + PO concept (+ adaptor gate)",
        fontsize=14,
        fontweight="bold",
    )
    b5 = results["B5"]
    ts = [r["t"] + 1 for r in b5]

    ax = fig.add_subplot(gs[0, 0])
    for m, c in zip(MODS, ["#C53030", "#2B6CB0", "#38A169"]):
        ax.plot(ts, [r["concept"][m] for r in b5], "o-", color=c, label=f"con_{m}")
        ax.plot(ts, [r["cov"][m] for r in b5], "s--", color=c, alpha=0.45, label=f"cov_{m}")
    ax.set_title("B5 sensors (PO concept / MMD cov)")
    ax.legend(frameon=False, fontsize=6, ncol=3)

    ax = fig.add_subplot(gs[0, 1])
    for m, c in zip(MODS, ["#C53030", "#2B6CB0", "#38A169"]):
        ax.plot(ts, [r["lr_mult"][m] for r in b5], "o-", color=c, label=f"LR* {m}")
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("B5 actuator (param-group LR)")
    ax.legend(frameon=False, fontsize=7)

    ax = fig.add_subplot(gs[1, 0])
    for pol, c, mk in zip(pols, ["#718096", "#DD6B20", "#805AD5", "#C53030"], ["o", "s", "D", "P"]):
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
        f"B5-B1 lift={summary['B5_minus_B1_lift']:+.3f}  "
        f"B5g flops={summary['B5g']['mean_flops_rel']:.2f}"
    )
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.012,
        f"portable MSR-VTT use-case | agod {agod.__version__} | "
        "telemetry: cov/con/alpha/LR/sel/flops/Acc",
        ha="center",
        fontsize=8.5,
        color="#333",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", default="packed", choices=["packed", "features_pack", "features_window"])
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device} agod={agod.__version__} pack={args.pack}", flush=True)
    feats, y, data_dir = load_pack(args.pack)
    print(
        f"n={len(y)} dims=video{feats['video'].shape[1]}/"
        f"text{feats['text'].shape[1]}/audio{feats['audio'].shape[1]}",
        flush=True,
    )
    global N_REF, N_CUR
    if len(y) < N_REF + N_CUR:
        N_REF = min(N_REF, max(80, len(y)//5))
        N_CUR = min(N_CUR, max(60, len(y)//6))
        print(f"small pack: N_REF={N_REF} N_CUR={N_CUR}", flush=True)

    results = {}
    for pol in POLICIES:
        print(f"\n===== {pol} =====", flush=True)
        stream = make_stream(y)
        results[pol] = run_policy(feats, y, stream, device, pol)

    summary = summarize(results)
    print("\n=== MSR-VTT AGOD SUMMARY ===", flush=True)
    print(json.dumps(summary, indent=2), flush=True)

    payload = {
        "dataset": f"MSR-VTT {args.pack} multimodal features",
        "agod_version": agod.__version__,
        "control_plane": {
            "sensor": "MMD2 cov + PO/residual concept",
            "state": "EMA Softmax alpha",
            "actuator": "AdamW param-group LR_m(alpha)",
            "adaptor": "hard-gate proj BWD if alpha_m < theta (B5g)",
        },
        "summary": summary,
        "trajectory": results,
    }
    jp = OUT / f"agod_msrvtt_{args.pack}_mmd_lr.json"
    jp.write_text(json.dumps(payload, indent=2, default=float))
    dash = OUT / f"AGOD_MSRVTT_{args.pack}_MMD_LR_Acc_Dashboard.png"
    plot_dash(results, summary, dash)

    note = DOCS / "AGOD_msrvtt_mmd_lr_prototype.md"
    rows = "\n".join(
        f"| {p} | {summary[p]['mean_acc_lift']:+.3f} | "
        f"{summary[p]['mean_flops_rel']:.3f} | "
        f"{summary[p]['mean_cost_utility']:+.3f} |"
        for p in POLICIES
    )
    note.write_text(
        "# MSR-VTT AGOD portable prototype\n\n"
        "Same control plane as Amazon (`agod/`): sensor -> EMA state -> "
        "param-group LR actuator -> optional hard-gate adaptor.\n\n"
        "| Policy | Acc lift | flops_rel | cost_utility |\n"
        "|---|---:|---:|---:|\n"
        f"{rows}\n\n"
        f"- B5-B1 Acc lift = **{summary['B5_minus_B1_lift']:+.3f}**\n"
        f"- B5g flops_rel = **{summary['B5g']['mean_flops_rel']:.3f}** "
        f"(vs B1 1.000)\n\n"
        "```bash\nPYTHONPATH=. python3 scripts/run_agod_msrvtt_mmd_lr.py\n```\n"
    )
    (OUT / "README.md").write_text(note.read_text())
    for p in [dash, jp, note]:
        (ART / Path(p).name).write_bytes(Path(p).read_bytes())
    print("dashboard", dash, flush=True)
    print(
        "B5-B1",
        summary["B5_minus_B1_lift"],
        "B5g-flops",
        summary["B5g"]["mean_flops_rel"],
        flush=True,
    )


if __name__ == "__main__":
    main()
