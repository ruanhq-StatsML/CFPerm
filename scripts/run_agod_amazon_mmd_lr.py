#!/usr/bin/env python3
"""Amazon Reviews: MMD concept↑LR / covariate↓LR (accuracy-focused).

MMD logic (FSDS spirit)
-----------------------
  covariate_m = MMD²(X_m^ref, X_m^cur)                 # P(X) shift in RKHS
  joint_m     = MMD²([X_m, Y]^ref, [X_m, Y]^cur)
  concept_m   = max(joint_m - covariate_m, 0)          # label-coupled excess
                + residual gap of f_ref: Y ~ X_m        # P(Y|X) change

  score_m = λ_c · z(concept) − λ_v · z(covariate)
  α       = Softmax(score / τ)   (EMA)
  LR_m    = lr0 · (β + (1−β)·α_m·|M|) · gain_m
  gain_m  = 1 + κ · relu(concept_z − covariate_z)      # B4 only

Policies
  B1: equal LR
  B2: RF AUC/PO concept−cov (previous Amazon rule)
  B3: MMD concept−cov → Softmax → LR
  B4: B3 + concept-intensity LR gain

  python3 scripts/run_agod_amazon_mmd_lr.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import accuracy_score, roc_auc_score
from torch.utils.data import DataLoader
from torchvision import transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_agod_amazon_modality_lr as base  # noqa: E402

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
LAM_CONCEPT = 1.0
LAM_COV = 1.0
KAPPA = 1.25
MMD_MAX_N = 96
RESID_TREES = 40


def _z(d: dict) -> dict:
    v = np.array([max(float(d[m]), 0.0) for m in MODS], float)
    if v.sum() <= 1e-12:
        return {m: 1.0 / len(MODS) for m in MODS}
    v = v / v.sum()
    return {m: float(v[i]) for i, m in enumerate(MODS)}


def _softmax(d: dict, tau: float) -> dict:
    z = np.array([d[m] for m in MODS], float) / max(tau, 1e-6)
    z = z - z.max()
    e = np.exp(z)
    e = e / e.sum()
    return {m: float(e[i]) for i, m in enumerate(MODS)}


def _whiten_pair(X0: np.ndarray, X1: np.ndarray):
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    pooled = np.concatenate([X0, X1], 0)
    mu = pooled.mean(0, keepdims=True)
    sd = pooled.std(0, keepdims=True) + 1e-6
    return (X0 - mu) / sd, (X1 - mu) / sd


def rbf_mmd2(X0, X1, *, max_n: int = MMD_MAX_N, seed: int = 0) -> float:
    """Unbiased RBF MMD² (median bandwidth) — same estimator as FSDS."""
    rng = np.random.default_rng(seed)
    X0 = np.asarray(X0, float)
    X1 = np.asarray(X1, float)
    if len(X0) > max_n:
        X0 = X0[rng.choice(len(X0), max_n, replace=False)]
    if len(X1) > max_n:
        X1 = X1[rng.choice(len(X1), max_n, replace=False)]
    if len(X0) < 4 or len(X1) < 4:
        return 0.0
    Z = np.vstack([X0, X1])
    idx = rng.choice(len(Z), size=min(256, len(Z)), replace=False)
    S = Z[idx]
    d2 = np.sum((S[:, None, :] - S[None, :, :]) ** 2, axis=-1)
    med = float(np.median(d2[d2 > 0])) if np.any(d2 > 0) else 1.0
    gamma = 1.0 / max(med, 1e-8)

    def k(A, B):
        dd = np.sum((A[:, None, :] - B[None, :, :]) ** 2, axis=-1)
        return np.exp(-gamma * dd)

    K00, K11, K01 = k(X0, X0), k(X1, X1), k(X0, X1)
    n0, n1 = len(X0), len(X1)
    mmd2 = (
        (K00.sum() - np.trace(K00)) / max(n0 * (n0 - 1), 1)
        + (K11.sum() - np.trace(K11)) / max(n1 * (n1 - 1), 1)
        - 2.0 * K01.mean()
    )
    return float(max(mmd2, 0.0))


def residual_concept(X0, Y0, X1, Y1, *, seed: int) -> float:
    """P(Y|X) proxy: f_ref fit on ref; error gap + residual MMD on cur."""
    if len(X0) < 12 or len(X1) < 12:
        return 0.0
    rf = RandomForestRegressor(
        n_estimators=RESID_TREES,
        max_depth=5,
        min_samples_leaf=3,
        random_state=seed,
        n_jobs=1,
    )
    rf.fit(X0, Y0)
    e0 = float(np.mean(np.abs(Y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(Y1 - rf.predict(X1))))
    r0 = (Y0 - rf.predict(X0)).reshape(-1, 1)
    r1 = (Y1 - rf.predict(X1)).reshape(-1, 1)
    mmd_r = rbf_mmd2(r0, r1, max_n=min(MMD_MAX_N, 64), seed=seed + 3)
    return max(e1 - e0, 0.0) + 0.5 * mmd_r


def decompose_mmd(b0, b1, y0, y1, *, seed: int) -> dict:
    y0 = np.asarray(y0, float)
    y1 = np.asarray(y1, float)
    cov_raw, con_raw, joint_raw, resid_raw = {}, {}, {}, {}
    y_all = np.concatenate([y0, y1])
    ys = y_all.std() + 1e-6
    y0n = ((y0 - y_all.mean()) / ys).reshape(-1, 1)
    y1n = ((y1 - y_all.mean()) / ys).reshape(-1, 1)
    for i, m in enumerate(MODS):
        X0, X1 = _whiten_pair(b0[m], b1[m])
        cov = rbf_mmd2(X0, X1, seed=seed + i)
        joint = rbf_mmd2(
            np.hstack([X0, y0n]), np.hstack([X1, y1n]), seed=seed + 17 + i
        )
        excess = max(joint - cov, 0.0)
        resid = residual_concept(X0, y0, X1, y1, seed=seed + 31 + i)
        cov_raw[m] = cov
        joint_raw[m] = joint
        resid_raw[m] = resid
        con_raw[m] = excess + resid
    cov, con = _z(cov_raw), _z(con_raw)
    score = {m: LAM_CONCEPT * con[m] - LAM_COV * cov[m] for m in MODS}
    return {
        "cov_raw": cov_raw,
        "con_raw": con_raw,
        "joint_raw": joint_raw,
        "resid_raw": resid_raw,
        "cov": cov,
        "con": con,
        "score": score,
    }


def decompose_rf(msg: base.MSG) -> dict:
    cov_raw = {
        m: max(float(msg.auc[m]) - 0.5, 0.0) * (1.0 + float(msg.vimp[m]))
        for m in MODS
    }
    con_raw = {m: max(float(msg.po[m]), 0.0) for m in MODS}
    cov, con = _z(cov_raw), _z(con_raw)
    score = {m: LAM_CONCEPT * con[m] - LAM_COV * cov[m] for m in MODS}
    return {
        "cov_raw": cov_raw,
        "con_raw": con_raw,
        "cov": cov,
        "con": con,
        "score": score,
    }




def decompose_hybrid(msg: base.MSG, mmd_d: dict) -> dict:
    """FSDS-style: MMD² for covariate (P(X)), PO for concept (P(Y|X))."""
    cov_raw = dict(mmd_d["cov_raw"])
    con_raw = {m: max(float(msg.po[m]), 0.0) for m in MODS}
    cov, con = _z(cov_raw), _z(con_raw)
    score = {m: LAM_CONCEPT * con[m] - LAM_COV * cov[m] for m in MODS}
    return {
        "cov_raw": cov_raw,
        "con_raw": con_raw,
        "cov": cov,
        "con": con,
        "score": score,
    }

def alpha_to_lr(alpha: dict, gain: dict | None = None) -> dict:
    inv = float(len(MODS))
    gain = gain or {m: 1.0 for m in MODS}
    out = {
        m: float((BETA + (1.0 - BETA) * alpha[m] * inv) * gain[m]) for m in MODS
    }
    out["shared"] = float(np.mean([out[m] for m in MODS]))
    return out


def intensity_gain(decomp: dict) -> dict:
    return {
        m: float(1.0 + KAPPA * max(decomp["con"][m] - decomp["cov"][m], 0.0))
        for m in MODS
    }


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
    """B1 equal | B2 RF | B3 MMD | B4 MMD+gain | B5 MMD-cov+PO-concept (FSDS)."""
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = ctor().to(device)
    opt = base.build_optim(model)
    ref = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref, device, tf)
    ema = {m: 1.0 / len(MODS) for m in MODS}
    traj = []

    for w in stream["windows"]:
        adapt, hold = split_rows(w["rows"], SEED + 13 * w["t"])
        cur_b, y_cur = base.extract_blocks(model, adapt, device, tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        rf_d = decompose_rf(msg)
        mmd_d = decompose_mmd(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])

        if policy == "B1":
            raw = {m: 1.0 / len(MODS) for m in MODS}
            decomp, gain = mmd_d, {m: 1.0 for m in MODS}
        elif policy == "B2":
            raw = _softmax(rf_d["score"], TAU)
            decomp, gain = rf_d, {m: 1.0 for m in MODS}
        elif policy == "B3":
            raw = _softmax(mmd_d["score"], TAU)
            decomp, gain = mmd_d, {m: 1.0 for m in MODS}
        elif policy == "B4":
            raw = _softmax(mmd_d["score"], TAU)
            decomp, gain = mmd_d, intensity_gain(mmd_d)
        else:  # B5 FSDS hybrid: MMD covariate + PO concept
            hy_d = decompose_hybrid(msg, mmd_d)
            raw = _softmax(hy_d["score"], TAU)
            decomp, gain = hy_d, {m: 1.0 for m in MODS}

        for m in MODS:
            ema[m] = EMA * ema[m] + (1.0 - EMA) * raw[m]
        s = sum(ema.values())
        alpha = {m: ema[m] / s for m in MODS}
        lr_mult = alpha_to_lr(alpha, gain=gain)

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
                "mmd_cov": mmd_d["cov"],
                "mmd_concept": mmd_d["con"],
                "mmd_cov_raw": mmd_d["cov_raw"],
                "mmd_con_raw": mmd_d["con_raw"],
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
            "sum_acc_lift": float(np.nansum(lifts)),
            "mean_acc_post": float(np.nanmean(posts)),
            "mean_auc_lift": float(np.nanmean([r["auc_lift"] for r in traj])),
            "mean_lr_image": float(np.mean([r["lr_mult"]["image"] for r in traj])),
            "mean_lr_text": float(np.mean([r["lr_mult"]["text"] for r in traj])),
            "wins_vs_zero": int(sum(1 for x in lifts if x > 0)),
            "n_windows": len(traj),
        }
    for a, b in [
        ("B3", "B1"), ("B3", "B2"),
        ("B4", "B1"), ("B4", "B2"), ("B4", "B3"),
        ("B5", "B1"), ("B5", "B2"), ("B5", "B3"), ("B5", "B4"),
    ]:
        out[f"{a}_minus_{b}_lift"] = out[a]["mean_acc_lift"] - out[b]["mean_acc_lift"]
        out[f"{a}_minus_{b}_post"] = out[a]["mean_acc_post"] - out[b]["mean_acc_post"]
    return out


def plot_dash(results, summary, path: Path):
    fig = plt.figure(figsize=(14.0, 9.0), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.28)
    fig.suptitle(
        "Amazon MMD: concept↑LR / covariate↓LR — Acc lifts",
        fontsize=14,
        fontweight="bold",
    )
    b4 = results["B4"]
    ts = [r["t"] + 1 for r in b4]

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(ts, [r["mmd_concept"]["image"] for r in b4], "o-", color="#C53030", label="MMD concept_image")
    ax.plot(ts, [r["mmd_cov"]["image"] for r in b4], "s--", color="#2B6CB0", label="MMD cov_image")
    ax.plot(ts, [r["mmd_concept"]["text"] for r in b4], "o-", color="#DD6B20", label="MMD concept_text")
    ax.plot(ts, [r["mmd_cov"]["text"] for r in b4], "s--", color="#38A169", label="MMD cov_text")
    ax.set_title("MMD shift decomposition (normalized)")
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(ts, [r["lr_mult"]["image"] for r in b4], "o-", color="#C05621", label="B4 LR× image")
    ax.plot(ts, [r["lr_mult"]["text"] for r in b4], "o-", color="#2B6CB0", label="B4 LR× text")
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("B4 LR (MMD score × intensity gain)")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 0])
    for pol, c, mk in zip(
        ["B1", "B2", "B3", "B4", "B5"],
        ["#718096", "#DD6B20", "#2B6CB0", "#C53030", "#805AD5"],
        ["o", "s", "^", "D", "P"],
    ):
        ax.plot(
            [r["t"] + 1 for r in results[pol]],
            [r["acc_lift"] for r in results[pol]],
            marker=mk,
            color=c,
            label=pol,
        )
    ax.axhline(0.0, color="#999", ls="--", lw=0.8)
    ax.set_title("Held-out Acc lift per window")
    ax.set_ylabel("Acc_post − Acc_pre")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    xs = np.arange(5)
    lifts = [summary[p]["mean_acc_lift"] for p in ["B1", "B2", "B3", "B4", "B5"]]
    posts = [summary[p]["mean_acc_post"] for p in ["B1", "B2", "B3", "B4", "B5"]]
    ax.bar(xs - 0.15, lifts, 0.3, color="#C53030", label="mean Acc lift")
    ax.bar(xs + 0.15, posts, 0.3, color="#4A5568", label="mean Acc post")
    ax.set_xticks(xs)
    ax.set_xticklabels(["B1", "B2 RF", "B3 MMD", "B4+gain", "B5 hyb"])
    ax.axhline(0.0, color="#999", ls="--", lw=0.7)
    ax.set_title(
        f"B5−B1={summary['B5_minus_B1_lift']:+.3f}  "
        f"B5−B2={summary['B5_minus_B2_lift']:+.3f}  "
        f"B5−B3={summary['B5_minus_B3_lift']:+.3f}"
    )
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.012,
        "MMD cov=MMD²(X); concept=max(MMD²([X,Y])−MMD²(X),0)+residual gap; "
        "LR∝Softmax(λc·con−λv·cov)·(1+κ·relu(con−cov))",
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
    print(f"device={device}", flush=True)
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
    for pol in ["B1", "B2", "B3", "B4", "B5"]:
        print(f"\n===== {pol} =====", flush=True)
        results[pol] = run_policy(base.AmazonAGOD, stream, device, tf, pol)

    summary = summarize(results)
    print("\n=== ACCURACY SUMMARY ===", flush=True)
    print(json.dumps(summary, indent=2), flush=True)

    payload = {
        "rule": "MMD cov + joint-excess/residual concept → Softmax(con−cov) → LR",
        "mmd": {
            "covariate": "unbiased RBF MMD²(X_m^ref, X_m^cur)",
            "concept": "max(MMD²([X,Y])−MMD²(X),0) + residual(f_ref) gap/MMD",
            "gain_B4": "1 + κ·relu(concept_z − covariate_z)",
        },
        "summary": summary,
        "trajectory": results,
    }
    jp = OUT / "agod_amazon_mmd_lr.json"
    jp.write_text(json.dumps(payload, indent=2, default=float))
    dash = OUT / "AGOD_Amazon_MMD_LR_Acc_Dashboard.png"
    plot_dash(results, summary, dash)

    lines = [
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Policy & Mean Acc lift & Mean Acc post & Mean AUC lift & Wins $>0$ \\",
        r"\midrule",
    ]
    labels = {
        "B1": "B1 equal",
        "B2": "B2 RF con$-$cov",
        "B3": "B3 MMD con$-$cov",
        "B4": "B4 MMD+gain",
        "B5": "B5 MMD-cov+PO",
    }
    for pol in ["B1", "B2", "B3", "B4", "B5"]:
        s = summary[pol]
        lines.append(
            f"{labels[pol]} & {s['mean_acc_lift']:+.3f} & {s['mean_acc_post']:.3f} & "
            f"{s['mean_auc_lift']:+.3f} & {s['wins_vs_zero']}/{s['n_windows']} \\\\"
        )
    lines += [
        r"\midrule",
        f"B3$-$B2 lift & {summary['B3_minus_B2_lift']:+.3f} & "
        f"{summary['B3_minus_B2_post']:+.3f} (post) & & \\\\",
        f"B4$-$B2 lift & {summary['B4_minus_B2_lift']:+.3f} & "
        f"{summary['B4_minus_B2_post']:+.3f} (post) & & \\\\",
        f"B4$-$B1 lift & {summary['B4_minus_B1_lift']:+.3f} & "
        f"{summary['B4_minus_B1_post']:+.3f} (post) & & \\\\",
        f"B5$-$B2 lift & {summary['B5_minus_B2_lift']:+.3f} & "
        f"{summary['B5_minus_B2_post']:+.3f} (post) & & \\\\",
        f"B5$-$B1 lift & {summary['B5_minus_B1_lift']:+.3f} & "
        f"{summary['B5_minus_B1_post']:+.3f} (post) & & \\\\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    tex = "\n".join(lines)
    (OUT / "AGOD_amazon_mmd_lr_tables_only.tex").write_text(tex)
    (DOCS / "AGOD_amazon_mmd_lr_tables_only.tex").write_text(tex)

    note = DOCS / "AGOD_mmd_concept_covariate_lr_amazon.md"
    note.write_text(
        "# MMD-informed concept↑ / covariate↓ LR on Amazon\n\n"
        "## Why MMD can improve the LR signal\n\n"
        "RF Domain AUC/VIMP asks *can a classifier separate domains?* — a useful but "
        "indirect `P(X)` proxy that saturates near 0.5/0.5 when both modalities drift.\n\n"
        "Unbiased RBF **MMD²** measures RKHS mean embedding distance directly "
        "(same estimator as FSDS MMD-LOCO). Per modality:\n\n"
        "- **covariate** = `MMD²(X_m^ref, X_m^cur)` → damp LR\n"
        "- **concept** = `max(MMD²([X,Y]) − MMD²(X), 0) + residual(f_ref)` → raise LR\n\n"
        "Joint-excess isolates shift that only appears when labels couple to features; "
        "residual gap of `f_ref: Y~X` catches `P(Y|X)` movement even when `P(X)` is stable.\n\n"
        "B4 multiplies LR by `1+κ·relu(concept−cov)` so concept-dominated modalities "
        "take **larger** steps (accuracy-seeking).\n\n"
        f"## Smoke Acc\n\n"
        f"| Policy | Mean Acc lift | Mean Acc post |\n"
        f"|---|---:|---:|\n"
        f"| B1 equal | {summary['B1']['mean_acc_lift']:+.3f} | {summary['B1']['mean_acc_post']:.3f} |\n"
        f"| B2 RF con−cov | {summary['B2']['mean_acc_lift']:+.3f} | {summary['B2']['mean_acc_post']:.3f} |\n"
        f"| B3 MMD con−cov | {summary['B3']['mean_acc_lift']:+.3f} | {summary['B3']['mean_acc_post']:.3f} |\n"
        f"| **B4 MMD+gain** | **{summary['B4']['mean_acc_lift']:+.3f}** | **{summary['B4']['mean_acc_post']:.3f}** |\n\n"
        f"- B3−B2 lift = **{summary['B3_minus_B2_lift']:+.3f}**\n"
        f"- B4−B2 lift = **{summary['B4_minus_B2_lift']:+.3f}**\n"
        f"- B4−B1 lift = **{summary['B4_minus_B1_lift']:+.3f}**\n"
    )

    for p in [dash, jp, note]:
        (ART / Path(p).name).write_bytes(Path(p).read_bytes())
    print("dashboard", dash, flush=True)
    print(
        "B5-B1", summary["B5_minus_B1_lift"],
        "B5-B2", summary["B5_minus_B2_lift"],
        "B3-B2", summary["B3_minus_B2_lift"],
        "B5-B3", summary["B5_minus_B3_lift"],
        flush=True,
    )


if __name__ == "__main__":
    main()
