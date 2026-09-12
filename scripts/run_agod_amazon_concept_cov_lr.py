#!/usr/bin/env python3
"""Amazon Reviews: concept↑ → LR↑, covariate↑ → LR↓ (accuracy-focused).

Per modality m, online window t vs reference 0:

  covariate_m = max(AUC_m - 0.5, 0) * (1 + VIMP_m)   # P(X) shift
  concept_m   = PO_m                                   # P(Y|X) shift proxy

  score_m = λ_c * z(concept) - λ_v * z(covariate)
  α_m     = Softmax(score / τ)   (EMA)
  LR_m    = lr0 * (β + (1-β) * α_m * |M|)

Policies
  B1: equal LR
  B2: old AGOD — high total shift → high LR   Softmax(AUC·VIMP+γPO)
  B3: THIS — concept↑ LR↑, covariate↑ LR↓

Metric of interest: held-out Acc lift = Acc_post - Acc_pre (and mean Acc).

  python3 scripts/run_agod_amazon_concept_cov_lr.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
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


def _z(d: dict) -> dict:
    """Shift dict values to non-neg then z-like normalize across mods (sum=1)."""
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


def decompose_shift(msg: base.MSG) -> dict:
    """Split MSG into covariate vs concept per modality."""
    cov_raw = {
        m: max(float(msg.auc[m]) - 0.5, 0.0) * (1.0 + float(msg.vimp[m]))
        for m in MODS
    }
    con_raw = {m: max(float(msg.po[m]), 0.0) for m in MODS}
    cov = _z(cov_raw)
    con = _z(con_raw)
    # score: concept boosts LR, covariate damps LR
    score = {
        m: LAM_CONCEPT * con[m] - LAM_COV * cov[m] for m in MODS
    }
    return {"cov_raw": cov_raw, "con_raw": con_raw, "cov": cov, "con": con, "score": score}


def alpha_policies(msg: base.MSG, decomp: dict) -> dict:
    # B1 uniform
    a1 = {m: 1.0 / len(MODS) for m in MODS}
    # B2 old: high total shift → high LR
    a2 = msg.alpha["B3"]  # MSG softmax from base
    # B3 new: Softmax(concept - covariate)
    a3 = _softmax(decomp["score"], TAU)
    return {"B1": a1, "B2": a2, "B3": a3}


def alpha_to_lr(alpha: dict, beta: float = BETA) -> dict:
    inv = float(len(MODS))
    out = {m: float(beta + (1.0 - beta) * alpha[m] * inv) for m in MODS}
    out["shared"] = float(np.mean([out[m] for m in MODS]))
    return out


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
        logits.append(model(imgs.to(device), txts, return_mods=False).cpu().numpy())
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
    ref = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref, device, tf)
    ema = {m: 1.0 / len(MODS) for m in MODS}
    traj = []

    for w in stream["windows"]:
        adapt, hold = split_rows(w["rows"], SEED + 13 * w["t"])
        cur_b, y_cur = base.extract_blocks(model, adapt, device, tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        decomp = decompose_shift(msg)
        alphas = alpha_policies(msg, decomp)
        raw = alphas[policy]
        for m in MODS:
            ema[m] = EMA * ema[m] + (1.0 - EMA) * raw[m]
        s = sum(ema.values())
        alpha = {m: ema[m] / s for m in MODS}
        lr_mult = alpha_to_lr(alpha)

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
            m: np.vstack([ref_b[m][-(base.N_REF // 2) :], cur_b[m][: base.N_REF // 2]])
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
            f"acc {pre['acc']:.3f}→{post['acc']:.3f} (Δ={dacc:+.3f})"
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
            "mean_concept_image": float(np.mean([r["concept"]["image"] for r in traj])),
            "mean_cov_image": float(np.mean([r["cov"]["image"] for r in traj])),
            "wins_vs_zero": int(sum(1 for x in lifts if x > 0)),
            "n_windows": len(traj),
        }
    out["B3_minus_B1_lift"] = out["B3"]["mean_acc_lift"] - out["B1"]["mean_acc_lift"]
    out["B3_minus_B2_lift"] = out["B3"]["mean_acc_lift"] - out["B2"]["mean_acc_lift"]
    out["B3_minus_B1_post"] = out["B3"]["mean_acc_post"] - out["B1"]["mean_acc_post"]
    return out


def plot_dash(results, summary, path: Path):
    fig = plt.figure(figsize=(13.2, 8.6), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.28)
    fig.suptitle(
        "Amazon: concept↑LR / covariate↓LR — accuracy lifts",
        fontsize=14,
        fontweight="bold",
    )
    b3 = results["B3"]
    ts = [r["t"] + 1 for r in b3]

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(ts, [r["concept"]["image"] for r in b3], "o-", color="#C53030", label="concept_image")
    ax.plot(ts, [r["cov"]["image"] for r in b3], "s--", color="#2B6CB0", label="covariate_image")
    ax.plot(ts, [r["concept"]["text"] for r in b3], "o-", color="#DD6B20", label="concept_text")
    ax.plot(ts, [r["cov"]["text"] for r in b3], "s--", color="#38A169", label="covariate_text")
    ax.set_title("Shift decomposition (normalized)")
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(ts, [r["lr_mult"]["image"] for r in b3], "o-", color="#C05621", label="LR× image")
    ax.plot(ts, [r["lr_mult"]["text"] for r in b3], "o-", color="#2B6CB0", label="LR× text")
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("B3 LR schedule (concept↑ / cov↓)")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 0])
    for pol, c, mk in zip(
        ["B1", "B2", "B3"], ["#718096", "#DD6B20", "#C53030"], ["o", "s", "D"]
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
    xs = np.arange(3)
    lifts = [summary[p]["mean_acc_lift"] for p in ["B1", "B2", "B3"]]
    posts = [summary[p]["mean_acc_post"] for p in ["B1", "B2", "B3"]]
    ax.bar(xs - 0.15, lifts, 0.3, color="#C53030", label="mean Acc lift")
    ax.bar(xs + 0.15, posts, 0.3, color="#4A5568", label="mean Acc post")
    ax.set_xticks(xs)
    ax.set_xticklabels(["B1 equal", "B2 shift↑LR↑", "B3 con↑/cov↓"])
    ax.axhline(0.0, color="#999", ls="--", lw=0.7)
    ax.set_title(
        f"B3−B1 lift={summary['B3_minus_B1_lift']:+.3f}  "
        f"B3−B2 lift={summary['B3_minus_B2_lift']:+.3f}"
    )
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.012,
        "Rule: score=λ_c·concept − λ_v·covariate → Softmax → LR; "
        "concept-drift adapts faster, pure covariate-shift is damped.",
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
    print(f"device={device}")
    shards = sorted(base.SHARD_DIR.glob("*.tar.gz"))
    samples = base.load_shards(shards)
    print(f"samples={len(samples)}")
    tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
        ]
    )
    stream = base.make_stream(samples)
    print(f"ref={stream['ref_cat']} wins={[w['category'] for w in stream['windows']]}")

    results = {}
    for pol in ["B1", "B2", "B3"]:
        print(f"\n===== {pol} =====")
        results[pol] = run_policy(base.AmazonAGOD, stream, device, tf, pol)

    summary = summarize(results)
    print("\n=== ACCURACY SUMMARY ===")
    print(json.dumps(summary, indent=2))

    payload = {
        "rule": "LR ∝ Softmax(λ_c·concept − λ_v·covariate)",
        "concept": "PO residual contrast (P(Y|X) shift proxy)",
        "covariate": "max(AUC-0.5,0)·(1+VIMP) (P(X) shift)",
        "summary": summary,
        "trajectory": results,
    }
    jp = OUT / "agod_amazon_concept_cov_lr.json"
    jp.write_text(json.dumps(payload, indent=2, default=float))
    dash = OUT / "AGOD_Amazon_ConceptCov_LR_Acc_Dashboard.png"
    plot_dash(results, summary, dash)

    lines = [
        r"\begin{tabular}{lcccc}",
        r"\toprule",
        r"Policy & Mean Acc lift & Mean Acc post & Mean AUC lift & Wins $>0$ \\",
        r"\midrule",
    ]
    for pol in ["B1", "B2", "B3"]:
        s = summary[pol]
        lines.append(
            f"{pol} & {s['mean_acc_lift']:+.3f} & {s['mean_acc_post']:.3f} & "
            f"{s['mean_auc_lift']:+.3f} & {s['wins_vs_zero']}/{s['n_windows']} \\\\"
        )
    lines += [
        r"\midrule",
        f"B3$-$B1 lift & {summary['B3_minus_B1_lift']:+.3f} & "
        f"{summary['B3_minus_B1_post']:+.3f} (post) & & \\\\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    tex = "\n".join(lines)
    (OUT / "AGOD_amazon_concept_cov_lr_tables_only.tex").write_text(tex)
    (DOCS / "AGOD_amazon_concept_cov_lr_tables_only.tex").write_text(tex)

    note = DOCS / "AGOD_concept_covariate_lr_amazon.md"
    note.write_text(
        "# Concept↑ LR↑ / Covariate↑ LR↓ on Amazon Reviews\n\n"
        "Per modality: `covariate = (AUC-0.5)+ · (1+VIMP)`, `concept = PO`.\n"
        "`score = λ_c·concept − λ_v·covariate` → Softmax → LR multipliers.\n\n"
        f"B3 − B1 mean Acc lift = **{summary['B3_minus_B1_lift']:+.3f}**\n"
        f"B3 − B2 mean Acc lift = **{summary['B3_minus_B2_lift']:+.3f}**\n"
        f"B3 mean Acc post = **{summary['B3']['mean_acc_post']:.3f}** "
        f"(B1={summary['B1']['mean_acc_post']:.3f})\n"
    )

    for p in [dash, jp, note]:
        (ART / Path(p).name).write_bytes(Path(p).read_bytes())
    print("dashboard", dash)
    print("B3-B1 lift", summary["B3_minus_B1_lift"])


if __name__ == "__main__":
    main()
