#!/usr/bin/env python3
"""Amazon AGOD: modality LR action + upstream vs downstream metrics.

Action (only): MSG → α_m → LR_m / grad-gate.
No sample / patch-token actions.

Upstream: AUC_m, VIMP_m, PO_m, g_m, α_m, LR×_m, adapt FLOPs
Downstream: holdout Acc/AUC pre→post; reference Acc delta (forgetting)

Also writes a logic note on sample / patch-token attribution grains
(measurement only — not wired into the optimizer).

  python3 scripts/run_agod_amazon_updownstream.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
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
EMA = base.EMA
BATCH = base.BATCH
HOLD_FRAC = 0.35

LOGIC_MD = r"""# Attribution grains: modality / sample / patch-token

This note separates **what each grain measures** from **what the Amazon
prototype actually controls**. Only the modality grain is closed-loop.

## 1. Three grains = three questions

| Grain | Question | Statistic | Matching action |
|---|---|---|---|
| **Modality** | Which channel drifts? | MSG \(g_m\) on block \(X_m\) | **LR\(_m\) / gate \(\\partial L/\\partial\\theta_m\)** ← implemented |
| **Sample** | Which examples carry the shift? | domain score / residual on row \(i\) | sample weight / replay — *not* LR\(_m\) |
| **Patch / token** | Where inside a modality? | LOO / occlusion \(\\Delta_k\) | patch mask / token distill weight — *not* LR\(_m\) |

## 2. Upstream vs downstream (different jobs)

- **Upstream** (shift diagnostics → control): AUC\(_m\), VIMP, PO, \(g_m\), \(\\alpha_m\), LR×, adapt FLOPs.
- **Downstream** (task utility): holdout Acc/AUC Δ after the update; reference Acc Δ (forgetting).

They need not move together. High \(g_{\\text{image}}\) does not guarantee larger holdout ΔAcc.

## 3. Modality layer (implemented control law)

```
g_m = Normalize(AUC_m · VIMP_m + γ · PO_m)
α   = Softmax(g / τ)          # EMA
LR_m = lr0 · (β + (1-β) · α_m · |M|)
if α_m < θ: zero grads of modality-m projection
```

## 4. Sample layer (logic only)

On a fixed modality representation (or concat):
- domain RF: reference vs current;
- score \(s_i = P(W{=}1\\mid x_i)-1/2\) (or residual risk).
High \(s_i\) = row looks most out-of-reference.

**Do not** fold \(s_i\) into Softmax→LR\(_m\). That mixes grains.
If ever wired, action = \(w_i=w(s_i)\) inside \(\\sum_i w_i \\ell_i\) or replay buffer.

## 5. Patch / token layer (logic only)

Inside the already-chosen modality:
- leave-one-patch-out / token occlusion → \(\\Delta_k\) of domain score or PO;
- rank / visualize top patches or tokens.

**Do not** re-decide LR\(_m\). If ever wired, action = mask / reweight that
modality's local distill terms.

## 6. Why this prototype stops at modality

The actuator is already \(\\alpha_m\\to\\mathrm{LR}_m\).
Sample/patch scores answer finer localization questions; without a matching
actuator they are diagnostics, not “secondary AGOD”.
"""


def split_rows(rows, frac=HOLD_FRAC, seed=SEED):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    n_hold = max(20, int(len(rows) * frac))
    n_hold = min(n_hold, len(rows) // 2)
    hold = [rows[i] for i in idx[:n_hold]]
    adapt = [rows[i] for i in idx[n_hold:]]
    return adapt, hold


@torch.no_grad()
def eval_downstream(model, rows, device, image_tf):
    model.eval()
    if len(rows) < 8:
        return {"acc": float("nan"), "auc": float("nan"), "n": len(rows)}
    loader = DataLoader(base.AmazonDS(rows, image_tf), batch_size=BATCH, shuffle=False)
    logits_all, y_all = [], []
    for imgs, txts, y in loader:
        imgs = imgs.to(device)
        logits = model(imgs, txts, return_mods=False)
        logits_all.append(logits.cpu().numpy())
        y_all.append(y.numpy())
    logits = np.concatenate(logits_all)
    y = np.concatenate(y_all)
    acc = float(accuracy_score(y, logits.argmax(1)))
    auc = float("nan")
    if len(np.unique(y)) > 1:
        try:
            prob = torch.softmax(torch.from_numpy(logits), 1).numpy()[:, 1]
            auc = float(roc_auc_score(y, prob))
        except Exception:
            pass
    return {"acc": acc, "auc": auc, "n": int(len(y))}


def run_policy(ctor, stream, device, image_tf, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = ctor().to(device)
    opt = base.build_optim(model)
    ref_rows = stream["ref"]
    ref_b, y_ref = base.extract_blocks(model, ref_rows, device, image_tf)
    ema = {m: 1.0 / len(MODS) for m in MODS}
    traj = []

    for w in stream["windows"]:
        adapt_rows, hold_rows = split_rows(w["rows"], seed=SEED + 7 * w["t"])
        cur_b, y_cur = base.extract_blocks(model, adapt_rows, device, image_tf)
        msg = base.compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        raw = msg.alpha[policy]
        for m in MODS:
            ema[m] = EMA * ema[m] + (1.0 - EMA) * raw[m]
        s = sum(ema.values())
        alpha = {m: ema[m] / s for m in MODS}

        down_pre = {
            "hold": eval_downstream(model, hold_rows, device, image_tf),
            "ref": eval_downstream(model, ref_rows[:120], device, image_tf),
        }
        stats = base.train_window(
            model,
            opt,
            adapt_rows,
            device,
            image_tf,
            alpha=alpha,
            hard_gate=(policy != "B1"),
        )
        down_post = {
            "hold": eval_downstream(model, hold_rows, device, image_tf),
            "ref": eval_downstream(model, ref_rows[:120], device, image_tf),
        }

        ref_b = {
            m: np.vstack(
                [ref_b[m][-(base.N_REF // 2) :], cur_b[m][: base.N_REF // 2]]
            )
            for m in MODS
        }
        y_ref = np.concatenate(
            [y_ref[-(base.N_REF // 2) :], y_cur[: base.N_REF // 2]]
        )

        hold_auc_delta = float("nan")
        if np.isfinite(down_pre["hold"]["auc"]) and np.isfinite(down_post["hold"]["auc"]):
            hold_auc_delta = down_post["hold"]["auc"] - down_pre["hold"]["auc"]

        row = {
            "t": w["t"],
            "category": w["category"],
            "n_adapt": len(adapt_rows),
            "n_hold": len(hold_rows),
            "upstream": {
                "auc": msg.auc,
                "vimp": msg.vimp,
                "po": msg.po,
                "g": msg.g,
                "alpha_raw": raw,
                "alpha_ema": alpha,
                "lr_mult": stats["lr_mult"],
                "active": {m: bool(stats["active"][m]) for m in MODS},
                "rel_adapt_flops": stats["rel_adapt_flops"],
            },
            "downstream": {
                "train_loss": stats["mean_loss"],
                "hold_acc_pre": down_pre["hold"]["acc"],
                "hold_acc_post": down_post["hold"]["acc"],
                "hold_acc_delta": down_post["hold"]["acc"] - down_pre["hold"]["acc"],
                "hold_auc_pre": down_pre["hold"]["auc"],
                "hold_auc_post": down_post["hold"]["auc"],
                "hold_auc_delta": hold_auc_delta,
                "ref_acc_pre": down_pre["ref"]["acc"],
                "ref_acc_post": down_post["ref"]["acc"],
                "ref_acc_delta": down_post["ref"]["acc"] - down_pre["ref"]["acc"],
            },
        }
        traj.append(row)
        u, d = row["upstream"], row["downstream"]
        print(
            f"[{policy}] t={w['t']} {w['category'][:22]:<22} "
            f"α_img={u['alpha_ema']['image']:.3f} lr×_img={u['lr_mult']['image']:.2f} "
            f"flops={u['rel_adapt_flops']:.2f} | "
            f"holdΔacc={d['hold_acc_delta']:+.3f} holdΔauc={d['hold_auc_delta']:+.3f} "
            f"refΔacc={d['ref_acc_delta']:+.3f}"
        )
    return traj


def summarize(results: dict) -> dict:
    out = {}
    for pol, traj in results.items():
        a = np.array([r["upstream"]["alpha_ema"]["image"] for r in traj], float)
        b = np.array([r["downstream"]["hold_acc_delta"] for r in traj], float)
        corr = float("nan")
        if len(a) >= 3 and np.std(a) > 1e-8 and np.std(b) > 1e-8:
            corr = float(np.corrcoef(a, b)[0, 1])
        out[pol] = {
            "mean_rel_flops": float(
                np.mean([r["upstream"]["rel_adapt_flops"] for r in traj])
            ),
            "mean_train_loss": float(
                np.mean([r["downstream"]["train_loss"] for r in traj])
            ),
            "mean_hold_acc_delta": float(
                np.nanmean([r["downstream"]["hold_acc_delta"] for r in traj])
            ),
            "mean_hold_auc_delta": float(
                np.nanmean([r["downstream"]["hold_auc_delta"] for r in traj])
            ),
            "mean_ref_acc_delta": float(
                np.nanmean([r["downstream"]["ref_acc_delta"] for r in traj])
            ),
            "mean_alpha_image": float(a.mean()),
            "mean_lr_image": float(
                np.mean([r["upstream"]["lr_mult"]["image"] for r in traj])
            ),
            "corr_alpha_img_hold_dacc": corr,
        }
    return out


def plot_updown(results, stream, path: Path):
    fig = plt.figure(figsize=(13.2, 9.0), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.28)
    fig.suptitle(
        "Amazon AGOD — Upstream MSG/LR vs Downstream Acc/AUC",
        fontsize=14,
        fontweight="bold",
    )
    b3 = results["B3"]
    ts = [r["t"] + 1 for r in b3]

    ax = fig.add_subplot(gs[0, 0])
    for m, c in zip(MODS, ["#2B6CB0", "#C05621"]):
        ax.plot(
            ts,
            [r["upstream"]["alpha_ema"][m] for r in b3],
            "o-",
            color=c,
            label=rf"$\alpha_{{{m}}}$",
        )
        ax.plot(
            ts,
            [r["upstream"]["g"][m] for r in b3],
            "x--",
            color=c,
            alpha=0.55,
            label=rf"$g_{{{m}}}$",
        )
    ax.set_ylim(0, 1.05)
    ax.set_title("Upstream: MSG → routing mass (B3)")
    ax.set_xlabel("window t")
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = fig.add_subplot(gs[0, 1])
    for m, c in zip(MODS + ["shared"], ["#2B6CB0", "#C05621", "#718096"]):
        ax.plot(
            ts, [r["upstream"]["lr_mult"][m] for r in b3], "o-", color=c, label=m
        )
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("Action: per-modality LR multipliers (B3)")
    ax.set_xlabel("window t")
    ax.set_ylabel(r"$\times$ lr$_0$")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 0])
    for pol, c, mk in zip(
        ["B1", "B2", "B3"], ["#718096", "#DD6B20", "#C53030"], ["o", "s", "D"]
    ):
        ax.plot(
            [r["t"] + 1 for r in results[pol]],
            [r["downstream"]["hold_acc_delta"] for r in results[pol]],
            marker=mk,
            color=c,
            label=f"{pol} hold ΔAcc",
        )
    ax.axhline(0.0, color="#999", ls="--", lw=0.8)
    ax.set_title("Downstream: held-out Acc change after update")
    ax.set_xlabel("window t")
    ax.set_ylabel("Δ Acc (post − pre)")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    width = 0.25
    xs = np.arange(3)
    flops = [
        np.mean([r["upstream"]["rel_adapt_flops"] for r in results[p]])
        for p in ["B1", "B2", "B3"]
    ]
    dacc = [
        np.nanmean([r["downstream"]["hold_acc_delta"] for r in results[p]])
        for p in ["B1", "B2", "B3"]
    ]
    ax.bar(xs - width / 2, flops, width, color="#4A5568", label="rel adapt FLOPs")
    ax.bar(xs + width / 2, dacc, width, color="#C53030", label="mean hold ΔAcc")
    ax.set_xticks(xs)
    ax.set_xticklabels(["B1", "B2", "B3"])
    ax.set_title("Upstream cost vs downstream gain")
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0.0, color="#999", ls="--", lw=0.7)

    cats = " → ".join(w["category"][:16] for w in stream["windows"])
    fig.text(
        0.5,
        0.012,
        f"ref={stream['ref_cat'][:22]} ‖ {cats}  |  "
        "Upstream drives LR only; sample/patch grains discussed in docs, not actuated.",
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
    if not shards:
        raise SystemExit(f"No shards in {base.SHARD_DIR}")
    t0 = time.time()
    samples = base.load_shards(shards)
    print(f"samples={len(samples)} in {time.time() - t0:.1f}s")

    image_tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
        ]
    )
    stream = base.make_stream(samples)
    print(
        f"ref={stream['ref_cat']} windows={[w['category'] for w in stream['windows']]}"
    )

    results = {}
    for pol in ["B1", "B2", "B3"]:
        print(f"\n===== policy {pol} =====")
        results[pol] = run_policy(base.AmazonAGOD, stream, device, image_tf, pol)

    summary = summarize(results)
    payload = {
        "dataset": "jingxiang11111/amazon_reviews_for_rec",
        "ref_category": stream["ref_cat"],
        "protocol": {
            "action": "modality LR + grad gate from MSG only",
            "upstream": [
                "AUC_m",
                "VIMP_m",
                "PO_m",
                "g_m",
                "alpha_m",
                "LR_mult",
                "adapt_flops",
            ],
            "downstream": [
                "holdout Acc/AUC pre/post",
                "ref Acc pre/post",
                "train loss",
            ],
            "not_implemented_actions": ["sample reweight", "patch/token mask"],
        },
        "summary": summary,
        "trajectory": results,
    }
    json_path = OUT / "agod_amazon_updownstream.json"
    json_path.write_text(json.dumps(payload, indent=2, default=float))

    dash = OUT / "AGOD_Amazon_Upstream_Downstream_Dashboard.png"
    plot_updown(results, stream, dash)

    logic_path = DOCS / "AGOD_attribution_hierarchy_modality_sample_token.md"
    logic_path.write_text(LOGIC_MD)
    (OUT / "AGOD_attribution_hierarchy_modality_sample_token.md").write_text(LOGIC_MD)

    lines = [
        r"\begin{tabular}{lccccc}",
        r"\toprule",
        r"Policy & Rel.\ FLOPs & Hold $\Delta$Acc & Hold $\Delta$AUC & Ref $\Delta$Acc & $\alpha_{\mathrm{img}}$ \\",
        r"\midrule",
    ]
    for pol in ["B1", "B2", "B3"]:
        s = summary[pol]
        lines.append(
            f"{pol} & {s['mean_rel_flops']:.3f} & {s['mean_hold_acc_delta']:+.3f} & "
            f"{s['mean_hold_auc_delta']:+.3f} & {s['mean_ref_acc_delta']:+.3f} & "
            f"{s['mean_alpha_image']:.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    tex = "\n".join(lines)
    (OUT / "AGOD_amazon_updownstream_tables_only.tex").write_text(tex)
    (DOCS / "AGOD_amazon_updownstream_tables_only.tex").write_text(tex)

    for p in [dash, json_path, logic_path]:
        (ART / Path(p).name).write_bytes(Path(p).read_bytes())

    print("\n=== UP/DOWN SUMMARY ===")
    print(json.dumps(summary, indent=2))
    print("dashboard", dash)
    print("logic", logic_path)


if __name__ == "__main__":
    main()
