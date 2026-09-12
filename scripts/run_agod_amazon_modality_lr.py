#!/usr/bin/env python3
"""AGOD on Amazon Reviews: modality-level MSG → per-modality learning rates.

Dataset shards (local): data/amazon_reviews/shards/*.tar.gz
  Source: jingxiang11111/amazon_reviews_for_rec

Sample fields:
  *.item.json / *.user.json  → text
  *.patch.bin + *.misc.json → 196×3×16×16 uint8 patches → 224² RGB
  *.label.json              → label_good / label_best

Core control loop
-----------------
  MSG_m = Normalize(AUC_m * meanVIMP_m + γ * PO_m)
  α = Softmax(MSG / τ)
  LR_m ← lr0 * (β + (1-β) * α_m * |M|)
  hard-gate: if α_m < θ, zero grads of modality-m projection
             (adaptation FLOPs ↓; NOT inference latency)

  python3 scripts/run_agod_amazon_modality_lr.py
"""
from __future__ import annotations

import json
import re
import tarfile
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from PIL import Image
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset
from torchvision import transforms
from torchvision.models import ResNet18_Weights, resnet18

ROOT = Path(__file__).resolve().parents[1]
SHARD_DIR = ROOT / "data" / "amazon_reviews" / "shards"
OUT = ROOT / "results" / "agod_amazon"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_amazon")

SEED = 2026
FUSION = 128
TEXT_DIM = 512
N_REF = 180
N_CUR = 120
T_STEPS = 6
GAMMA = 1.0
TAU = 0.28
GATE_TH = 0.40  # with |M|=2, gate minority modality when α < 0.40
LR0 = 3e-4
EMA = 0.55
RF_TREES = 80
PO_TREES = 40
STEPS_PER_WIN = 8
BATCH = 16
MODS = ["text", "image"]


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def parse_category(item_text: str) -> str:
    m = re.search(r"Item category:\s*([^.]+)", item_text or "")
    return m.group(1).strip() if m else "Unknown"


def patches_to_image(raw: bytes, shape) -> Image.Image:
    arr = np.frombuffer(raw, dtype=np.uint8).reshape(tuple(shape))
    n, c, h, w = arr.shape
    g = int(np.sqrt(n))
    if g * g != n:
        raise ValueError(f"non-square patch grid n={n}")
    canvas = np.zeros((g * h, g * w, c), dtype=np.uint8)
    for i in range(n):
        r, c0 = divmod(i, g)
        canvas[r * h : (r + 1) * h, c0 * w : (c0 + 1) * w] = arr[i].transpose(1, 2, 0)
    return Image.fromarray(canvas)


def load_shards(paths: list[Path]) -> list[dict]:
    need = {"item.json", "user.json", "patch.bin", "misc.json", "label.json"}
    samples: list[dict] = []
    for path in paths:
        with tarfile.open(path, "r:gz") as tf:
            by_key: dict[str, dict] = defaultdict(dict)
            for mem in tf.getmembers():
                if not mem.isfile() or "." not in mem.name:
                    continue
                key, ext = mem.name.split(".", 1)
                by_key[key][ext] = tf.extractfile(mem).read()
            for key, parts in by_key.items():
                if not need <= set(parts):
                    continue
                misc = json.loads(parts["misc.json"])
                if not misc.get("has_image", 0):
                    continue
                label = json.loads(parts["label.json"])
                item = parts["item.json"].decode("utf-8", "ignore")
                user = parts["user.json"].decode("utf-8", "ignore")
                samples.append(
                    {
                        "text": f"{item}\n{user}",
                        "image": patches_to_image(parts["patch.bin"], misc["shape"]),
                        "y": int(label.get("label_good", 0)),
                        "category": parse_category(item),
                        "key": key,
                    }
                )
    return samples


class AmazonDS(Dataset):
    def __init__(self, rows, image_tf):
        self.rows = rows
        self.image_tf = image_tf

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, i):
        r = self.rows[i]
        return self.image_tf(r["image"]), r["text"], torch.tensor(r["y"], dtype=torch.long)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class AmazonAGOD(nn.Module):
    def __init__(self):
        super().__init__()
        self.hasher = HashingVectorizer(
            n_features=TEXT_DIM, alternate_sign=False, norm="l2", ngram_range=(1, 2)
        )
        bb = resnet18(weights=ResNet18_Weights.DEFAULT)
        bb.fc = nn.Identity()
        for p in bb.parameters():
            p.requires_grad = False
        self.image_encoder = bb
        self.text_proj = nn.Sequential(
            nn.Linear(TEXT_DIM, FUSION),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(FUSION, FUSION),
        )
        self.image_proj = nn.Sequential(
            nn.Linear(512, FUSION),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(FUSION, FUSION),
        )
        self.gate = nn.Linear(FUSION * 2, 2)
        self.head = nn.Sequential(
            nn.Linear(FUSION, 64),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(64, 2),
        )

    def encode_text(self, texts, device):
        x = self.hasher.transform(list(texts)).toarray().astype(np.float32)
        return torch.from_numpy(x).to(device)

    def forward(self, images, texts, return_mods: bool = False):
        with torch.no_grad():
            img_raw = self.image_encoder(images)
        txt_raw = self.encode_text(texts, images.device)
        t_feat = self.text_proj(txt_raw)
        i_feat = self.image_proj(img_raw)
        gates = torch.softmax(self.gate(torch.cat([t_feat, i_feat], dim=1)), dim=1)
        fused = gates[:, :1] * t_feat + gates[:, 1:] * i_feat
        logits = self.head(fused)
        if return_mods:
            return logits, {"text": t_feat, "image": i_feat, "gates": gates}
        return logits

    def param_groups(self):
        return {
            "text": list(self.text_proj.parameters()),
            "image": list(self.image_proj.parameters()),
            "shared": list(self.gate.parameters()) + list(self.head.parameters()),
        }


# ---------------------------------------------------------------------------
# MSG attribution
# ---------------------------------------------------------------------------


def rf_auc_vimp(X0, X1, *, seed: int):
    X = np.vstack([X0, X1])
    W = np.array([0] * len(X0) + [1] * len(X1))
    if len(np.unique(W)) < 2 or len(X) < 40:
        return 0.5, 0.0
    clf = RandomForestClassifier(
        n_estimators=RF_TREES,
        max_depth=max(3, int(np.sqrt(X.shape[1]))),
        min_samples_leaf=2,
        n_jobs=-1,
        random_state=seed,
        oob_score=True,
    )
    clf.fit(X, W)
    vimp = float(clf.feature_importances_.mean())
    try:
        Xtr, Xte, Wtr, Wte = train_test_split(
            X, W, test_size=0.3, random_state=seed, stratify=W
        )
        clf2 = RandomForestClassifier(
            n_estimators=60,
            max_depth=max(3, int(np.sqrt(X.shape[1]))),
            min_samples_leaf=2,
            n_jobs=-1,
            random_state=seed + 1,
        )
        clf2.fit(Xtr, Wtr)
        auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    except Exception:
        auc = float(getattr(clf, "oob_score_", 0.5))
    return auc, vimp


def po_risk(X0, Y0, X1, Y1, *, seed: int):
    X = np.vstack([X0, X1])
    Y = np.concatenate([Y0, Y1]).astype(float)
    if len(X) < 40:
        return 0.0
    model = RandomForestRegressor(
        n_estimators=PO_TREES,
        max_depth=6,
        min_samples_leaf=2,
        n_jobs=-1,
        random_state=seed,
    )
    model.fit(X, Y)
    pred = model.predict(X)
    r0 = float(np.mean((Y[: len(X0)] - pred[: len(X0)]) ** 2))
    r1 = float(np.mean((Y[len(X0) :] - pred[len(X0) :]) ** 2))
    return abs(r1 - r0)


def _norm(v):
    v = np.maximum(np.asarray(v, float), 0.0)
    s = v.sum()
    return v / s if s > 1e-12 else np.ones_like(v) / len(v)


def _softmax(z, tau):
    z = np.asarray(z, float) / max(tau, 1e-6)
    z = z - z.max()
    e = np.exp(z)
    return e / e.sum()


@dataclass
class MSG:
    auc: dict
    vimp: dict
    po: dict
    g: dict
    alpha: dict


def compute_msg(b0, b1, y0, y1, *, seed: int = SEED) -> MSG:
    auc, vimp, po, raw = {}, {}, {}, {}
    for i, m in enumerate(MODS):
        a, v = rf_auc_vimp(b0[m], b1[m], seed=seed + i)
        p = po_risk(b0[m], y0, b1[m], y1, seed=seed + 17 + i)
        auc[m], vimp[m], po[m] = a, v, p
        raw[m] = a * v + GAMMA * p
    gvec = _norm([raw[m] for m in MODS])
    g = {m: float(gvec[i]) for i, m in enumerate(MODS)}
    alpha = {
        "B1": {m: 1.0 / len(MODS) for m in MODS},
        "B2": {
            m: float(x) for m, x in zip(MODS, _softmax([auc[m] for m in MODS], TAU))
        },
        "B3": {m: float(x) for m, x in zip(MODS, _softmax(gvec, TAU))},
    }
    return MSG(auc=auc, vimp=vimp, po=po, g=g, alpha=alpha)


def alpha_to_lr_mult(alpha: dict, beta: float = 0.15) -> dict:
    inv = float(len(MODS))
    out = {m: float(beta + (1.0 - beta) * alpha[m] * inv) for m in MODS}
    out["shared"] = float(np.mean([out[m] for m in MODS]))
    return out


# ---------------------------------------------------------------------------
# Stream / features
# ---------------------------------------------------------------------------


@torch.no_grad()
def extract_blocks(model: AmazonAGOD, rows, device, image_tf):
    model.eval()
    loader = DataLoader(AmazonDS(rows, image_tf), batch_size=BATCH, shuffle=False)
    texts, images, ys = [], [], []
    for imgs, txts, y in loader:
        imgs = imgs.to(device)
        _, feats = model(imgs, txts, return_mods=True)
        texts.append(feats["text"].cpu().numpy())
        images.append(feats["image"].cpu().numpy())
        ys.append(y.numpy())
    return (
        {"text": np.concatenate(texts), "image": np.concatenate(images)},
        np.concatenate(ys).astype(float),
    )


def make_stream(samples: list[dict]):
    by_cat = defaultdict(list)
    for s in samples:
        by_cat[s["category"]].append(s)
    ref_cat = "Tools & Home Improvement"
    if len(by_cat[ref_cat]) < N_REF:
        ref_cat = max(by_cat, key=lambda k: len(by_cat[k]))
    ref = by_cat[ref_cat][: max(N_REF, min(250, len(by_cat[ref_cat])))]
    others = sorted(
        [(c, rows) for c, rows in by_cat.items() if c != ref_cat and len(rows) >= 40],
        key=lambda x: -len(x[1]),
    )
    windows = []
    for t in range(T_STEPS):
        if not others:
            break
        c, rows = others[t % len(others)]
        n = min(N_CUR, len(rows))
        cur = rows[:n]
        mix_n = int((0.35 * (1 - t / max(T_STEPS - 1, 1))) * n)
        if mix_n > 0:
            cur = ref[:mix_n] + cur[mix_n:]
        windows.append({"t": t, "category": c, "rows": cur})
    return {"ref_cat": ref_cat, "ref": ref, "windows": windows}


# ---------------------------------------------------------------------------
# Train with modality LR
# ---------------------------------------------------------------------------


def build_optim(model: AmazonAGOD):
    g = model.param_groups()
    return torch.optim.AdamW(
        [
            {"params": g["text"], "lr": LR0, "name": "text"},
            {"params": g["image"], "lr": LR0, "name": "image"},
            {"params": g["shared"], "lr": LR0, "name": "shared"},
        ]
    )


def set_lrs(opt, mult: dict):
    for g in opt.param_groups:
        g["lr"] = LR0 * float(mult.get(g.get("name", "shared"), 1.0))


def train_window(model, opt, rows, device, image_tf, *, alpha, hard_gate=True):
    model.train()
    model.image_encoder.eval()
    mult = alpha_to_lr_mult(alpha)
    set_lrs(opt, mult)
    active = {m: (alpha[m] >= GATE_TH) for m in MODS}
    if hard_gate and not any(active.values()):
        mstar = max(MODS, key=lambda m: alpha[m])
        active = {m: m == mstar for m in MODS}

    loader = DataLoader(AmazonDS(rows, image_tf), batch_size=BATCH, shuffle=True)
    crit = nn.CrossEntropyLoss()
    losses = []
    it = iter(loader)
    for _ in range(STEPS_PER_WIN):
        try:
            imgs, txts, y = next(it)
        except StopIteration:
            it = iter(loader)
            imgs, txts, y = next(it)
        imgs, y = imgs.to(device), y.to(device)
        logits, _ = model(imgs, txts, return_mods=True)
        loss = crit(logits, y)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        if hard_gate:
            groups = model.param_groups()
            for m in MODS:
                if not active[m]:
                    for p in groups[m]:
                        p.grad = None
        opt.step()
        losses.append(float(loss.item()))
    return {
        "mean_loss": float(np.mean(losses)),
        "lr_mult": mult,
        "active": active,
        "rel_adapt_flops": float(sum(active.values()) / len(MODS)),
    }


def run_policy(ctor, stream, device, image_tf, policy: str):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    model = ctor().to(device)
    opt = build_optim(model)
    ref_b, y_ref = extract_blocks(model, stream["ref"], device, image_tf)
    ema = {m: 1.0 / len(MODS) for m in MODS}
    traj = []
    for w in stream["windows"]:
        cur_b, y_cur = extract_blocks(model, w["rows"], device, image_tf)
        msg = compute_msg(ref_b, cur_b, y_ref, y_cur, seed=SEED + 10 * w["t"])
        raw = msg.alpha[policy]
        for m in MODS:
            ema[m] = EMA * ema[m] + (1.0 - EMA) * raw[m]
        s = sum(ema.values())
        alpha = {m: ema[m] / s for m in MODS}
        stats = train_window(
            model,
            opt,
            w["rows"],
            device,
            image_tf,
            alpha=alpha,
            hard_gate=(policy != "B1"),
        )
        ref_b = {
            m: np.vstack([ref_b[m][-(N_REF // 2) :], cur_b[m][: N_REF // 2]])
            for m in MODS
        }
        y_ref = np.concatenate([y_ref[-(N_REF // 2) :], y_cur[: N_REF // 2]])
        traj.append(
            {
                "t": w["t"],
                "category": w["category"],
                "auc": msg.auc,
                "vimp": msg.vimp,
                "po": msg.po,
                "g": msg.g,
                "alpha_raw": raw,
                "alpha_ema": alpha,
                **stats,
            }
        )
        print(
            f"[{policy}] t={w['t']} cat={w['category'][:28]:<28} "
            f"α={ {m: round(alpha[m], 3) for m in MODS} } "
            f"lr×={ {k: round(v, 3) for k, v in stats['lr_mult'].items()} } "
            f"flops={stats['rel_adapt_flops']:.2f} loss={stats['mean_loss']:.3f}"
        )
    return traj


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def plot_dashboard(results, stream, path: Path):
    fig = plt.figure(figsize=(12.8, 8.4), facecolor="#f7f5f1")
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)
    fig.suptitle(
        "Amazon Reviews AGOD — Modality MSG → Per-modality LR",
        fontsize=14,
        fontweight="bold",
    )
    traj = results["B3"]
    ts = [r["t"] + 1 for r in traj]

    ax = fig.add_subplot(gs[0, 0])
    for m, c in zip(MODS, ["#2B6CB0", "#C05621"]):
        ax.plot(ts, [r["alpha_ema"][m] for r in traj], "o-", color=c, label=rf"$\alpha_{{{m}}}$")
        ax.plot(
            ts, [r["g"][m] for r in traj], "x--", color=c, alpha=0.55, label=rf"$g_{{{m}}}$"
        )
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("online window t")
    ax.set_ylabel("mass")
    ax.set_title(f"MSG & EMA routing (ref={stream['ref_cat'][:28]})")
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = fig.add_subplot(gs[0, 1])
    for m, c in zip(MODS + ["shared"], ["#2B6CB0", "#C05621", "#718096"]):
        ax.plot(ts, [r["lr_mult"][m] for r in traj], "o-", color=c, label=m)
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_xlabel("online window t")
    ax.set_ylabel(r"LR multiplier ($\times$ lr$_0$)")
    ax.set_title("Per-modality learning-rate schedule (B3)")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 0])
    for pol, c, mk in zip(
        ["B1", "B2", "B3"], ["#718096", "#DD6B20", "#C53030"], ["o", "s", "D"]
    ):
        ax.plot(
            [r["t"] + 1 for r in results[pol]],
            [r["rel_adapt_flops"] for r in results[pol]],
            marker=mk,
            color=c,
            label=pol,
        )
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("online window t")
    ax.set_ylabel("relative adaptation FLOPs")
    ax.set_title("Hard-gate update cost (not inference latency)")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    for pol, c, mk in zip(
        ["B1", "B2", "B3"], ["#718096", "#DD6B20", "#C53030"], ["o", "s", "D"]
    ):
        ax.plot(
            [r["t"] + 1 for r in results[pol]],
            [r["mean_loss"] for r in results[pol]],
            marker=mk,
            color=c,
            label=pol,
        )
    ax.set_xlabel("online window t")
    ax.set_ylabel("window train loss")
    ax.set_title("B1 static / B2 AUC-only / B3 AGOD-LR")
    ax.legend(frameon=False, fontsize=8)

    cats = " → ".join(w["category"][:18] for w in stream["windows"])
    fig.text(
        0.5,
        0.015,
        f"Stream: {stream['ref_cat'][:24]} ‖ {cats}  |  "
        r"Control: $g_m\to\alpha_m\to$LR$_m$; skip $\partial L/\partial\theta_m$ if $\alpha_m<\theta$",
        ha="center",
        fontsize=8.5,
        color="#333",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def to_latex(summary: dict) -> str:
    lines = [
        r"\begin{tabular}{llccc}",
        r"\toprule",
        r"Dataset & Policy & Rel.\ adapt.\ FLOPs & Mean loss & Mean $\alpha_{\mathrm{text}}$ \\",
        r"\midrule",
    ]
    for pol in ["B1", "B2", "B3"]:
        s = summary[pol]
        lines.append(
            f"Amazon Reviews & {pol} & {s['mean_rel_flops']:.3f} & "
            f"{s['mean_loss']:.3f} & {s['mean_alpha_text']:.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}")
    shards = sorted(SHARD_DIR.glob("*.tar.gz"))
    if not shards:
        raise SystemExit(f"No shards in {SHARD_DIR}")
    print(f"loading {len(shards)} shards…")
    t0 = time.time()
    samples = load_shards(shards)
    print(f"samples={len(samples)} in {time.time() - t0:.1f}s")
    print("categories:", Counter(s["category"] for s in samples).most_common(8))

    image_tf = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
        ]
    )
    stream = make_stream(samples)
    print(
        f"ref={stream['ref_cat']} n_ref={len(stream['ref'])} "
        f"windows={[w['category'] for w in stream['windows']]}"
    )

    results = {}
    for pol in ["B1", "B2", "B3"]:
        print(f"\n===== policy {pol} =====")
        results[pol] = run_policy(AmazonAGOD, stream, device, image_tf, pol)

    summary = {}
    for pol, traj in results.items():
        summary[pol] = {
            "mean_rel_flops": float(np.mean([r["rel_adapt_flops"] for r in traj])),
            "mean_loss": float(np.mean([r["mean_loss"] for r in traj])),
            "mean_alpha_text": float(np.mean([r["alpha_ema"]["text"] for r in traj])),
            "mean_alpha_image": float(np.mean([r["alpha_ema"]["image"] for r in traj])),
            "mean_lr_text": float(np.mean([r["lr_mult"]["text"] for r in traj])),
            "mean_lr_image": float(np.mean([r["lr_mult"]["image"] for r in traj])),
        }

    payload = {
        "dataset": "jingxiang11111/amazon_reviews_for_rec",
        "ref_category": stream["ref_cat"],
        "windows": [
            {"t": w["t"], "category": w["category"], "n": len(w["rows"])}
            for w in stream["windows"]
        ],
        "config": {
            "gamma": GAMMA,
            "tau": TAU,
            "gate_th": GATE_TH,
            "lr0": LR0,
            "ema": EMA,
            "note": "adaptation FLOPs via gated modality updates; NOT inference latency",
        },
        "summary": summary,
        "trajectory": results,
    }
    (OUT / "agod_amazon_modality_lr.json").write_text(
        json.dumps(payload, indent=2, default=float)
    )
    tex = to_latex(summary)
    (OUT / "AGOD_amazon_modality_lr_tables_only.tex").write_text(tex)
    (DOCS / "AGOD_amazon_modality_lr_tables_only.tex").write_text(tex)

    dash = OUT / "AGOD_Amazon_Modality_LR_Dashboard.png"
    plot_dashboard(results, stream, dash)
    for p in [dash, OUT / "agod_amazon_modality_lr.json"]:
        (ART / p.name).write_bytes(p.read_bytes())

    (OUT / "README.md").write_text(
        "# AGOD Amazon — modality MSG → per-modality LR\n\n"
        "Not inference latency. Online adaptation: attribute text vs image "
        "distribution-shift contribution, rescale modality projection LRs, "
        "hard-gate inactive modality grads.\n\n"
        f"Ref category: `{stream['ref_cat']}`\n\n"
        + "\n".join(
            f"- {pol}: flops={summary[pol]['mean_rel_flops']:.3f}, "
            f"loss={summary[pol]['mean_loss']:.3f}, "
            f"α_text={summary[pol]['mean_alpha_text']:.3f}, "
            f"lr×_text/image={summary[pol]['mean_lr_text']:.2f}/"
            f"{summary[pol]['mean_lr_image']:.2f}"
            for pol in ["B1", "B2", "B3"]
        )
        + "\n"
    )
    print("\n=== SUMMARY ===")
    print(json.dumps(summary, indent=2))
    print("dashboard", dash)


if __name__ == "__main__":
    main()
