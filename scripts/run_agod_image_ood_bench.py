#!/usr/bin/env python3
"""Image-OOD: PO-risk μ vs classical embedding OOD scores (no DRE).

  PYTHONPATH=. python3 scripts/run_agod_image_ood_bench.py \\
    --packs food101_vit coco_outdoor_indoor coco_time_order indiana_cxr fashion_iq
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from agod.hard_rank_metrics import hard_rank_metrics
from agod.image_ood import (
    eval_scores_on_split,
    fit_mahalanobis,
    fit_po_classifier,
    fit_po_regressor,
    list_available_packs,
    load_image_ood_pack,
    score_centroid,
    score_knn,
    score_mahalanobis,
    score_po_energy,
    score_po_msp,
    score_po_nll,
    score_po_resid,
)

METHODS = ("po_msp", "po_energy", "po_nll", "maha", "knn", "centroid")
SCORE_COLORS = {
    "po_msp": "#5E81AC",
    "po_energy": "#EBCB8B",
    "po_nll": "#88C0D0",
    "maha": "#B48EAD",
    "knn": "#A3BE8C",
    "centroid": "#D08770",
}


def run_pack(root: Path, pack: str, *, pca_d: int, seed: int, extract_max_n: int) -> dict:
    kwargs = {"pca_d": pca_d, "seed": seed}
    if pack == "food101_vit":
        kwargs["extract_max_n"] = extract_max_n
    split = load_image_ood_pack(root, pack, **kwargs)

    clf = fit_po_classifier(split.X_id_train, split.y_id_train, seed=seed)
    reg = fit_po_regressor(split.X_id_train, split.y_id_train, seed=seed)
    maha = fit_mahalanobis(split.X_id_train, split.y_id_train)

    X_eval = np.vstack([split.X_id_test, split.X_ood])
    y_eval = np.concatenate([split.y_id_test, split.y_ood])

    scores = {
        "po_msp": np.concatenate(
            [score_po_msp(clf, split.X_id_test), score_po_msp(clf, split.X_ood)]
        ),
        "po_energy": np.concatenate(
            [score_po_energy(clf, split.X_id_test), score_po_energy(clf, split.X_ood)]
        ),
        "po_nll": np.concatenate(
            [
                score_po_nll(clf, split.X_id_test, split.y_id_test),
                score_po_nll(clf, split.X_ood, split.y_ood),
            ]
        ),
        "po_resid": np.concatenate(
            [
                score_po_resid(reg, split.X_id_test, split.y_id_test),
                score_po_resid(reg, split.X_ood, split.y_ood),
            ]
        ),
        "maha": np.concatenate(
            [
                score_mahalanobis(maha, split.X_id_test),
                score_mahalanobis(maha, split.X_ood),
            ]
        ),
        "knn": score_knn(split.X_id_train, X_eval, k=5),
        "centroid": score_centroid(split.X_id_train, X_eval, split.y_id_train),
    }

    class_holdout = split.meta.get("protocol") == "class_holdout"
    metrics = eval_scores_on_split(split, scores)
    if class_holdout:
        metrics["po_resid"] = {
            "auroc": float("nan"),
            "aupr": float("nan"),
            "fpr95": float("nan"),
            "score_mean_id": metrics["po_resid"]["score_mean_id"],
            "score_mean_ood": float("nan"),
            "note": "undefined under class-holdout",
        }

    if class_holdout:
        truth = np.concatenate(
            [np.zeros(len(split.X_id_test), float), np.ones(len(split.X_ood), float)]
        )
        rank_names = ("po_msp", "po_energy", "maha", "knn", "centroid")
    else:
        truth = score_po_resid(reg, X_eval, y_eval)
        rank_names = METHODS
    hard = {name: hard_rank_metrics(scores[name], truth) for name in rank_names}

    return {
        "pack": pack,
        "id_domain": split.id_domain,
        "ood_domain": split.ood_domain,
        "meta": split.meta,
        "metrics": metrics,
        "hard_rank": {
            k: {
                "spearman": v.get("spearman"),
                "precision_at_k": v.get("precision_at_k"),
                "lift_at_k": v.get("lift_at_k"),
                "ndcg_at_k": v.get("ndcg_at_k"),
                "auroc_topk": v.get("auroc_topk"),
            }
            for k, v in hard.items()
        },
    }


def _bar(all_res: dict, out_dir: Path, key: str, ylabel: str, title: str, fname: str, ylim=None):
    packs = list(all_res.keys())
    fig, ax = plt.subplots(figsize=(max(9, 1.8 * len(packs)), 4.8))
    x = np.arange(len(packs))
    w = 0.13
    mid = (len(METHODS) - 1) / 2.0
    for i, m in enumerate(METHODS):
        vals = [all_res[p]["metrics"][m][key] for p in packs]
        ax.bar(x + (i - mid) * w, vals, w, label=m, color=SCORE_COLORS[m])
    if key == "auroc":
        ax.axhline(0.5, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(packs, rotation=15, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if ylim:
        ax.set_ylim(*ylim)
    ax.legend(ncol=3, fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / fname
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def report(all_res: dict) -> str:
    def _fmt(v):
        if v is None or (isinstance(v, float) and (v != v)):
            return "—"
        return f"{v:.3f}"

    lines = [
        "# Image-OOD: PO-risk μ vs classical embedding scores (no DRE)",
        "",
        "Freeze backbone → embedding → **pseudo-outcome μ** (PO-risk) vs Mahalanobis / kNN / centroid.",
        "",
        "- `food101_vit`: frozen ViT + **class-holdout** (far-OOD) — primary.",
        "- CLIP packs: cached img embeddings + domain shift (near-OOD).",
        "",
        "### Scores",
        "",
        "| score | definition |",
        "|---|---|",
        "| `po_msp` | `1 − max_c p_μ(c\\|x)` — PO-risk (label-free) |",
        "| `po_energy` | Shannon entropy of μ — PO-risk (label-free) |",
        "| `po_nll` | `1 − p_μ(y\\|x)` (MSP fallback if y unseen) |",
        "| `maha` | min class-conditional Mahalanobis (Lee et al.) |",
        "| `knn` | mean L2 to 5-NN in ID-train (Sun et al.) |",
        "| `centroid` | L2 to nearest ID class centroid |",
        "",
        "## AUROC (ID vs OOD)",
        "",
        "| pack | shift | po_msp | po_energy | po_nll | maha | knn | centroid | best |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    wins = {m: 0 for m in METHODS}
    for pack, blob in all_res.items():
        m = blob["metrics"]
        best = max(
            METHODS,
            key=lambda k: (m[k]["auroc"] if m[k]["auroc"] == m[k]["auroc"] else -1.0),
        )
        if m[best]["auroc"] == m[best]["auroc"]:
            wins[best] += 1
        shift = f"{blob['id_domain']}→{blob['ood_domain']}"
        lines.append(
            f"| `{pack}` | {shift} | {_fmt(m['po_msp']['auroc'])} | {_fmt(m['po_energy']['auroc'])} | "
            f"{_fmt(m['po_nll']['auroc'])} | {_fmt(m['maha']['auroc'])} | {_fmt(m['knn']['auroc'])} | "
            f"{_fmt(m['centroid']['auroc'])} | `{best}` |"
        )
    lines += [
        "",
        "**AUROC wins:** " + ", ".join(f"`{k}`={v}" for k, v in wins.items()),
        "",
        "## FPR95 (↓ better)",
        "",
        "| pack | po_msp | po_energy | po_nll | maha | knn | centroid |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for pack, blob in all_res.items():
        m = blob["metrics"]
        lines.append(
            f"| `{pack}` | {_fmt(m['po_msp']['fpr95'])} | {_fmt(m['po_energy']['fpr95'])} | "
            f"{_fmt(m['po_nll']['fpr95'])} | {_fmt(m['maha']['fpr95'])} | {_fmt(m['knn']['fpr95'])} | "
            f"{_fmt(m['centroid']['fpr95'])} |"
        )
    lines += [
        "",
        "### Takeaway",
        "",
        "Compare **PO-risk on a frozen embedding** to standard embedding OOD detectors.",
        "No DRE — domain classifiers are not part of this comparison.",
        "",
        "See `docs/agod/AGOD_image_ood_bench.md`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--packs", nargs="+", default=None)
    ap.add_argument("--pca-d", type=int, default=64)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--extract-max-n", type=int, default=8000)
    ap.add_argument("--out", type=Path, default=Path("results/agod_image_ood"))
    args = ap.parse_args()

    available = list_available_packs(args.root)
    packs = args.packs or available
    packs = [p for p in packs if p in available or p == "food101_vit"]
    if not packs:
        raise SystemExit(f"no packs available under {args.root}; found {available}")

    args.out.mkdir(parents=True, exist_ok=True)
    all_res: Dict[str, dict] = {}
    for pack in packs:
        print(f"=== {pack} ===", flush=True)
        try:
            blob = run_pack(
                args.root,
                pack,
                pca_d=args.pca_d,
                seed=args.seed,
                extract_max_n=args.extract_max_n,
            )
        except Exception as e:
            print(f"[skip] {pack}: {e}", flush=True)
            continue
        all_res[pack] = blob
        m = blob["metrics"]
        print(
            f"  AUROC po_msp={m['po_msp']['auroc']:.3f} po_energy={m['po_energy']['auroc']:.3f} "
            f"maha={m['maha']['auroc']:.3f} knn={m['knn']['auroc']:.3f} "
            f"centroid={m['centroid']['auroc']:.3f}",
            flush=True,
        )

    payload = {
        "pca_d": args.pca_d,
        "seed": args.seed,
        "extract_max_n": args.extract_max_n,
        "methods": list(METHODS),
        "packs": list(all_res.keys()),
        "results": all_res,
        "note": "PO-risk vs maha/knn/centroid — no DRE",
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_res)
    (args.out / "IMAGE_OOD_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_image_ood_bench.md").write_text(md, encoding="utf-8")
    plots = []
    if all_res:
        plots.append(
            _bar(
                all_res,
                args.out,
                "auroc",
                "AUROC",
                "Image-OOD AUROC: PO-risk vs maha / knn / centroid (no DRE)",
                "image_ood_auroc.png",
                ylim=(0.35, 1.02),
            )
        )
        plots.append(
            _bar(
                all_res,
                args.out,
                "fpr95",
                "FPR95 (↓)",
                "Image-OOD FPR95 (no DRE)",
                "image_ood_fpr95.png",
            )
        )
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
