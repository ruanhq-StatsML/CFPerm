#!/usr/bin/env python3
"""Image-OOD benchmark redesign: ViT embedding + pseudo-outcome PO-risk.

Protocol
--------
1. Freeze ViT (or reuse cached CLIP img feats as embeddings).
2. Fit pseudo-outcome μ on **ID train** embeddings only.
3. Score PO-risk on ID-holdout ∪ OOD:
     po_msp / po_energy   — label-free (far-OOD / class-holdout)
     po_nll / po_resid    — need y (near-OOD / domain shift)
     dre                  — domain classifier baseline
4. Primary metrics: AUROC / FPR95 / AUPR; secondary: hard-rank vs residual.

  PYTHONPATH=. python3 scripts/run_agod_image_ood_bench.py \\
    --packs food101_vit coco_outdoor_indoor coco_time_order indiana_cxr fashion_iq
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from agod.hard_rank_metrics import hard_rank_metrics
from agod.image_ood import (
    eval_scores_on_split,
    fit_po_classifier,
    fit_po_regressor,
    list_available_packs,
    load_image_ood_pack,
    score_dre_domain,
    score_energy,
    score_po_msp,
    score_po_nll,
    score_po_resid,
)

SCORE_COLORS = {
    "po_msp": "#5E81AC",
    "po_nll": "#88C0D0",
    "po_resid": "#A3BE8C",
    "po_energy": "#EBCB8B",
    "dre": "#BF616A",
}
METHODS = ("po_msp", "po_nll", "po_resid", "po_energy", "dre")


def run_pack(root: Path, pack: str, *, pca_d: int, seed: int, extract_max_n: int) -> dict:
    kwargs = {"pca_d": pca_d, "seed": seed}
    if pack == "food101_vit":
        kwargs["extract_max_n"] = extract_max_n
    split = load_image_ood_pack(root, pack, **kwargs)

    clf = fit_po_classifier(split.X_id_train, split.y_id_train, seed=seed)
    reg = fit_po_regressor(split.X_id_train, split.y_id_train, seed=seed)

    X_eval = np.vstack([split.X_id_test, split.X_ood])
    y_eval = np.concatenate([split.y_id_test, split.y_ood])

    scores = {
        "po_msp": np.concatenate(
            [score_po_msp(clf, split.X_id_test), score_po_msp(clf, split.X_ood)]
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
        "po_energy": np.concatenate(
            [score_energy(clf, split.X_id_test), score_energy(clf, split.X_ood)]
        ),
        "dre": score_dre_domain(split.X_id_train, X_eval, seed=seed),
    }

    metrics = eval_scores_on_split(split, scores)

    # Hard-rank diagnostic:
    # - class-holdout: truth = binary OOD indicator (far-OOD ranking)
    # - domain-shift: truth = residual hardness (resid self-rank is tautological → still shown)
    if split.meta.get("protocol") == "class_holdout":
        truth = np.concatenate(
            [np.zeros(len(split.X_id_test), float), np.ones(len(split.X_ood), float)]
        )
        rank_names = ("po_msp", "po_nll", "po_resid", "po_energy")
    else:
        truth = score_po_resid(reg, X_eval, y_eval)
        rank_names = ("po_msp", "po_nll", "po_energy", "po_resid")
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
    fig, ax = plt.subplots(figsize=(max(8, 1.7 * len(packs)), 4.6))
    x = np.arange(len(packs))
    w = 0.15
    for i, m in enumerate(METHODS):
        vals = [all_res[p]["metrics"][m][key] for p in packs]
        ax.bar(x + (i - 2) * w, vals, w, label=m, color=SCORE_COLORS[m])
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
    lines = [
        "# Image-OOD: ViT / CLIP embedding + pseudo-outcome PO-risk",
        "",
        "Benchmark redesign: **freeze backbone → embedding → pseudo-outcome μ → PO-risk**.",
        "",
        "- `food101_vit`: frozen ViT embeddings + **class-holdout** (far-OOD).",
        "- CLIP packs: cached img embeddings + **domain shift** (near-OOD).",
        "",
        "### Scores",
        "",
        "| score | definition | needs y at test? |",
        "|---|---|---|",
        "| `po_msp` | `1 − max_c p_μ(c\\|x)` | no |",
        "| `po_energy` | Shannon entropy of μ (RF has no logits) | no |",
        "| `po_nll` | `1 − p_μ(y\\|x)` | yes |",
        "| `po_resid` | `|y − μ_reg(x)|` | yes |",
        "| `dre` | logistic domain score | no (mix) |",
        "",
        "## AUROC (ID vs OOD)",
        "",
        "| pack | shift | po_msp | po_nll | **po_resid** | po_energy | dre | best |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    wins = {m: 0 for m in METHODS}
    for pack, blob in all_res.items():
        m = blob["metrics"]
        best = max(METHODS, key=lambda k: m[k]["auroc"])
        wins[best] += 1
        shift = f"{blob['id_domain']}→{blob['ood_domain']}"
        lines.append(
            f"| `{pack}` | {shift} | {m['po_msp']['auroc']:.3f} | {m['po_nll']['auroc']:.3f} | "
            f"**{m['po_resid']['auroc']:.3f}** | {m['po_energy']['auroc']:.3f} | "
            f"{m['dre']['auroc']:.3f} | `{best}` |"
        )
    lines += [
        "",
        f"**AUROC wins:** " + ", ".join(f"`{k}`={v}" for k, v in wins.items()),
        "",
        "## FPR95 (↓ better)",
        "",
        "| pack | po_msp | po_nll | po_resid | po_energy | dre |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for pack, blob in all_res.items():
        m = blob["metrics"]
        lines.append(
            f"| `{pack}` | {m['po_msp']['fpr95']:.3f} | {m['po_nll']['fpr95']:.3f} | "
            f"{m['po_resid']['fpr95']:.3f} | {m['po_energy']['fpr95']:.3f} | {m['dre']['fpr95']:.3f} |"
        )
    lines += [
        "",
        "## Hard-rank (vs residual truth) — label-free scores",
        "",
        "| pack | spearman msp / energy / resid | P@20% msp / energy / resid |",
        "|---|---|---|",
    ]
    for pack, blob in all_res.items():
        h = blob["hard_rank"]
        lines.append(
            f"| `{pack}` | "
            f"{h['po_msp']['spearman']:.2f} / {h['po_energy']['spearman']:.2f} / {h['po_resid']['spearman']:.2f} | "
            f"{h['po_msp']['precision_at_k']:.2f} / {h['po_energy']['precision_at_k']:.2f} / {h['po_resid']['precision_at_k']:.2f} |"
        )
    lines += [
        "",
        "### Takeaway",
        "",
        "PO-risk on a **frozen ViT embedding** is enough for an image-OOD bench:",
        "no generative model, no fine-tune — just μ on embeddings, then residual / MSP.",
        "On **Food101 class-holdout** (far-OOD), prefer label-free `po_msp` / `po_energy`;",
        "DRE is a domain-separability upper bound on near-OOD packs, not a PO-risk substitute.",
        "`po_resid` / `po_nll` need labels (or label-free fallbacks on held-out classes).",
        "",
        "See `docs/agod/AGOD_image_ood_bench.md`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--packs",
        nargs="+",
        default=None,
        help="default: all available (food101_vit + CLIP packs)",
    )
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
            f"  AUROC msp={m['po_msp']['auroc']:.3f} resid={m['po_resid']['auroc']:.3f} "
            f"energy={m['po_energy']['auroc']:.3f} dre={m['dre']['auroc']:.3f} "
            f"meta={blob['meta']}",
            flush=True,
        )

    payload = {
        "pca_d": args.pca_d,
        "seed": args.seed,
        "extract_max_n": args.extract_max_n,
        "packs": list(all_res.keys()),
        "results": all_res,
        "note": "ViT/CLIP embedding + pseudo-outcome PO-risk image-OOD bench",
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
                "Image-OOD AUROC: PO-risk on frozen embeddings",
                "image_ood_auroc.png",
                ylim=(0.4, 1.02),
            )
        )
        plots.append(
            _bar(
                all_res,
                args.out,
                "fpr95",
                "FPR95 (↓)",
                "Image-OOD FPR95",
                "image_ood_fpr95.png",
            )
        )
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
