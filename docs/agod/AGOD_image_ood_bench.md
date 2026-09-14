# Image-OOD: ViT / CLIP embedding + pseudo-outcome PO-risk

Benchmark redesign: **freeze backbone → embedding → pseudo-outcome μ → PO-risk**.

- `food101_vit`: frozen ViT embeddings + **class-holdout** (far-OOD) — primary.
- CLIP packs: cached img embeddings + **domain shift** (near-OOD).

### Scores

| score | definition | needs y at test? |
|---|---|---|
| `po_msp` | `1 − max_c p_μ(c\|x)` | no |
| `po_energy` | Shannon entropy of μ (RF has no logits) | no |
| `po_nll` | `1 − p_μ(y\|x)` (MSP fallback if y unseen) | yes* |
| `po_resid` | `|y − μ_reg(x)|` (n/a on class-holdout OOD) | yes |
| `dre` | logistic domain score | no (mix) |

## AUROC (ID vs OOD)

| pack | shift | po_msp | po_nll | po_resid | po_energy | dre | best PO |
|---|---|---:|---:|---:|---:|---:|---|
| `food101_vit` | id_classes→heldout_classes | 0.630 | 0.508 | — | 0.610 | 0.706 | `po_msp` |
| `coco_outdoor_indoor` | outdoor→indoor | 0.430 | 0.440 | 0.445 | 0.442 | 0.991 | `po_resid` |
| `coco_time_order` | early_id→late_id | 0.494 | 0.505 | 0.507 | 0.498 | 0.579 | `po_resid` |
| `indiana_cxr` | Frontal→Lateral | 0.708 | 0.600 | 0.649 | 0.844 | 1.000 | `po_energy` |
| `fashion_iq` | train→test | 0.518 | 0.530 | 0.522 | 0.578 | 0.567 | `po_energy` |

**PO AUROC wins (excl. dre):** `po_msp`=1, `po_nll`=0, `po_resid`=2, `po_energy`=2

## FPR95 (↓ better)

| pack | po_msp | po_nll | po_resid | po_energy | dre |
|---|---:|---:|---:|---:|---:|
| `food101_vit` | 0.748 | 0.756 | — | 0.784 | 0.862 |
| `coco_outdoor_indoor` | 0.933 | 0.943 | 0.943 | 0.936 | 0.021 |
| `coco_time_order` | 0.944 | 0.940 | 0.943 | 0.947 | 0.917 |
| `indiana_cxr` | 0.708 | 0.773 | 0.908 | 0.467 | 0.000 |
| `fashion_iq` | 0.921 | 0.909 | 0.936 | 0.919 | 0.938 |

## Hard-rank

Class-holdout truth = binary OOD; domain-shift truth = residual hardness.

| pack | spearman msp / energy / resid | P@20% msp / energy / resid |
|---|---|---|
| `food101_vit` | 0.169 / 0.129 / — | 0.162 / 0.159 / — |
| `coco_outdoor_indoor` | 0.469 / 0.455 / 1.000 | 0.380 / 0.367 / 1.000 |
| `coco_time_order` | 0.217 / 0.200 / 1.000 | 0.288 / 0.272 / 1.000 |
| `indiana_cxr` | 0.128 / 0.179 / 1.000 | 0.211 / 0.226 / 1.000 |
| `fashion_iq` | 0.220 / 0.171 / 1.000 | 0.188 / 0.198 / 1.000 |

### Takeaway

PO-risk on a **frozen ViT embedding** is enough for an image-OOD bench:
no generative model, no fine-tune — just μ on embeddings, then MSP / entropy.
On **Food101 class-holdout** (far-OOD), label-free `po_msp` / `po_energy` are the PO scores;
`dre` is a domain-separability upper bound on near-OOD packs, not a PO-risk substitute.

See `docs/agod/AGOD_image_ood_bench.md`.
