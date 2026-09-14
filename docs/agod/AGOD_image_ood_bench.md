# Image-OOD: ViT / CLIP embedding + pseudo-outcome PO-risk

Benchmark redesign: **freeze backbone → embedding → pseudo-outcome μ → PO-risk**.

- `food101_vit`: frozen ViT embeddings + **class-holdout** (far-OOD).
- CLIP packs: cached img embeddings + **domain shift** (near-OOD).

### Scores

| score | definition | needs y at test? |
|---|---|---|
| `po_msp` | `1 − max_c p_μ(c\|x)` | no |
| `po_energy` | energy of class probs | no |
| `po_nll` | `1 − p_μ(y\|x)` | yes |
| `po_resid` | `|y − μ_reg(x)|` | yes |
| `dre` | logistic domain score | no (mix) |

## AUROC (ID vs OOD)

| pack | shift | po_msp | po_nll | **po_resid** | po_energy | dre | best |
|---|---|---:|---:|---:|---:|---:|---|
| `coco_outdoor_indoor` | outdoor→indoor | 0.430 | 0.440 | **0.445** | 0.459 | 0.991 | `dre` |
| `coco_time_order` | early_id→late_id | 0.494 | 0.505 | **0.507** | 0.490 | 0.579 | `dre` |
| `indiana_cxr` | Frontal→Lateral | 0.708 | 0.600 | **0.649** | 0.468 | 1.000 | `dre` |
| `fashion_iq` | train→test | 0.518 | 0.530 | **0.522** | 0.507 | 0.567 | `dre` |

**AUROC wins:** `po_msp`=0, `po_nll`=0, `po_resid`=0, `po_energy`=0, `dre`=4

## FPR95 (↓ better)

| pack | po_msp | po_nll | po_resid | po_energy | dre |
|---|---:|---:|---:|---:|---:|
| `coco_outdoor_indoor` | 0.933 | 0.943 | 0.943 | 0.992 | 0.021 |
| `coco_time_order` | 0.944 | 0.940 | 0.943 | 0.997 | 0.917 |
| `indiana_cxr` | 0.708 | 0.773 | 0.908 | 0.939 | 0.000 |
| `fashion_iq` | 0.921 | 0.909 | 0.936 | 0.983 | 0.938 |

## Hard-rank (vs residual truth) — label-free scores

| pack | spearman msp / energy / resid | P@20% msp / energy / resid |
|---|---|---|
| `coco_outdoor_indoor` | 0.47 / -0.10 / 1.00 | 0.38 / 0.17 / 1.00 |
| `coco_time_order` | 0.22 / -0.05 / 1.00 | 0.29 / 0.17 / 1.00 |
| `indiana_cxr` | 0.13 / -0.03 / 1.00 | 0.21 / 0.19 / 1.00 |
| `fashion_iq` | 0.22 / -0.06 / 1.00 | 0.19 / 0.20 / 1.00 |

### Takeaway

PO-risk on a **frozen ViT embedding** is enough for an image-OOD bench:
no generative model, no fine-tune — just μ on embeddings, then residual / MSP.
Prefer label-free `po_msp` / `po_energy` for far-OOD (class-holdout);
`po_resid` / `po_nll` for near-OOD domain shift when labels exist.

See `docs/agod/AGOD_image_ood_bench.md`.
