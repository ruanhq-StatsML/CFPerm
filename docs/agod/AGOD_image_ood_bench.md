# Image-OOD: PO-risk μ vs classical embedding scores (no DRE)

Freeze backbone → embedding → **pseudo-outcome μ** (PO-risk) vs Mahalanobis / kNN / centroid.

- `food101_vit`: frozen ViT + **class-holdout** (far-OOD) — primary.
- CLIP packs: cached img embeddings + domain shift (near-OOD).

### Scores

| score | definition |
|---|---|
| `po_msp` | `1 − max_c p_μ(c\|x)` — PO-risk (label-free) |
| `po_energy` | Shannon entropy of μ — PO-risk (label-free) |
| `po_nll` | `1 − p_μ(y\|x)` (MSP fallback if y unseen) |
| `maha` | min class-conditional Mahalanobis (Lee et al.) |
| `knn` | mean L2 to 5-NN in ID-train (Sun et al.) |
| `centroid` | L2 to nearest ID class centroid |

## AUROC (ID vs OOD)

| pack | shift | po_msp | po_energy | po_nll | maha | knn | centroid | best |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `food101_vit` | id_classes→heldout_classes | 0.630 | 0.610 | 0.508 | 0.697 | 0.639 | 0.706 | `centroid` |
| `coco_outdoor_indoor` | outdoor→indoor | 0.430 | 0.442 | 0.440 | 0.387 | 0.669 | 0.051 | `knn` |
| `coco_time_order` | early_id→late_id | 0.494 | 0.498 | 0.505 | 0.501 | 0.534 | 0.500 | `knn` |
| `indiana_cxr` | Frontal→Lateral | 0.708 | 0.844 | 0.600 | 0.994 | 0.997 | 0.985 | `knn` |
| `fashion_iq` | train→test | 0.518 | 0.578 | 0.530 | 0.499 | 0.645 | 0.485 | `knn` |

**AUROC wins:** `po_msp`=0, `po_energy`=0, `po_nll`=0, `maha`=0, `knn`=4, `centroid`=1

## FPR95 (↓ better)

| pack | po_msp | po_energy | po_nll | maha | knn | centroid |
|---|---:|---:|---:|---:|---:|---:|
| `food101_vit` | 0.748 | 0.784 | 0.756 | 0.843 | 0.909 | 0.820 |
| `coco_outdoor_indoor` | 0.933 | 0.936 | 0.943 | 0.963 | 0.585 | 1.000 |
| `coco_time_order` | 0.944 | 0.947 | 0.940 | 0.952 | 0.919 | 0.945 |
| `indiana_cxr` | 0.708 | 0.467 | 0.773 | 0.016 | 0.015 | 0.065 |
| `fashion_iq` | 0.921 | 0.919 | 0.909 | 0.963 | 0.929 | 0.935 |

### Takeaway

Compare **PO-risk on a frozen embedding** to standard embedding OOD detectors.
No DRE — domain classifiers are not part of this comparison.

See `docs/agod/AGOD_image_ood_bench.md`.
