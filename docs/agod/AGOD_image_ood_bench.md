# Image-OOD: PO-risk is batch / hard-rank — not obs-level detection

**Point:** observation-level PO-risk AUROC is too noisy for detection.
If you want a binary OOD detector, use an **RF binary classifier** (oracle ceiling).
PO-risk belongs to AGOD as a **batch / hard-sample** score after OnlineRFPerm.

- `food101_vit`: frozen ViT + class-holdout (far-OOD) — primary.
- CLIP packs: cached embeddings + domain shift.
- `rf_binary`: oracle ID vs OOD RF (**uses OOD labels at train**) — ceiling, not a fair PO rival.

## Primary: batch-mean AUROC (batch_size=32)

| pack | po_msp | po_energy | maha | knn | centroid | rf_binary† | best deploy |
|---|---:|---:|---:|---:|---:|---:|---|
| `food101_vit` | 1.000 | 0.996 | 0.997 | 0.972 | 0.997 | 1.000 | `po_msp` |
| `coco_outdoor_indoor` | 0.167 | 0.219 | 0.036 | 1.000 | 0.000 | 1.000 | `knn` |
| `coco_time_order` | 0.467 | 0.512 | 0.497 | 0.710 | 0.464 | 0.962 | `knn` |
| `indiana_cxr` | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | `po_msp` |
| `fashion_iq` | 0.643 | 0.887 | 0.564 | 0.963 | 0.531 | 1.000 | `knn` |

**Deployable batch-AUROC wins:** `po_msp`=2, `po_energy`=0, `maha`=0, `knn`=3, `centroid`=0

† `rf_binary` = oracle ceiling (OOD labels at train).

## Diagnostic: obs-level AUROC (noisy — do not prefer)

| pack | po_msp | po_energy | maha | knn | centroid | rf_binary† |
|---|---:|---:|---:|---:|---:|---:|
| `food101_vit` | 0.630 | 0.610 | 0.697 | 0.639 | 0.706 | 0.753 |
| `coco_outdoor_indoor` | 0.430 | 0.442 | 0.387 | 0.669 | 0.051 | 0.996 |
| `coco_time_order` | 0.494 | 0.498 | 0.501 | 0.534 | 0.500 | 0.606 |
| `indiana_cxr` | 0.708 | 0.844 | 0.994 | 0.997 | 0.985 | 1.000 |
| `fashion_iq` | 0.518 | 0.578 | 0.499 | 0.645 | 0.485 | 0.739 |

## Hard-rank (PO job): spearman / P@20%

| pack | po_msp | po_energy | maha | knn | centroid | rf_binary† |
|---|---|---|---|---|---|---|
| `food101_vit` | 0.169/0.162 | 0.129/0.159 | 0.333/0.347 | 0.236/0.340 | 0.325/0.326 | 0.369/0.359 |
| `coco_outdoor_indoor` | 0.469/0.380 | 0.455/0.367 | 0.211/0.259 | 0.173/0.251 | 0.134/0.228 | -0.109/0.140 |
| `coco_time_order` | 0.217/0.288 | 0.200/0.272 | 0.129/0.212 | 0.064/0.207 | 0.027/0.184 | -0.011/0.156 |
| `indiana_cxr` | 0.128/0.211 | 0.179/0.226 | 0.256/0.270 | 0.257/0.276 | 0.220/0.257 | 0.224/0.237 |
| `fashion_iq` | 0.220/0.188 | 0.171/0.198 | 0.039/0.201 | -0.058/0.217 | 0.075/0.230 | -0.014/0.176 |

### Takeaway

1. Obs-level PO-risk AUROC is a bad primary — too noisy; RF-binary dominates if OOD labels exist.
2. Batch-mean AUROC is closer to the OnlineRFPerm reject unit.
3. Hard-rank is the AGOD-native PO question: does μ-risk surface the hard rows?

See `docs/agod/AGOD_image_ood_bench.md`.
