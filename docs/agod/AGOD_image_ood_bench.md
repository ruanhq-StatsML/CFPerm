# Verdict: PO-risk is **not** suited to image-OOD

**Conclusion:** freeze ViT → embedding → pseudo-outcome μ → PO-risk does **not**
make a good image-OOD detector. Observation-level PO-risk is too noisy; if the
goal is binary ID vs OOD, an RF binary classifier (or kNN / Mahalanobis on the
same embedding) is the right tool. Stop here for image-OOD.

## Why

1. **Obs-level PO AUROC is weak / noisy** — Food101 class-holdout: `po_msp=0.630`
   vs `knn=0.639` / `centroid=0.706` / oracle `rf_binary=0.753`.
2. **Near-OOD domain packs** — PO stays near chance (`~0.43–0.58`); kNN or
   oracle RF-binary dominate when domains separate in embedding space.
3. **PO-risk’s job in AGOD** is batch reject → hard-sample reweight on streams,
   not sample-level image OOD detection.

## Evidence (obs-level AUROC, no DRE)

| pack | po_msp | po_energy | maha | knn | centroid | rf_binary† |
|---|---:|---:|---:|---:|---:|---:|
| `food101_vit` (far-OOD) | 0.630 | 0.610 | 0.697 | 0.639 | 0.706 | **0.753** |
| `coco_outdoor_indoor` | 0.430 | 0.442 | 0.387 | 0.669 | 0.051 | **0.996** |
| `coco_time_order` | 0.494 | 0.498 | 0.501 | 0.534 | 0.500 | **0.606** |
| `indiana_cxr` | 0.708 | 0.844 | 0.994 | 0.997 | 0.985 | **1.000** |
| `fashion_iq` | 0.518 | 0.578 | 0.499 | 0.645 | 0.485 | **0.739** |

† `rf_binary` = oracle ceiling (trained with OOD labels). Fair deployable
baselines already beat PO on the primary far-OOD pack.

## What not to do

- Do not sell PO-MSP / PO-energy as an image-OOD method.
- Do not use obs-level image-OOD AUROC as a success metric for AGOD PO-risk.

## What PO-risk is for

Stream packs after OnlineRFPerm reject: recent-control μ₀ → PO on the OOD
batch → soft reweight / hard-rank. That pipeline stays; image-OOD does not.

Code/bench kept under `agod/image_ood.py` + `scripts/run_agod_image_ood_bench.py`
only as a **negative result** / ablation.
