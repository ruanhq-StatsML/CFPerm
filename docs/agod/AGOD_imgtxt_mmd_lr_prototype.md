# Image/Text AGOD portable prototype

## Gate semantics

- **Yes, dropout-like**: hard-gate zeros modality **proj gradients** when `α_m < θ` → structured **adapt-dropout** / sparse update.
- **Not** "OOD large ⇒ skip inference": **FWD still runs**; only adapt BWD for low-α modalities is skipped. High covariate often lowers α (Acc rule: cov↑ → LR↓), so gate fires more often under strong OOD.

## Smoke board

| Dataset | mods | B1 | B5 | B5g | B5r | B5g flops | B5r flops | B5g−B5r |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| coco_outdoor_indoor | image,text,bbox | +0.034 | +0.019 | +0.031 | +0.018 | 0.810 | 0.841 | +0.013 |
| coco_time_order | image,text,bbox | +0.012 | +0.028 | +0.007 | +0.001 | 0.778 | 0.841 | +0.006 |
| coco_center_split | image,text,bbox | +0.012 | +0.028 | +0.007 | +0.001 | 0.778 | 0.841 | +0.006 |
| mm_train_test | image,text | -0.004 | -0.022 | +0.010 | -0.007 | 0.778 | 0.911 | +0.018 |
| fashion_iq | image,text | -0.004 | -0.022 | +0.010 | -0.007 | 0.778 | 0.911 | +0.018 |
| indiana_cxr | image,text | +0.024 | +0.016 | +0.007 | +0.025 | 0.778 | 0.911 | -0.018 |
| microscopy_clip | image,text | +0.040 | +0.027 | +0.018 | +0.039 | 0.822 | 0.911 | -0.021 |

```bash
PYTHONPATH=. python3 scripts/run_agod_imgtxt_mmd_lr.py
```
