# Image/Text AGOD portable prototype

## Gate semantics

- **Yes, dropout-like**: hard-gate zeros modality **proj gradients** when `α_m < θ` → structured **adapt-dropout** / sparse update.
- **Not** "OOD large ⇒ skip inference": **FWD still runs**; only adapt BWD for low-α modalities is skipped. High covariate often lowers α (Acc rule: cov↑ → LR↓), so gate fires more often under strong OOD.

## Smoke board

| Dataset | mods | B1 lift | B5 lift | B5g lift | B5g flops | B5−B1 |
|---|---|---:|---:|---:|---:|---:|
| coco_outdoor_indoor | image,text,bbox | +0.034 | +0.019 | +0.031 | 0.810 | -0.015 |
| mm_train_test | image,text | -0.004 | -0.022 | +0.010 | 0.778 | -0.018 |
| indiana_cxr | image,text | +0.024 | +0.016 | +0.007 | 0.778 | -0.007 |

```bash
PYTHONPATH=. python3 scripts/run_agod_imgtxt_mmd_lr.py
```

LaTeX tables (Amazon / MSR-VTT / ImgTxt board):
`docs/agod/AGOD_portable_prototype_tables_only.tex`
