# Smooth drift-vs-noise router (Amazon)

## Rule incorporated

1. **Drift vs noise gate** (per modality):
   - boost only if `PO` high **and** unimodal Acc not bad **and** VIMP/Fisher not low
   - else damp PO (noise path) / near adapter-only LR floor
2. **Smooth control law**:
   ```
   g = normalize(w1·PO_gated + w2·MMD + w3·VIMP)
   # combine: score = normalize(PO_gated) − normalize(MMD·(1+VIMP))
   α = EMA(Softmax(g / τ))
   w = w0 · (β + (1-β)·α·|M|)
   ```

## Smoke

| policy | MSE pre | MSE post | MSE drop↑ | Acc↑ | Acc post | frac signal |
|---|---:|---:|---:|---:|---:|---:|
| `equal` | 0.2580 | 0.2496 | +0.0083 | +0.069 | 0.494 | nan |
| `msg_softmax` | 0.2588 | 0.2484 | +0.0104 | +0.059 | 0.494 | nan |
| `concept_minus_cov` | 0.2587 | 0.2481 | +0.0106 | +0.062 | 0.494 | nan |
| `smooth_drift_noise` | 0.2581 | 0.2487 | +0.0094 | +0.066 | 0.490 | 0.17 |
| `smooth_concept_cov` | 0.2569 | 0.2476 | +0.0092 | +0.066 | 0.490 | 0.33 |

Best MSE-drop: **`concept_minus_cov`** (+0.0106).
Best Acc↑: **`equal`** (+0.069).

Calibration note: Fisher/VIMP uses **relative** cross-mod quantile first;
flat VIMP (common on small Amazon RF probes) does not auto-damp high-PO
modalities. Absolute `vimp_floor` is soft. Noise-gated LR keeps
`adapter_lr_keep` of free mass above the β floor (not a hard kill).

```bash
PYTHONPATH=. python3 scripts/run_agod_smooth_drift_noise.py
```
