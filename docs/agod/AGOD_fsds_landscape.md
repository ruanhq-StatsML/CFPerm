# FSDS + instance-disc + proto-drift + grad-memory (MSR-VTT)

## Logic

1. **FSDS hybrid**: PO (concept) - MMD (covariate) as the base attribution score.
2. **Nonparametric instance discrimination**: kNN domain purity on each modality block
   -> feature-level separability sensor (no deep NCE).
3. **Prototype drift**: EMA modality centroids; d=1-cos(proto, batch_mean);
   high PO + high proto-drift => true concept move; high MMD + low proto-drift => mush.
4. **Grad-projection memory + landscape**: bank of common-dim grad signatures;
   residual energy + Gram erank -> LR gains (new direction up, collinear basin down).
5. **Proportion + next-step**: alpha_video/alpha_text/alpha_audio now; a-hat from EMA+sensor delta;
   fused policy uses a-hat * landscape_gain as the adapt step.

## Smoke (MSR-VTT packed, 3-mod)

| policy | Acc up | Acc post | prop MAE | erank | alpha (v/t/a) |
|---|---:|---:|---:|---:|---|
| `equal` | +0.032 | 0.524 | 0.000 | 3.00 | video:0.33 / text:0.33 / audio:0.33 |
| `fsds` | +0.007 | 0.509 | 0.051 | 3.00 | video:0.54 / text:0.35 / audio:0.11 |
| `fsds_disc_proto` | +0.020 | 0.524 | 0.057 | 3.00 | video:0.63 / text:0.27 / audio:0.10 |
| `fused_landscape` | +0.012 | 0.503 | 0.057 | 3.00 | video:0.63 / text:0.27 / audio:0.10 |

Best Acc up: **`equal`** (+0.032).

```bash
PYTHONPATH=. python3 scripts/run_agod_fsds_landscape_smoke.py
```
