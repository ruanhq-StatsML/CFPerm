# Feature / sample / EWMA weights — FULL accuracy

EWMA decay w=exp(-λ·n_idx) puts CE mass on earliest novel-category hits; compare FULL overall/novel accuracy and incremental Δacc(k) vs uniform.

KPI = **overall / known / novel accuracy on full held-out LIVE** (not first-n). n_per_novel≈**400**, seeds=[2026, 2027, 2028], loss=`CE (log_loss)`, α_feat=3.0, EWMA `w=exp(-λ·n_idx)` with λ=**0.02**.

## Overall accuracy (full eval)

| dataset | U | sample | feature | both | **ewma** | ewma+feat | Δewma | Δfeat |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.479 | 0.485 | 0.484 | 0.488 | 0.372 | 0.379 | -0.107 | 0.005 |
| news20_rec | 0.546 | 0.538 | 0.540 | 0.542 | 0.431 | 0.474 | -0.115 | -0.006 |
| news20_mix | 0.492 | 0.493 | 0.487 | 0.491 | 0.397 | 0.381 | -0.094 | -0.005 |
| cnae9_text | 0.891 | 0.894 | 0.899 | 0.894 | 0.884 | 0.889 | -0.007 | 0.008 |
| amazon_cat | 0.744 | 0.759 | 0.761 | 0.771 | 0.679 | 0.713 | -0.065 | 0.016 |

mean overall acc: uniform=0.630  sample_w=0.634  feature_w=0.634  both=0.637  ewma=0.553  ewma_feat=0.567
mean Δacc vs U: sample_w=0.003  feature_w=0.004  both=0.007  ewma=-0.078  ewma_feat=-0.063

## Novel-class accuracy (full eval)

| dataset | U | sample | feature | both | **ewma** | ewma+feat | Δewma | Δfeat |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.510 | 0.521 | 0.509 | 0.512 | 0.379 | 0.380 | -0.131 | -0.001 |
| news20_rec | 0.584 | 0.584 | 0.570 | 0.579 | 0.445 | 0.482 | -0.138 | -0.013 |
| news20_mix | 0.538 | 0.543 | 0.528 | 0.537 | 0.419 | 0.398 | -0.118 | -0.009 |
| cnae9_text | 0.895 | 0.908 | 0.907 | 0.901 | 0.871 | 0.886 | -0.024 | 0.012 |
| amazon_cat | 0.798 | 0.816 | 0.813 | 0.828 | 0.694 | 0.731 | -0.105 | 0.015 |

mean novel_acc: uniform=0.665  sample_w=0.674  feature_w=0.666  both=0.671  ewma=0.562  ewma_feat=0.575
mean Δnovel vs U: sample_w=0.009  feature_w=0.001  both=0.006  ewma=-0.103  ewma_feat=-0.090

## Incremental Δacc(k) — EWMA early-stage lift

Train on known + first-k novel hits/category; score fixed FULL held-out eval.

| k | Δewma | Δewma+feat | Δfeature | acc U | acc ewma |
| --- | --- | --- | --- | --- | --- |
| 10 | -0.005 | -0.001 | 0.008 | 0.422 | 0.418 |
| 25 | -0.009 | -0.011 | -0.000 | 0.501 | 0.491 |
| 50 | -0.006 | -0.002 | 0.009 | 0.549 | 0.542 |
| 100 | -0.018 | -0.010 | 0.013 | 0.499 | 0.481 |
| 150 | -0.023 | -0.010 | 0.023 | 0.514 | 0.492 |
| 200 | -0.040 | -0.037 | 0.007 | 0.536 | 0.496 |

```bash
python3 scripts/run_text_feature_sample_weights.py --multiseed
```
