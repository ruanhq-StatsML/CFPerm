# Feature-level vs sample-level weights — FULL accuracy

Compare feature-level VIMP scales vs sample-level CE weights on FULL held-out LIVE accuracy for text multiclass category shift.

KPI = **overall / known / novel accuracy on full held-out LIVE** (not first-n subset). n_per_novel≈**400**, seeds=[2026, 2027, 2028], loss=`CE (log_loss)`, α_feat=3.0.

## Overall accuracy (full eval)

| dataset | uniform | sample_w | feature_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.449 | 0.454 | 0.460 | 0.452 | 0.011 | 0.006 | 0.004 |
| news20_rec | 0.553 | 0.554 | 0.562 | 0.557 | 0.008 | 0.001 | 0.004 |
| news20_mix | 0.488 | 0.498 | 0.493 | 0.496 | 0.005 | 0.010 | 0.008 |
| cnae9_text | 0.884 | 0.888 | 0.905 | 0.912 | 0.021 | 0.004 | 0.028 |
| amazon_cat | 0.751 | 0.747 | 0.763 | 0.755 | 0.012 | -0.004 | 0.004 |

mean overall acc: U=0.625  sample=0.628  feature=0.636  both=0.635
mean Δacc vs U: feature=0.012  sample=0.003  both=0.010

## Novel-class accuracy (full eval)

| dataset | uniform | sample_w | feature_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.479 | 0.485 | 0.490 | 0.483 | 0.011 | 0.006 | 0.004 |
| news20_rec | 0.587 | 0.596 | 0.597 | 0.597 | 0.010 | 0.009 | 0.010 |
| news20_mix | 0.532 | 0.546 | 0.535 | 0.540 | 0.003 | 0.014 | 0.008 |
| cnae9_text | 0.883 | 0.896 | 0.904 | 0.922 | 0.021 | 0.013 | 0.039 |
| amazon_cat | 0.772 | 0.785 | 0.782 | 0.785 | 0.010 | 0.013 | 0.013 |

mean novel_acc: U=0.651  sample=0.662  feature=0.662  both=0.666
mean Δnovel vs U: feature=0.011  sample=0.011  both=0.015

## Known retention (full eval)

| dataset | uniform | sample_w | feature_w | both |
| --- | --- | --- | --- | --- |
| news20_comp | 0.188 | 0.193 | 0.197 | 0.192 |
| news20_rec | 0.282 | 0.213 | 0.277 | 0.240 |
| news20_mix | 0.080 | 0.054 | 0.109 | 0.090 |
| cnae9_text | 0.889 | 0.870 | 0.906 | 0.889 |
| amazon_cat | 0.637 | 0.541 | 0.661 | 0.593 |

```bash
python3 scripts/run_text_feature_sample_weights.py --multiseed
```
