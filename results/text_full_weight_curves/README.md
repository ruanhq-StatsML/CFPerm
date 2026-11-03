# Feature vs sample weights — FULL holdout acc / AUC curves

Cold-start empowerment only: after a novel category reaches 200 samples, stop weighting it. Feature-level = X·(1+α·vimp); sample-level = membership w_i. FULL LIVE holdout. No EWMA.

KPI = **overall accuracy + mean novel OVA-AUC on full held-out LIVE** (not first-n). budget=novel/category [20, 40, 60, 80, 100, 150, 200, 250, 300, 400], cold_start_cap=**200** (stop weighting past cap), seeds=[2026, 2027, 2028], loss=`CE (log_loss)`, α_feat=3.0, **EWMA=off**.

## At cold-start cap (budget=200/cat) — overall acc

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.493 | 0.486 | 0.484 | 0.497 | -0.007 | -0.009 | 0.004 |
| news20_rec | 0.534 | 0.528 | 0.522 | 0.539 | -0.006 | -0.013 | 0.005 |
| news20_mix | 0.476 | 0.474 | 0.464 | 0.473 | -0.002 | -0.012 | -0.004 |
| cnae9_text | 0.892 | 0.919 | 0.882 | 0.913 | 0.027 | -0.010 | 0.021 |
| amazon_cat | 0.760 | 0.764 | 0.741 | 0.754 | 0.004 | -0.019 | -0.006 |

## Endpoint overall accuracy (full eval, last budget)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.502 | 0.499 | 0.504 | 0.510 | -0.002 | 0.002 | 0.008 |
| news20_rec | 0.545 | 0.550 | 0.550 | 0.570 | 0.006 | 0.006 | 0.026 |
| news20_mix | 0.504 | 0.507 | 0.502 | 0.510 | 0.003 | -0.003 | 0.006 |
| cnae9_text | 0.894 | 0.916 | 0.889 | 0.892 | 0.022 | -0.004 | -0.001 |
| amazon_cat | 0.737 | 0.751 | 0.747 | 0.744 | 0.014 | 0.010 | 0.007 |

mean overall acc: uniform=0.636  feature_w=0.645  sample_w=0.638  both=0.645
mean Δacc vs U: feature_w=0.009  sample_w=0.002  both=0.009

## Endpoint mean novel OVA-AUC (full eval)

| dataset | uniform | feature_w | sample_w | both | Δfeat | Δsample | Δboth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| news20_comp | 0.783 | 0.772 | 0.763 | 0.761 | -0.010 | -0.020 | -0.022 |
| news20_rec | 0.749 | 0.756 | 0.749 | 0.757 | 0.008 | 0.000 | 0.008 |
| news20_mix | 0.841 | 0.832 | 0.842 | 0.837 | -0.009 | 0.001 | -0.003 |
| cnae9_text | 0.968 | 0.976 | 0.973 | 0.965 | 0.008 | 0.005 | -0.003 |
| amazon_cat | 0.859 | 0.853 | 0.852 | 0.847 | -0.007 | -0.007 | -0.012 |

mean OVA-AUC: uniform=0.840  feature_w=0.838  sample_w=0.836  both=0.833
mean ΔOVA-AUC vs U: feature_w=-0.002  sample_w=-0.004  both=-0.006

## Mean curves (acc / OVA-AUC vs novel budget)

| budget | acc U | acc feat | acc sample | acc both | AUC U | AUC feat | AUC sample | AUC both |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20 | 0.456 | 0.461 | 0.468 | 0.483 | 0.729 | 0.722 | 0.727 | 0.729 |
| 40 | 0.534 | 0.522 | 0.531 | 0.541 | 0.769 | 0.758 | 0.766 | 0.759 |
| 60 | 0.567 | 0.584 | 0.572 | 0.574 | 0.781 | 0.784 | 0.781 | 0.784 |
| 80 | 0.568 | 0.574 | 0.567 | 0.563 | 0.779 | 0.783 | 0.780 | 0.771 |
| 100 | 0.595 | 0.597 | 0.588 | 0.587 | 0.806 | 0.802 | 0.796 | 0.791 |
| 150 | 0.624 | 0.619 | 0.612 | 0.617 | 0.820 | 0.818 | 0.814 | 0.812 |
| 200 | 0.631 | 0.634 | 0.619 | 0.635 | 0.829 | 0.826 | 0.821 | 0.824 |
| 250 | 0.631 | 0.641 | 0.636 | 0.640 | 0.833 | 0.836 | 0.834 | 0.832 |
| 300 | 0.643 | 0.647 | 0.632 | 0.645 | 0.837 | 0.836 | 0.830 | 0.829 |
| 400 | 0.636 | 0.645 | 0.638 | 0.645 | 0.840 | 0.838 | 0.836 | 0.833 |

```bash
python3 scripts/run_text_full_weight_curves.py --multiseed
```
