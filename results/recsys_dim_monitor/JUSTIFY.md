# Justify: recsys stream × grain — OnlineRFPerm / RFPerm / FSDS

WHEN = frozen RF, T=MSE−E_ref, last-two hop; rank-p into ADDIS (primary) / SAFFRON (contrast).
Delay = first rejection − onset. FAR = mark before labeled onset. miss = never marked.
WHICH columns = FSDS + RFPerm ΔMSE + CFPerm φ-VIMP. Ranking metric = Kendall-τ vs planted magnitude.
WHICH accounts = south vs north own-ref, each with MMD / CMean / PO-risk. Y never a feature. Localization, not unique decomp. Not a graph method.

onset_true=2. n_ref=400, n_new=120, batches=6.

| kind | dim | ADDIS | FSDS recovered | south MMD | south ‖ΔX‖ | south PO | south ΔE[Y] | north MMD | north PO |
|---|---|---|---|---|---|---|---|---|---|
| covariate_south | order | 2 (d=0) | amount,channel | 0.35 | 2.55 | 7.59e-05 | 0.349 | 0.00109 | 0.000898 |
| covariate_south | merchant | 5 (d=3) | merchant_gmv | 0.291 | 1.7 | 0.000475 | 0.349 | -0.00538 | 0.000359 |
| covariate_south | user | 2 (d=0) | — | -0.00133 | 0.0659 | 0.00329 | 0.349 | -0.000738 | 0.000486 |
| covariate_south | all | FAR@0 | amount,channel,merchant_gmv | 0.296 | 3.07 | 0.000224 | 0.349 | 0.000566 | 0.00116 |
| concept_south | order | 2 (d=0) | amount | -0.00137 | 0.275 | 0.00307 | -0.0412 | 0.00109 | 0.000898 |
| concept_south | merchant | miss | — | -0.00969 | 0.059 | 0.000936 | -0.0412 | -0.00538 | 0.000359 |
| concept_south | user | 2 (d=0) | — | -0.00133 | 0.0659 | 0.000132 | -0.0412 | -0.000738 | 0.000486 |
| concept_south | all | FAR@0 | amount | -0.00413 | 0.289 | 0.00376 | -0.0412 | 0.000566 | 0.00116 |
| both | order | 2 (d=0) | amount,channel | 0.35 | 2.55 | 0.0015 | -0.262 | 0.00109 | 0.000898 |
| both | merchant | miss | merchant_gmv | 0.291 | 1.7 | 0.000157 | -0.262 | -0.00538 | 0.000359 |
| both | user | 2 (d=0) | — | -0.00133 | 0.0659 | 0.00183 | -0.262 | -0.000738 | 0.000486 |
| both | all | FAR@0 | amount,channel,merchant_gmv | 0.296 | 3.07 | 0.000552 | -0.262 | 0.000566 | 0.00116 |
