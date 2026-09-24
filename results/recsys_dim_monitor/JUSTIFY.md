# Justify: recsys stream × grain — OnlineRFPerm / RFPerm / FSDS

WHEN = frozen RF, T=MSE−E_ref, last-two hop; rank-p into ADDIS (primary) / SAFFRON (contrast).
Delay = first rejection − onset. FAR = mark before labeled onset. miss = never marked.
WHICH columns = FSDS + RFPerm ΔMSE + CFPerm φ-VIMP. Ranking metric = Kendall-τ vs planted magnitude.
WHICH accounts = south vs north own-ref, each with MMD / CMean / PO-risk. Y never a feature. Localization, not unique decomp. Not a graph method.

onset_true=2. n_ref=400, n_new=120, batches=6.

| kind | dim | ADDIS | FSDS recovered | south MMD | south ‖ΔX‖ | south PO | south ΔE[Y] | north MMD | north PO |
|---|---|---|---|---|---|---|---|---|---|
| covariate_south | order | 2 (d=0) | amount,channel | 0.35 | 2.55 | 7.59e-05 | 0.349 | 0.00109 | 0.000898 |
| covariate_south | merchant | 2 (d=0) | merchant_gmv | 0.465 | 1.45 | 0.00107 | 0.349 | 0.00333 | 0.000141 |
| covariate_south | user | 2 (d=0) | — | -0.00133 | 0.0659 | 0.00329 | 0.349 | -0.000738 | 0.000486 |
| covariate_south | all | 2 (d=0) | merchant_gmv,amount,channel | 0.323 | 2.94 | 0.000263 | 0.349 | 0.00155 | 0.000828 |
| concept_south | order | 2 (d=0) | amount | -0.00137 | 0.275 | 0.00307 | -0.0412 | 0.00109 | 0.000898 |
| concept_south | merchant | 2 (d=0) | — | 0.00546 | 0.063 | 0.00029 | -0.0412 | 0.00333 | 0.000141 |
| concept_south | user | 2 (d=0) | — | -0.00133 | 0.0659 | 0.000132 | -0.0412 | -0.000738 | 0.000486 |
| concept_south | all | 2 (d=0) | amount | -0.00302 | 0.29 | 0.00295 | -0.0412 | 0.00155 | 0.000828 |
| both | order | 2 (d=0) | amount,channel | 0.35 | 2.55 | 0.0015 | -0.262 | 0.00109 | 0.000898 |
| both | merchant | 2 (d=0) | merchant_gmv | 0.465 | 1.45 | 7.9e-06 | -0.262 | 0.00333 | 0.000141 |
| both | user | 2 (d=0) | — | -0.00133 | 0.0659 | 0.00183 | -0.262 | -0.000738 | 0.000486 |
| both | all | 2 (d=0) | merchant_gmv,amount,channel | 0.323 | 2.94 | 0.000139 | -0.262 | 0.00155 | 0.000828 |
