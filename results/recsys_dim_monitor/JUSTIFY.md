# Justify: recsys stream × grain — OnlineRFPerm / RFPerm / FSDS

WHEN = frozen RF, T=MSE−E_ref, last-two hop; rank-p into ADDIS (primary) / SAFFRON (contrast).
Delay = first rejection − onset. FAR = mark before labeled onset. miss = never marked.
WHICH columns = FSDS + RFPerm ΔMSE + CFPerm φ-VIMP. Ranking metric = Kendall-τ vs planted magnitude.
WHICH accounts = south vs north own-ref. Y never a feature. Localization, not unique decomp. Not a graph method.

onset_true=2. n_ref=400, n_new=120, batches=6.

| kind | dim | ADDIS | rank-p | hop | FSDS recovered | τ_FSDS | τ_RFPerm | τ_CFPerm | south MMD | north MMD |
|---|---|---|---|---|---|---|---|---|---|---|
| covariate_south | order | 2 (d=0) | 2 (d=0) | miss | amount,channel | 0.913 | 0.548 | -0.183 | 0.35 | 0.00109 |
| covariate_south | merchant | 5 (d=3) | 2 (d=0) | miss | merchant_gmv | 0.816 | -0.816 | -0.816 | 0.291 | -0.00538 |
| covariate_south | user | 2 (d=0) | FAR@0 | miss | — |  |  |  | -0.00133 | -0.000738 |
| covariate_south | all | FAR@0 | FAR@0 | miss | amount,channel,merchant_gmv | 0.691 | 0.182 | -0.182 | 0.296 | 0.000566 |
| concept_south | order | 2 (d=0) | 2 (d=0) | miss | amount | 0.236 | 0.236 | 0.707 | -0.00137 | 0.00109 |
| concept_south | merchant | miss | 2 (d=0) | miss | — |  |  |  | -0.00969 | -0.00538 |
| concept_south | user | 2 (d=0) | FAR@0 | miss | — |  |  |  | -0.00133 | -0.000738 |
| concept_south | all | FAR@0 | FAR@0 | miss | amount | 0.354 | 0.354 | 0.471 | -0.00413 | 0.000566 |
| both | order | 2 (d=0) | 2 (d=0) | miss | amount,channel | 0.913 | -0.183 | 0.913 | 0.35 | 0.00109 |
| both | merchant | miss | 2 (d=0) | miss | merchant_gmv | 0.816 | 0 | -0.816 | 0.291 | -0.00538 |
| both | user | 2 (d=0) | FAR@0 | miss | — |  |  |  | -0.00133 | -0.000738 |
| both | all | FAR@0 | FAR@0 | miss | amount,channel,merchant_gmv | 0.691 | 0.0364 | 0.255 | 0.296 | 0.000566 |
