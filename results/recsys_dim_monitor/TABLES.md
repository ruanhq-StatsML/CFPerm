# Recsys order-stream by grain — OnlineRFPerm × RFPerm/CFPerm × FSDS

Data: synthetic order / merchant / user stream (recommendation-shaped table). Y = conversion, never a feature. T = batch label. OnlineRFPerm (Algorithm 1) marks WHEN. RFPerm ΔMSE and CFPerm/PermuCATE VIMP plus FSDS mark WHICH columns; Kendall-τ is the ranking recovery in the MetaLearner paper. Post-hoc region split marks WHICH accounts (south vs north). Localization, not a unique decomposition. Not a graph method.

n_ref=400, n_new=120, batches=6, onset_true=2, merchants=16. Planted subset = south. Delay = first rejection − onset_true. Negative / FAR = mark before labeled onset.

## Table 1 · OnlineRFPerm (WHEN) — first rejection and delay

| kind | dim | onset_true | hop | rank-p (p<0.05) | ADDIS | SAFFRON | last T | last MMD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | 2 | miss | 2 (d=0) | 2 (d=0) | 2 (d=0) | 0.067 | 0.0952 |
| covariate_south | merchant | 2 | miss | FAR@1 | 2 (d=0) | 2 (d=0) | 0.131 | 0.121 |
| covariate_south | user | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | 0.0569 | 0.000627 |
| covariate_south | all | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | 0.0602 | 0.0865 |
| concept_south | order | 2 | miss | 2 (d=0) | 2 (d=0) | 2 (d=0) | 0.113 | -0.000532 |
| concept_south | merchant | 2 | miss | FAR@1 | 2 (d=0) | 2 (d=0) | 0.0338 | 0.00216 |
| concept_south | user | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | 0.00947 | 0.000627 |
| concept_south | all | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | 0.116 | -0.000467 |
| both | order | 2 | miss | 2 (d=0) | 2 (d=0) | 2 (d=0) | 0.157 | 0.0952 |
| both | merchant | 2 | miss | FAR@1 | 2 (d=0) | 2 (d=0) | -0.0165 | 0.121 |
| both | user | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | -0.000541 | 0.000627 |
| both | all | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | 0.134 | 0.0865 |

## Table 2 · Ranking recovery (WHICH columns) — Kendall-τ vs planted magnitude

Planted order on covariate/both: amount ≻ merchant_gmv ≻ channel ≻ noise. Concept plants amount in Y|X only. τ is Kendall’s τ between scores and that magnitude. CFPerm reject = max φ-VIMP vs 95% of T-permuted nulls (B=12). Empty reject is a miss, not a quiet stream.

| kind | dim | planted in grain | FSDS recovered | τ_FSDS | RFPerm ΔMSE top-3 | τ_RFPerm | CFPerm φ top-3 | τ_CFPerm | CFPerm reject |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | amount,channel | amount,channel | 0.913 | amount,n_items,channel | 0.548 | n_items,channel,amount | -0.183 |  |
| covariate_south | merchant | merchant_gmv | merchant_gmv | 0.816 | merchant_cat,n_skus,merchant_gmv | -0.816 | merchant_gmv,n_skus,merchant_cat | 0.816 |  |
| covariate_south | user | — | — |  | user_tenure,user_hist_freq |  | user_hist_freq,user_tenure |  |  |
| covariate_south | all | amount,channel,merchant_gmv | merchant_gmv,amount,channel | 0.691 | amount,merchant_cat,n_items | 0.0364 | n_skus,n_items,user_hist_freq | -0.109 |  |
| concept_south | order | amount | amount | 0.236 | channel,amount,hour | 0.236 | amount,hour,channel | 0.707 |  |
| concept_south | merchant | — | — |  | merchant_gmv,merchant_cat,n_skus |  | merchant_gmv,merchant_cat,n_skus |  |  |
| concept_south | user | — | — |  | user_hist_freq,user_tenure |  | user_hist_freq,user_tenure |  |  |
| concept_south | all | amount | amount | 0.354 | amount,merchant_cat,merchant_gmv | 0.471 | amount,n_items,hour | 0.471 |  |
| both | order | amount,channel | amount,channel | 0.913 | channel,hour,n_items | -0.183 | amount,channel,n_items | 0.913 |  |
| both | merchant | merchant_gmv | merchant_gmv | 0.816 | merchant_gmv,merchant_cat,n_skus | 0.816 | merchant_gmv,merchant_cat,n_skus | 0.816 |  |
| both | user | — | — |  | user_hist_freq,user_tenure |  | user_tenure,user_hist_freq |  |  |
| both | all | amount,channel,merchant_gmv | merchant_gmv,amount,channel | 0.691 | merchant_gmv,merchant_cat,hour | -0.109 | user_hist_freq,channel,n_items | -0.0364 |  |

## Table 3 · Post-hoc localization (south vs north, last batch)

Own-ref clock: this region's new bag vs this region's D_ref. Three readouts together: MMD (P(X)), CMean (||ΔE[X]|| and ΔE[Y]), PO-risk (P(Y|X)). Subset key is region, not Y.

| kind | dim | south MMD | south ‖ΔE[X]‖ | south PO | south ΔE[Y] | north MMD | north ‖ΔE[X]‖ | north PO | north ΔE[Y] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | 0.35 | 2.55 | 7.59e-05 | 0.349 | 0.00109 | 0.302 | 0.000898 | -0.0489 |
| covariate_south | merchant | 0.465 | 1.45 | 0.00107 | 0.349 | 0.00333 | 0.212 | 0.000141 | -0.0489 |
| covariate_south | user | -0.00133 | 0.0659 | 0.00329 | 0.349 | -0.000738 | 0.189 | 0.000486 | -0.0489 |
| covariate_south | all | 0.323 | 2.94 | 0.000263 | 0.349 | 0.00155 | 0.415 | 0.000828 | -0.0489 |
| concept_south | order | -0.00137 | 0.275 | 0.00307 | -0.0412 | 0.00109 | 0.302 | 0.000898 | -0.0489 |
| concept_south | merchant | 0.00546 | 0.063 | 0.00029 | -0.0412 | 0.00333 | 0.212 | 0.000141 | -0.0489 |
| concept_south | user | -0.00133 | 0.0659 | 0.000132 | -0.0412 | -0.000738 | 0.189 | 0.000486 | -0.0489 |
| concept_south | all | -0.00302 | 0.29 | 0.00295 | -0.0412 | 0.00155 | 0.415 | 0.000828 | -0.0489 |
| both | order | 0.35 | 2.55 | 0.0015 | -0.262 | 0.00109 | 0.302 | 0.000898 | -0.0489 |
| both | merchant | 0.465 | 1.45 | 7.9e-06 | -0.262 | 0.00333 | 0.212 | 0.000141 | -0.0489 |
| both | user | -0.00133 | 0.0659 | 0.00183 | -0.262 | -0.000738 | 0.189 | 0.000486 | -0.0489 |
| both | all | 0.323 | 2.94 | 0.000139 | -0.262 | 0.00155 | 0.415 | 0.000828 | -0.0489 |

## Table 4 · Sequential T_t / p_t on the order grain (OnlineRFPerm Algorithm 1)

| kind | t | onset | T = MSE−E_ref | p | hop | MMD | ΔE[Y] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | 0 |  | 0.0274 | 0.0547 |  | 0.0024 | -0.0375 |
| covariate_south | 1 |  | 0.0233 | 0.0796 |  | 0.00252 | 0.0208 |
| covariate_south | 2 | yes | 0.0545 | 0.00498 |  | -0.00198 | 0.0792 |
| covariate_south | 3 | yes | 0.0775 | 0.00498 |  | 0.031 | 0.0542 |
| covariate_south | 4 | yes | 0.0479 | 0.00498 |  | 0.0498 | 0.104 |
| covariate_south | 5 | yes | 0.067 | 0.00498 |  | 0.0952 | 0.146 |
| concept_south | 0 |  | 0.0274 | 0.0547 |  | 0.0024 | -0.0375 |
| concept_south | 1 |  | 0.0233 | 0.0796 |  | 0.00252 | 0.0208 |
| concept_south | 2 | yes | 0.0526 | 0.00498 |  | -0.00321 | 0.0375 |
| concept_south | 3 | yes | 0.0658 | 0.00498 |  | 0.00227 | -0.0708 |
| concept_south | 4 | yes | 0.0926 | 0.00498 |  | 0.000621 | -0.0375 |
| concept_south | 5 | yes | 0.113 | 0.00498 |  | -0.000532 | -0.0458 |
| both | 0 |  | 0.0274 | 0.0547 |  | 0.0024 | -0.0375 |
| both | 1 |  | 0.0233 | 0.0796 |  | 0.00252 | 0.0208 |
| both | 2 | yes | 0.0501 | 0.00498 |  | -0.00198 | 0.0708 |
| both | 3 | yes | 0.0983 | 0.00498 |  | 0.031 | -0.0542 |
| both | 4 | yes | 0.105 | 0.00498 |  | 0.0498 | -0.0708 |
| both | 5 | yes | 0.157 | 0.00498 |  | 0.0952 | -0.154 |
