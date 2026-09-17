# Recsys order-stream by grain — OnlineRFPerm × RFPerm/CFPerm × FSDS

Data: synthetic order / merchant / user stream (recommendation-shaped table). Y = conversion, never a feature. T = batch label. OnlineRFPerm (Algorithm 1) marks WHEN. RFPerm ΔMSE and CFPerm/PermuCATE VIMP plus FSDS mark WHICH columns; Kendall-τ is the ranking recovery in the MetaLearner paper. Post-hoc region split marks WHICH accounts (south vs north). Localization, not a unique decomposition. Not a graph method.

n_ref=400, n_new=120, batches=6, onset_true=2, merchants=16. Planted subset = south. Delay = first rejection − onset_true. Negative / FAR = mark before labeled onset.

## Table 1 · OnlineRFPerm (WHEN) — first rejection and delay

| kind | dim | onset_true | hop | rank-p (p<0.05) | ADDIS | SAFFRON | last T | last MMD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | 2 | miss | 2 (d=0) | 2 (d=0) | 2 (d=0) | 0.067 | 0.0952 |
| covariate_south | merchant | 2 | miss | 2 (d=0) | 5 (d=3) | 5 (d=3) | 0.066 | 0.0645 |
| covariate_south | user | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | 0.0569 | 0.000627 |
| covariate_south | all | 2 | miss | FAR@0 | FAR@0 | FAR@0 | 0.0531 | 0.0773 |
| concept_south | order | 2 | miss | 2 (d=0) | 2 (d=0) | 2 (d=0) | 0.113 | -0.000532 |
| concept_south | merchant | 2 | miss | 2 (d=0) | miss | miss | 0.0188 | -0.00488 |
| concept_south | user | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | 0.00947 | 0.000627 |
| concept_south | all | 2 | miss | FAR@0 | FAR@0 | FAR@0 | 0.11 | -0.00121 |
| both | order | 2 | miss | 2 (d=0) | 2 (d=0) | 2 (d=0) | 0.157 | 0.0952 |
| both | merchant | 2 | miss | 2 (d=0) | miss | miss | -0.0187 | 0.0645 |
| both | user | 2 | miss | FAR@0 | 2 (d=0) | 2 (d=0) | -0.000541 | 0.000627 |
| both | all | 2 | miss | FAR@0 | FAR@0 | FAR@0 | 0.17 | 0.0773 |

## Table 2 · Ranking recovery (WHICH columns) — Kendall-τ vs planted magnitude

Planted order on covariate/both: amount ≻ merchant_gmv ≻ channel ≻ noise. Concept plants amount in Y|X only. τ is Kendall’s τ between scores and that magnitude. CFPerm reject = max φ-VIMP vs 95% of T-permuted nulls (B=12). Empty reject is a miss, not a quiet stream.

| kind | dim | planted in grain | FSDS recovered | τ_FSDS | RFPerm ΔMSE top-3 | τ_RFPerm | CFPerm φ top-3 | τ_CFPerm | CFPerm reject |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | amount,channel | amount,channel | 0.913 | amount,n_items,channel | 0.548 | n_items,channel,amount | -0.183 |  |
| covariate_south | merchant | merchant_gmv | merchant_gmv | 0.816 | merchant_cat,n_skus,merchant_gmv | -0.816 | merchant_cat,n_skus,merchant_gmv | -0.816 |  |
| covariate_south | user | — | — |  | user_tenure,user_hist_freq |  | user_hist_freq,user_tenure |  |  |
| covariate_south | all | amount,channel,merchant_gmv | amount,channel,merchant_gmv | 0.691 | amount,merchant_cat,user_hist_freq | 0.182 | n_items,channel,user_hist_freq | -0.182 |  |
| concept_south | order | amount | amount | 0.236 | channel,amount,hour | 0.236 | amount,hour,channel | 0.707 |  |
| concept_south | merchant | — | — |  | merchant_cat,merchant_gmv,n_skus |  | n_skus,merchant_cat,merchant_gmv |  |  |
| concept_south | user | — | — |  | user_hist_freq,user_tenure |  | user_hist_freq,user_tenure |  |  |
| concept_south | all | amount | amount | 0.354 | merchant_cat,amount,merchant_gmv | 0.354 | amount,n_items,user_hist_freq | 0.471 |  |
| both | order | amount,channel | amount,channel | 0.913 | channel,hour,n_items | -0.183 | amount,channel,n_items | 0.913 |  |
| both | merchant | merchant_gmv | merchant_gmv | 0.816 | merchant_cat,merchant_gmv,n_skus | 0 | n_skus,merchant_cat,merchant_gmv | -0.816 |  |
| both | user | — | — |  | user_hist_freq,user_tenure |  | user_tenure,user_hist_freq |  |  |
| both | all | amount,channel,merchant_gmv | amount,channel,merchant_gmv | 0.691 | merchant_cat,merchant_gmv,channel | 0.0364 | channel,user_hist_freq,amount | 0.255 |  |

## Table 3 · Post-hoc localization (south vs north, last batch)

| kind | dim | n_south | south MMD vs own-ref | north MMD vs own-ref | pair MMD | south ΔE[Y] | north ΔE[Y] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | 59 | 0.35 | 0.00109 | 0.316 | 0.349 | -0.0489 |
| covariate_south | merchant | 59 | 0.291 | -0.00538 | 0.286 | 0.349 | -0.0489 |
| covariate_south | user | 59 | -0.00133 | -0.000738 | -0.00304 | 0.349 | -0.0489 |
| covariate_south | all | 59 | 0.296 | 0.000566 | 0.261 | 0.349 | -0.0489 |
| concept_south | order | 59 | -0.00137 | 0.00109 | 0.00159 | -0.0412 | -0.0489 |
| concept_south | merchant | 59 | -0.00969 | -0.00538 | 0.179 | -0.0412 | -0.0489 |
| concept_south | user | 59 | -0.00133 | -0.000738 | -0.00304 | -0.0412 | -0.0489 |
| concept_south | all | 59 | -0.00413 | 0.000566 | 0.0325 | -0.0412 | -0.0489 |
| both | order | 59 | 0.35 | 0.00109 | 0.316 | -0.262 | -0.0489 |
| both | merchant | 59 | 0.291 | -0.00538 | 0.286 | -0.262 | -0.0489 |
| both | user | 59 | -0.00133 | -0.000738 | -0.00304 | -0.262 | -0.0489 |
| both | all | 59 | 0.296 | 0.000566 | 0.261 | -0.262 | -0.0489 |

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
