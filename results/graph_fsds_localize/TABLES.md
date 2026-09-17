# Graph localization → FSDS unify → two-layer subset (order grain)

Shares are localization proxies, not a unique decomposition. T is the batch label. Y is never a feature and never the subset key.

n_ref=360, n_new=100, batches=3, onset=1. Planted subset is always **south**.

## Per-batch order-grain portraits (MMD + PO-risk + Conditional Mean)

| kind | t | onset | FSDS selected | loud subset | south π_MMD | south π_PO | south π_CMean | south mix | south MMD | south PO | south ΔE[Y] | south ‖ΔE[X]‖ | gap MMD vs other | fingerprint | LOGO π_MMD o/m/u |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | 0 |  | hour,channel | south | 0 | 0 | 0.742 | 0/0/1 | -0.00405 | 6.41e-05 | -0.0288 | 0.259 | -0.00804 | x_shift | 0/0/0 |
| covariate_south | 1 | yes | amount,hour | south | 1 | 0 | 0.816 | 0.551/0/0.449 | 0.0226 | 0.000109 | 0.172 | 0.575 | 0.0257 | x_shift | 0/0/0 |
| covariate_south | 2 | yes | amount,channel,merchant_gmv | south | 1 | 0 | 0.789 | 0.559/0/0.441 | 0.421 | 2.85e-05 | 0.325 | 2.45 | 0.45 | x_shift | 0/1/0 |
| concept_south | 0 |  | hour,channel | south | 0 | 0 | 0.742 | 0/0/1 | -0.00405 | 6.41e-05 | -0.0288 | 0.259 | -0.00804 | x_shift | 0/0/0 |
| concept_south | 1 | yes | amount,channel | south | 0 | 0 | 0.687 | 0/0/1 | -0.00247 | 0.000796 | -0.0832 | 0.324 | 0.00166 | x_shift | 0/0/0 |
| concept_south | 2 | yes | amount | south | 0 | 0.665 | 0.691 | 0/0.49/0.51 | 0.0157 | 0.00158 | -0.259 | 0.278 | 0.0151 | concept | 0/0/0 |
| both | 0 |  | hour,channel | south | 0 | 0 | 0.742 | 0/0/1 | -0.00405 | 6.41e-05 | -0.0288 | 0.259 | -0.00804 | x_shift | 0/0/0 |
| both | 1 | yes | n_items,amount | south | 1 | 0.943 | 0.716 | 0.376/0.355/0.269 | 0.0247 | 0.00118 | -0.0407 | 0.575 | 0.0265 | both | 0/0/0 |
| both | 2 | yes | amount,channel,merchant_gmv | south | 1 | 0 | 0.78 | 0.562/0/0.438 | 0.421 | 0.000859 | -0.3 | 2.45 | 0.45 | x_shift | 0/1/0 |

## FSDS ranking (last batch, top of each grain)

| kind | grain | feature | score | mmd | cmean_x | cmean_y | loud |
| --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | amount | 5.35 | 0.0802 | 0.774 | 0.156 | yes |
| covariate_south | order | channel | 3.13 | 0.047 | 0.685 | 0.136 | yes |
| covariate_south | order | hour | 1.59 | 0.00405 | 0.106 | 0.159 |  |
| covariate_south | merchant | merchant_gmv | 2.06 | -0.0418 | 0.722 | 0 | yes |
| covariate_south | merchant | n_skus | 1.41e-15 | -0.143 | 5.55e-16 | 0 |  |
| covariate_south | merchant | merchant_cat | 1.11e-16 | -0.152 | 2.78e-17 | 0 |  |
| covariate_south | user | user_hist_freq | 0.799 | -7.65e-06 | 0.0712 | 0.0799 |  |
| covariate_south | user | user_tenure | 0.581 | 0.00191 | 0.145 | 0.0429 |  |
| concept_south | order | amount | 3.05 | 0.0062 | 0.0917 | 0.305 | yes |
| concept_south | order | channel | 0.789 | 0.00508 | 0.122 | 0.0789 |  |
| concept_south | order | hour | 0.425 | 0.00405 | 0.106 | 0.0411 |  |
| concept_south | merchant | n_skus | 1.41e-15 | -0.143 | 5.55e-16 | 0 |  |
| concept_south | merchant | merchant_gmv | 6.35e-16 | -0.145 | 2.22e-16 | 0 |  |
| concept_south | merchant | merchant_cat | 1.11e-16 | -0.152 | 2.78e-17 | 0 |  |
| concept_south | user | user_hist_freq | 0.601 | -7.65e-06 | 0.0712 | 0.0601 |  |
| concept_south | user | user_tenure | 0.581 | 0.00191 | 0.145 | 0.0293 |  |
| both | order | amount | 5.51 | 0.0802 | 0.774 | 0.551 | yes |
| both | order | channel | 3.47 | 0.047 | 0.685 | 0.347 | yes |
| both | order | n_items | 1.46 | -0.00287 | 0.0635 | 0.146 |  |
| both | merchant | merchant_gmv | 2.06 | -0.0418 | 0.722 | 0 | yes |
| both | merchant | n_skus | 1.41e-15 | -0.143 | 5.55e-16 | 0 |  |
| both | merchant | merchant_cat | 1.11e-16 | -0.152 | 2.78e-17 | 0 |  |
| both | user | user_hist_freq | 0.675 | -7.65e-06 | 0.0712 | 0.0675 |  |
| both | user | user_tenure | 0.581 | 0.00191 | 0.145 | 0.0638 |  |
