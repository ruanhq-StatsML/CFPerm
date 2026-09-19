# Feature library → FSDS unify → subset scan (order grain)

Default cut = **own-ref subset scan** (level set `{φ≥τ}` / coverage prefix). `merchant_id` / `user_id` lift Ŝ onto orders. No graph is built. Y is never a feature. Board: `library.html`. Brief for stakeholders: `report.html` / `REPORT.md`.

n_ref=480, n_new=160, merchants=16, batches=3, onset=1. Planted region is **south** (second half of merchant ids).

## Subset scan by grain (lift to orders)

| kind | t | onset | cut | loud south_frac | J(scan mer,south) | J(coverage mer,south) | J(mass mer,south) | J(user scan,south) | n_loud mer | FSDS selected | fingerprint |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | 0 |  | level_set/{phi>=tau} |  | 0.422 | 0.57 | 0.442 | 0 | 13 | amount,channel |  |
| covariate_south | 1 | yes | level_set/{phi>=tau} | 0.565 | 0.473 | 0.473 | 0.488 | 0 | 10 | channel,amount | x_shift |
| covariate_south | 2 | yes | level_set/{phi>=tau} | 0.852 | 0.767 | 0.885 | 0.885 | 0 | 7 | amount,channel,merchant_gmv | x_shift |
| concept_south | 0 |  | level_set/{phi>=tau} |  | 0.422 | 0.57 | 0.442 | 0 | 13 | amount,channel |  |
| concept_south | 1 | yes | level_set/{phi>=tau} |  | 0.31 | 0.38 | 0.348 | 0 | 8 | amount,n_items |  |
| concept_south | 2 | yes | level_set/{phi>=tau} | 0.521 | 0.455 | 0.41 | 0.54 | 0 | 10 | amount,channel | x_shift |
| both | 0 |  | level_set/{phi>=tau} |  | 0.422 | 0.57 | 0.442 | 0 | 13 | amount,channel |  |
| both | 1 | yes | level_set/{phi>=tau} | 0.605 | 0.558 | 0.488 | 0.488 | 0 | 11 | channel,amount | x_shift |
| both | 2 | yes | level_set/{phi>=tau} | 0.852 | 0.767 | 0.885 | 0.885 | 0 | 7 | amount,channel,merchant_gmv | x_shift |

## Feature library (native-grain catalog · last batch)

| kind | grain | feature | role | planted | selected | loud | score | mmd | cmean_x | cmean_y |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | amount | 订单金额（进 logit） | x | yes | yes | 6.74 | 0.101 | 0.935 | 0.167 |
| covariate_south | order | hour | 下单时刻（噪声列） |  |  |  | 0.709 | 0.00979 | 0.0908 | 0.0844 |
| covariate_south | order | n_items | 件数（噪声列） |  |  |  | 0.247 | -0.003 | 0.0617 | 0.00157 |
| covariate_south | order | channel | 渠道 | x | yes | yes | 2.34 | 0.035 | 0.555 | 0.169 |
| covariate_south | merchant | merchant_cat | 商户类目（进 logit，不种 shift） |  |  |  | 4.44e-16 | -0.0736 | 1.11e-16 | 0 |
| covariate_south | merchant | merchant_gmv | 商户 GMV | x | yes | yes | 3.77 | 0.0566 | 0.725 | 0 |
| covariate_south | merchant | n_skus | SKU 数（噪声列） |  |  |  | 4.44e-16 | -0.0719 | 1.11e-16 | 0 |
| covariate_south | user | user_tenure | 用户 tenure（进 logit，不种 shift） |  |  |  | 0.584 | -0.00146 | 0.0788 | 0.0584 |
| covariate_south | user | user_hist_freq | 历史频次（噪声列） |  |  |  | 0.124 | -0.00378 | 0.00194 | 0.0124 |
| concept_south | order | amount | 订单金额（进 logit） | y|x | yes |  | 1.7 | -0.00388 | 0.0556 | 0.17 |
| concept_south | order | hour | 下单时刻（噪声列） |  |  |  | 0.653 | 0.00979 | 0.0908 | 0.0271 |
| concept_south | order | n_items | 件数（噪声列） |  |  |  | 0.288 | -0.003 | 0.0617 | 0.0288 |
| concept_south | order | channel | 渠道 |  | yes |  | 0.97 | 0.00255 | 0.064 | 0.097 |
| concept_south | merchant | merchant_cat | 商户类目（进 logit，不种 shift） |  |  |  | 4.44e-16 | -0.0736 | 1.11e-16 | 0 |
| concept_south | merchant | merchant_gmv | 商户 GMV |  |  |  | 6.18e-17 | -0.0715 | 5.55e-17 | 0 |
| concept_south | merchant | n_skus | SKU 数（噪声列） |  |  |  | 4.44e-16 | -0.0719 | 1.11e-16 | 0 |
| concept_south | user | user_tenure | 用户 tenure（进 logit，不种 shift） |  |  |  | 1 | -0.00146 | 0.0788 | 0.1 |
| concept_south | user | user_hist_freq | 历史频次（噪声列） |  |  |  | 0.354 | -0.00378 | 0.00194 | 0.0354 |
| both | order | amount | 订单金额（进 logit） | both | yes | yes | 6.74 | 0.101 | 0.935 | 0.232 |
| both | order | hour | 下单时刻（噪声列） |  |  |  | 0.653 | 0.00979 | 0.0908 | 0.0323 |
| both | order | n_items | 件数（噪声列） |  |  |  | 0.5 | -0.003 | 0.0617 | 0.05 |
| both | order | channel | 渠道 | x | yes | yes | 2.34 | 0.035 | 0.555 | 0.119 |
| both | merchant | merchant_cat | 商户类目（进 logit，不种 shift） |  |  |  | 4.44e-16 | -0.0736 | 1.11e-16 | 0 |
| both | merchant | merchant_gmv | 商户 GMV | x | yes | yes | 3.77 | 0.0566 | 0.725 | 0 |
| both | merchant | n_skus | SKU 数（噪声列） |  |  |  | 4.44e-16 | -0.0719 | 1.11e-16 | 0 |
| both | user | user_tenure | 用户 tenure（进 logit，不种 shift） |  |  |  | 1.3 | -0.00146 | 0.0788 | 0.13 |
| both | user | user_hist_freq | 历史频次（噪声列） |  |  |  | 0.892 | -0.00378 | 0.00194 | 0.0892 |

## Own-ref vs full-ref (node clock) and multi-layer lift-to-order

| kind | t | south own MMD | north own MMD | south own ‖ΔX‖ | north own ‖ΔX‖ | J(scan mer,south) | J(coverage mer,south) | J(user scan,south) | MMD-slice south_frac | ΔY-slice south_frac |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | 1 | 0.0291 | -0.0059 | 0.929 | 0.716 | 0.473 | 0.473 | 0 | 1 | 0.556 |
| covariate_south | 2 | 0.257 | 0.016 | 2.28 | 0.676 | 0.767 | 0.885 | 0 | 1 | 0.75 |
| concept_south | 1 | -0.00763 | -0.0059 | 0.644 | 0.716 | 0.31 | 0.38 | 0 | 1 | 0.571 |
| concept_south | 2 | -0.0106 | 0.016 | 0.598 | 0.676 | 0.455 | 0.41 | 0 | 0.333 | 0.625 |
| both | 1 | 0.0291 | -0.0059 | 0.929 | 0.716 | 0.558 | 0.488 | 0 | 1 | 0.571 |
| both | 2 | 0.257 | 0.016 | 2.28 | 0.676 | 0.767 | 0.885 | 0 | 1 | 0.714 |

## Subset portraits (scan loud vs other · MMD + PO + CMean)

| kind | t | subset | n | south_frac | π_MMD | π_PO | π_CMean | MMD | PO | ΔE[Y] | ‖ΔE[X]‖ | fingerprint |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | 1 | loud | 108 | 0.565 | 1 | 0 | 0.969 | 0.0295 | 6.17e-05 | 0.0188 | 0.538 | x_shift |
| covariate_south | 1 | other | 52 | 0.404 | 0 | 0 | 0.0307 | -0.00286 | 0.000141 | 0.0147 | 0.0171 | x_shift |
| covariate_south | 2 | loud | 81 | 0.852 | 1 | 0 | 0.893 | 0.375 | 6.47e-05 | 0.261 | 2.34 | x_shift |
| covariate_south | 2 | other | 79 | 0.114 | 0 | 0 | 0.107 | -0.0027 | 0.000215 | 0.0263 | 0.325 | x_shift |
| concept_south | 1 | loud | 87 | 0.46 | 0 | 0 | 0 | 0.00602 | 8.59e-05 | 0.00725 | 0.112 | x_shift |
| concept_south | 1 | other | 73 | 0.575 | 0 | 0 | 0 | -0.0072 | 6.26e-05 | 0.0456 | 0.13 | x_shift |
| concept_south | 2 | loud | 117 | 0.521 | 0 | 0 | 0.712 | 0.00497 | 0.000704 | -0.069 | 0.154 | x_shift |
| concept_south | 2 | other | 43 | 0.395 | 0 | 0 | 0.288 | -0.00832 | 0.000537 | -0.0279 | 0.141 | x_shift |
| both | 1 | loud | 119 | 0.605 | 1 | 0 | 0.738 | 0.0294 | 0.000154 | 0.0423 | 0.552 | x_shift |
| both | 1 | other | 41 | 0.244 | 0 | 0 | 0.262 | -0.000303 | 9.59e-05 | 0.0192 | 0.196 | x_shift |
| both | 2 | loud | 81 | 0.852 | 1 | 0 | 0.87 | 0.375 | 0.000213 | -0.233 | 2.34 | x_shift |
| both | 2 | other | 79 | 0.114 | 0 | 0 | 0.13 | -0.0027 | 0.000286 | -0.037 | 0.325 | x_shift |

## FSDS ranking (last batch, top of each grain)

| kind | grain | feature | score | mmd | cmean_x | cmean_y | loud |
| --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | order | amount | 6.74 | 0.101 | 0.935 | 0.167 | yes |
| covariate_south | order | channel | 2.34 | 0.035 | 0.555 | 0.169 | yes |
| covariate_south | order | hour | 0.709 | 0.00979 | 0.0908 | 0.0844 |  |
| covariate_south | merchant | merchant_gmv | 3.77 | 0.0566 | 0.725 | 0 | yes |
| covariate_south | merchant | merchant_cat | 4.44e-16 | -0.0736 | 1.11e-16 | 0 |  |
| covariate_south | merchant | n_skus | 4.44e-16 | -0.0719 | 1.11e-16 | 0 |  |
| covariate_south | user | user_tenure | 0.584 | -0.00146 | 0.0788 | 0.0584 |  |
| covariate_south | user | user_hist_freq | 0.124 | -0.00378 | 0.00194 | 0.0124 |  |
| concept_south | order | amount | 1.7 | -0.00388 | 0.0556 | 0.17 |  |
| concept_south | order | channel | 0.97 | 0.00255 | 0.064 | 0.097 |  |
| concept_south | order | hour | 0.653 | 0.00979 | 0.0908 | 0.0271 |  |
| concept_south | merchant | merchant_cat | 4.44e-16 | -0.0736 | 1.11e-16 | 0 |  |
| concept_south | merchant | n_skus | 4.44e-16 | -0.0719 | 1.11e-16 | 0 |  |
| concept_south | merchant | merchant_gmv | 6.18e-17 | -0.0715 | 5.55e-17 | 0 |  |
| concept_south | user | user_tenure | 1 | -0.00146 | 0.0788 | 0.1 |  |
| concept_south | user | user_hist_freq | 0.354 | -0.00378 | 0.00194 | 0.0354 |  |
| both | order | amount | 6.74 | 0.101 | 0.935 | 0.232 | yes |
| both | order | channel | 2.34 | 0.035 | 0.555 | 0.119 | yes |
| both | order | hour | 0.653 | 0.00979 | 0.0908 | 0.0323 |  |
| both | merchant | merchant_gmv | 3.77 | 0.0566 | 0.725 | 0 | yes |
| both | merchant | merchant_cat | 4.44e-16 | -0.0736 | 1.11e-16 | 0 |  |
| both | merchant | n_skus | 4.44e-16 | -0.0719 | 1.11e-16 | 0 |  |
| both | user | user_tenure | 1.3 | -0.00146 | 0.0788 | 0.13 |  |
| both | user | user_hist_freq | 0.892 | -0.00378 | 0.00194 | 0.0892 |  |
