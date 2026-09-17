# Graph localization → FSDS unify → two-layer subset (order grain)

Package: **networkx** `louvain_communities`. Not GraphRAG, not PyG.
Default cut = **bundled-shift** Louvain (changing-subset objective). Structural Louvain is reported only as a contrast — modularity of co-order/kNN is not the shift object.
Shares are localization proxies, not a unique decomposition. Y is never a feature.

n_ref=480, n_new=160, merchants=16, batches=3, onset=1. Planted region is **south** (second half of merchant ids).

## Graph cuts (networkx Louvain, two objectives)

| kind | t | onset | package | bundled loud | bundled south_frac | layer2 loud | bundled n_comm | bundled n_edges | struct n_comm | FSDS selected | fingerprint |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | 0 |  | networkx | C0 | 0.462 |  | 1 | 77 | 4 | amount,channel |  |
| covariate_south | 1 | yes | networkx | C0 | 0.6 | C0 | 5 | 45 | 3 | channel,amount | x_shift |
| covariate_south | 2 | yes | networkx | C0 | 0.857 | C0 | 6 | 21 | 3 | amount,channel,merchant_gmv | x_shift |
| concept_south | 0 |  | networkx | C0 | 0.462 |  | 1 | 77 | 4 | amount,channel |  |
| concept_south | 1 | yes | networkx | C0 | 0.5 |  | 7 | 28 | 4 | amount,n_items |  |
| concept_south | 2 | yes | networkx | C1 | 0.25 | C0 | 4 | 45 | 4 | amount,channel | concept |
| both | 0 |  | networkx | C0 | 0.462 |  | 1 | 77 | 4 | amount,channel |  |
| both | 1 | yes | networkx | C0 | 0.636 | C0 | 4 | 55 | 3 | channel,amount | x_shift |
| both | 2 | yes | networkx | C0 | 0.857 | C0 | 6 | 21 | 3 | amount,channel,merchant_gmv | x_shift |

## Community portraits (bundled cut · MMD + PO + CMean)

| kind | t | community | n | south_frac | π_MMD | π_PO | π_CMean | MMD | PO | ΔE[Y] | ‖ΔE[X]‖ | fingerprint |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| covariate_south | 1 | C0 | 108 | 0.6 | 1 | 0 | 1 | 0.0295 | 6.17e-05 | 0.0188 | 0.538 | x_shift |
| covariate_south | 2 | C0 | 81 | 0.857 | 0.797 | 0 | 0.667 | 0.375 | 6.47e-05 | 0.261 | 2.34 | x_shift |
| covariate_south | 2 | C? | 22 | 0 | 0.203 | 0 | 0.333 | 0.0954 |  | 0.118 | 1.28 | x_shift |
| concept_south | 1 | C0 | 87 | 0.5 | 0 | 0 | 0 | 0.00602 | 8.59e-05 | 0.00725 | 0.112 | x_shift |
| concept_south | 2 | C0 | 72 | 0.667 | 0 | 0.951 | 0.344 | -0.00386 | 0.00168 | -0.159 | 0.149 | concept |
| concept_south | 2 | C1 | 45 | 0.25 | 0.718 | 0.0489 | 0.406 | 0.0388 | 8.63e-05 | 0.0698 | 0.592 | x_shift |
| concept_south | 2 | C? | 22 | 0 | 0.282 | 0 | 0.251 | 0.0153 |  | -0.0634 | 0.294 | x_shift |
| both | 1 | C0 | 119 | 0.636 | 1 | 0 | 1 | 0.0294 | 0.000154 | 0.0423 | 0.552 | x_shift |
| both | 2 | C0 | 81 | 0.857 | 0.797 | 0 | 0.663 | 0.375 | 0.000213 | -0.233 | 2.34 | x_shift |
| both | 2 | C? | 22 | 0 | 0.203 | 0 | 0.337 | 0.0954 |  | -0.109 | 1.28 | x_shift |

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
