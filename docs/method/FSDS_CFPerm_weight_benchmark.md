# FSDS/CFPerm weight benchmark

FSDS RF-domain VIMP → membership → soft sample_w + hard observation subset_w on amazon/digits/news20/synth_cgcd. Pure membership is diagnostic; fused (membership ∪ private-live) is the cold-start weight channel.

## Aggregate (vs uniform)

| channel | mean Δ novel | mean Δ acc |
| --- | --- | --- |
| sample_w (membership) | -0.007 | — |
| subset_w (hard obs) | 0.031 | — |
| sample × subset | 0.117 | -0.017 |
| sample_w fused | — | 0.005 |
| subset_w fused | — | -0.001 |
| both fused | 0.165 | -0.040 |

## Per dataset

| dataset | n_ref/n_live | subset_frac | vimp_top | Δacc both | Δacc both_fused | Δnovel both_fused |
| --- | --- | --- | --- | --- | --- | --- |
| wine | 90/62 | 0.403 | 0.673 | 0.050 | 0.050 | 0.333 |
| covtype | 800/360 | 0.400 | 0.587 | -0.008 | -0.058 | 0.000 |
| news20b | 500/220 | 0.400 | 0.309 | -0.110 | -0.151 | 0.125 |
| cancer | 212/124 | 0.403 | 0.313 | 0.000 | 0.000 | 0.200 |

### Scheme detail (acc / novel / forgetting)

**wine** — UCI wine; known {0,1} → sparse inject {2}

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.900 | 0.667 | 0.941 | 0.022 | 0.000 |
| sample_w | 0.900 | 0.667 | 0.941 | 0.033 | 0.261 |
| subset_w | 0.900 | 0.667 | 0.941 | 0.044 | 0.304 |
| sample_x_subset | 0.950 | 1.000 | 0.941 | 0.056 | 0.432 |
| sample_w_fused | 0.900 | 0.667 | 0.941 | 0.033 | 0.261 |
| subset_w_fused | 0.900 | 0.667 | 0.941 | 0.044 | 0.304 |
| both_fused | 0.950 | 1.000 | 0.941 | 0.056 | 0.432 |

**covtype** — covertype subsample; known {1,2,3} → sparse {4,5}

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.592 | 0.524 | 0.606 | 0.363 | 0.000 |
| sample_w | 0.617 | 0.619 | 0.616 | 0.427 | 0.272 |
| subset_w | 0.608 | 0.571 | 0.616 | 0.421 | 0.316 |
| sample_x_subset | 0.583 | 0.333 | 0.636 | 0.396 | 0.450 |
| sample_w_fused | 0.600 | 0.619 | 0.596 | 0.441 | 0.258 |
| subset_w_fused | 0.575 | 0.619 | 0.566 | 0.476 | 0.317 |
| both_fused | 0.533 | 0.524 | 0.535 | 0.506 | 0.437 |

**news20b** — TF-IDF 4-topic; known atheism/religion → sparse med/hardware

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.630 | 0.125 | 0.692 | 0.286 | 0.000 |
| sample_w | 0.507 | 0.000 | 0.569 | 0.310 | 0.325 |
| subset_w | 0.562 | 0.000 | 0.631 | 0.332 | 0.309 |
| sample_x_subset | 0.521 | 0.250 | 0.554 | 0.318 | 0.467 |
| sample_w_fused | 0.521 | 0.375 | 0.538 | 0.376 | 0.311 |
| subset_w_fused | 0.521 | 0.250 | 0.554 | 0.402 | 0.311 |
| both_fused | 0.479 | 0.250 | 0.508 | 0.406 | 0.434 |

**cancer** — breast-cancer; REF benign → LIVE + sparse malignant

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.805 | 0.800 | 0.806 | 0.094 | 0.000 |
| sample_w | 0.927 | 0.800 | 0.944 | 0.052 | 0.280 |
| subset_w | 0.927 | 1.000 | 0.917 | 0.104 | 0.311 |
| sample_x_subset | 0.805 | 1.000 | 0.778 | 0.222 | 0.444 |
| sample_w_fused | 0.927 | 0.800 | 0.944 | 0.057 | 0.280 |
| subset_w_fused | 0.927 | 1.000 | 0.917 | 0.104 | 0.311 |
| both_fused | 0.805 | 1.000 | 0.778 | 0.222 | 0.442 |

Weights: `sample_w = 1+λ·membership`; `subset_w` boosts observations with membership ≥ q.

```bash
python3 scripts/run_fsds_cfperm_weight_benchmark.py --quick
```
