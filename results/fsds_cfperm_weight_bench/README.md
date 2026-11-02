# FSDS/CFPerm weight benchmark

FSDS RF-domain VIMP → membership → soft sample_w + hard observation subset_w on har/mice/cnae9/jvowel (non-routine category shift). Primary metric = novel_acc on injected new categories; fused (membership ∪ private-live) is the cold-start weight channel.

## Aggregate (vs uniform)

| channel | mean Δ novel | mean Δ acc |
| --- | --- | --- |
| sample_w (membership) | 0.058 | — |
| subset_w (hard obs) | 0.058 | — |
| sample × subset | -0.065 | -0.043 |
| sample_w fused | — | -0.024 |
| subset_w fused | — | -0.007 |
| both fused | 0.090 | -0.054 |

## Per dataset

| dataset | n_ref/n_live | subset_frac | vimp_top | Δacc both | Δacc both_fused | Δnovel both_fused |
| --- | --- | --- | --- | --- | --- | --- |
| har | 800/375 | 0.400 | 0.132 | -0.064 | -0.064 | -0.045 |
| mice | 320/208 | 0.399 | 0.195 | -0.130 | -0.101 | 0.167 |
| cnae9 | 350/228 | 0.399 | 0.407 | -0.039 | -0.092 | 0.188 |
| jvowel | 600/300 | 0.400 | 0.505 | 0.060 | 0.040 | 0.053 |

### Scheme detail (acc / novel / forgetting)

**har** — HAR smartphone; locomotion {1,2,3} → sparse posture {4,5,6}

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.864 | 0.636 | 0.913 | 0.065 | 0.000 |
| sample_w | 0.824 | 0.682 | 0.854 | 0.142 | 0.293 |
| subset_w | 0.880 | 0.955 | 0.864 | 0.104 | 0.309 |
| sample_x_subset | 0.800 | 0.682 | 0.825 | 0.134 | 0.465 |
| sample_w_fused | 0.776 | 0.682 | 0.796 | 0.137 | 0.293 |
| subset_w_fused | 0.880 | 0.955 | 0.864 | 0.104 | 0.309 |
| both_fused | 0.800 | 0.591 | 0.845 | 0.142 | 0.465 |

**mice** — MiceProtein; control c-* → sparse trisomic t-* classes

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.812 | 0.667 | 0.863 | 0.128 | 0.000 |
| sample_w | 0.783 | 0.667 | 0.824 | 0.116 | 0.278 |
| subset_w | 0.696 | 0.444 | 0.784 | 0.134 | 0.317 |
| sample_x_subset | 0.681 | 0.278 | 0.824 | 0.134 | 0.448 |
| sample_w_fused | 0.768 | 0.722 | 0.784 | 0.134 | 0.269 |
| subset_w_fused | 0.768 | 0.667 | 0.804 | 0.147 | 0.317 |
| both_fused | 0.710 | 0.833 | 0.667 | 0.200 | 0.437 |

**cnae9** — CNAE-9 company docs; known cats 1–5 → sparse novel 6–9

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.803 | 0.438 | 0.900 | 0.051 | 0.000 |
| sample_w | 0.763 | 0.625 | 0.800 | 0.086 | 0.329 |
| subset_w | 0.763 | 0.625 | 0.800 | 0.089 | 0.314 |
| sample_x_subset | 0.763 | 0.625 | 0.800 | 0.089 | 0.455 |
| sample_w_fused | 0.789 | 0.625 | 0.833 | 0.097 | 0.327 |
| subset_w_fused | 0.763 | 0.625 | 0.800 | 0.097 | 0.312 |
| both_fused | 0.711 | 0.625 | 0.733 | 0.111 | 0.439 |

**jvowel** — JapaneseVowels; known speakers 1–5 → sparse novel 6–9

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.780 | 0.579 | 0.827 | 0.128 | 0.000 |
| sample_w | 0.820 | 0.579 | 0.877 | 0.100 | 0.273 |
| subset_w | 0.760 | 0.526 | 0.815 | 0.145 | 0.317 |
| sample_x_subset | 0.840 | 0.474 | 0.926 | 0.077 | 0.449 |
| sample_w_fused | 0.830 | 0.684 | 0.864 | 0.117 | 0.239 |
| subset_w_fused | 0.820 | 0.632 | 0.864 | 0.158 | 0.315 |
| both_fused | 0.820 | 0.632 | 0.864 | 0.167 | 0.420 |

Weights: `sample_w = 1+λ·membership`; `subset_w` boosts observations with membership ≥ q.

```bash
python3 scripts/run_fsds_cfperm_weight_benchmark.py --quick
```
