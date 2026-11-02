# Downsample vs FSDS sample-weight ablation

Label downsample (ds_known) is a blunt rebalance; FSDS sample weights are selective reweighting. Incremental value = best_weight − ds_known. Need downsample only when weights cannot recover novel mass and label balance itself is the bottleneck.

## Aggregate (novel acc deltas)

| contrast | mean Δ novel |
| --- | --- |
| ds_known − uniform | 0.127 |
| ds_random − uniform | -0.023 |
| ds_mem − uniform | 0.041 |
| best_weight − uniform | 0.178 |
| best_weight − ds_known (**incremental**) | 0.052 |
| ds_known+sample_f − ds_known | 0.021 |

Weights beat ds_known on **5 / 9** boards.

## Per board

| dataset | uniform | ds_known | ds_mem | best_w | best−ds | knownΔ@best | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wine | 0.667 | 1.000 | 0.667 | 1.000 (both_fused) | 0.000 | 0.000 | ds enough |
| cancer | 0.800 | 0.800 | 1.000 | 1.000 (subset_w_fused) | 0.200 | 0.111 | weights>ds |
| har | 0.636 | 0.818 | 0.818 | 0.955 (subset_w_fused) | 0.136 | -0.049 | weights>ds |
| mice | 0.667 | 0.833 | 0.500 | 0.833 (both_fused) | 0.000 | -0.196 | ds enough |
| cnae9 | 0.438 | 0.562 | 0.562 | 0.625 (sample_w_fused) | 0.062 | -0.067 | weights>ds |
| jvowel | 0.579 | 0.579 | 0.684 | 0.684 (sample_w_fused) | 0.105 | 0.037 | weights>ds |
| covtype | 0.524 | 0.857 | 0.619 | 0.619 (sample_w_fused) | -0.238 | -0.010 | ds enough |
| news20b | 0.125 | 0.125 | 0.000 | 0.375 (sample_w_fused) | 0.250 | -0.154 | weights>ds |
| isolet | 0.350 | 0.350 | 0.300 | 0.300 (both_fused) | -0.050 | -0.104 | neither |

### Scheme detail

**wine** — UCI wine; known {0,1} → sparse inject {2}

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.667 | 0.941 | 0.022 | 42 | 0.214 |
| ds_known | 1.000 | 0.941 | 0.044 | 18 | 0.500 |
| ds_random | 0.667 | 0.882 | 0.056 | 18 | 0.278 |
| ds_mem | 0.667 | 0.941 | 0.044 | 18 | 0.500 |
| sample_w_fused | 0.667 | 0.941 | 0.033 | 42 | 0.214 |
| subset_w_fused | 0.667 | 0.941 | 0.044 | 42 | 0.214 |
| both_fused | 1.000 | 0.941 | 0.056 | 42 | 0.214 |
| ds_known+sample_f | 1.000 | 1.000 | 0.022 | 18 | 0.500 |

**cancer** — breast-cancer; REF benign → LIVE + sparse malignant

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.800 | 0.806 | 0.094 | 83 | 0.229 |
| ds_known | 0.800 | 0.889 | 0.066 | 38 | 0.500 |
| ds_random | 0.600 | 0.889 | 0.061 | 38 | 0.263 |
| ds_mem | 1.000 | 0.833 | 0.208 | 38 | 0.500 |
| sample_w_fused | 0.800 | 0.944 | 0.057 | 83 | 0.229 |
| subset_w_fused | 1.000 | 0.917 | 0.104 | 83 | 0.229 |
| both_fused | 1.000 | 0.778 | 0.222 | 83 | 0.229 |
| ds_known+sample_f | 0.600 | 0.806 | 0.123 | 38 | 0.500 |

**har** — HAR smartphone; locomotion {1,2,3} → sparse posture {4,5,6}

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.636 | 0.913 | 0.065 | 250 | 0.212 |
| ds_known | 0.818 | 0.864 | 0.114 | 106 | 0.500 |
| ds_random | 0.818 | 0.942 | 0.054 | 106 | 0.255 |
| ds_mem | 0.818 | 0.767 | 0.284 | 106 | 0.500 |
| sample_w_fused | 0.682 | 0.796 | 0.137 | 250 | 0.212 |
| subset_w_fused | 0.955 | 0.864 | 0.104 | 250 | 0.212 |
| both_fused | 0.591 | 0.845 | 0.142 | 250 | 0.212 |
| ds_known+sample_f | 0.818 | 0.699 | 0.266 | 106 | 0.500 |

**mice** — MiceProtein; control c-* → sparse trisomic t-* classes

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.667 | 0.863 | 0.128 | 139 | 0.216 |
| ds_known | 0.833 | 0.725 | 0.241 | 60 | 0.500 |
| ds_random | 0.389 | 0.765 | 0.116 | 60 | 0.217 |
| ds_mem | 0.500 | 0.784 | 0.166 | 60 | 0.383 |
| sample_w_fused | 0.722 | 0.784 | 0.134 | 139 | 0.216 |
| subset_w_fused | 0.667 | 0.804 | 0.147 | 139 | 0.216 |
| both_fused | 0.833 | 0.667 | 0.200 | 139 | 0.216 |
| ds_known+sample_f | 0.833 | 0.686 | 0.406 | 60 | 0.500 |

**cnae9** — CNAE-9 company docs; known cats 1–5 → sparse novel 6–9

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.438 | 0.900 | 0.051 | 152 | 0.211 |
| ds_known | 0.562 | 0.783 | 0.103 | 64 | 0.500 |
| ds_random | 0.562 | 0.867 | 0.057 | 64 | 0.281 |
| ds_mem | 0.562 | 0.867 | 0.071 | 64 | 0.469 |
| sample_w_fused | 0.625 | 0.833 | 0.097 | 152 | 0.211 |
| subset_w_fused | 0.625 | 0.800 | 0.097 | 152 | 0.211 |
| both_fused | 0.625 | 0.733 | 0.111 | 152 | 0.211 |
| ds_known+sample_f | 0.562 | 0.767 | 0.106 | 64 | 0.500 |

**jvowel** — JapaneseVowels; known speakers 1–5 → sparse novel 6–9

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.579 | 0.827 | 0.128 | 200 | 0.205 |
| ds_known | 0.579 | 0.605 | 0.300 | 82 | 0.500 |
| ds_random | 0.684 | 0.778 | 0.240 | 82 | 0.232 |
| ds_mem | 0.684 | 0.901 | 0.103 | 82 | 0.171 |
| sample_w_fused | 0.684 | 0.864 | 0.117 | 200 | 0.205 |
| subset_w_fused | 0.632 | 0.864 | 0.158 | 200 | 0.205 |
| both_fused | 0.632 | 0.864 | 0.167 | 200 | 0.205 |
| ds_known+sample_f | 0.579 | 0.593 | 0.365 | 82 | 0.500 |

**covtype** — covertype subsample; known {1,2,3} → sparse {4,5}

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.524 | 0.606 | 0.363 | 240 | 0.163 |
| ds_known | 0.857 | 0.465 | 0.488 | 78 | 0.500 |
| ds_random | 0.429 | 0.535 | 0.404 | 78 | 0.192 |
| ds_mem | 0.619 | 0.566 | 0.379 | 78 | 0.231 |
| sample_w_fused | 0.619 | 0.596 | 0.441 | 240 | 0.163 |
| subset_w_fused | 0.619 | 0.566 | 0.476 | 240 | 0.163 |
| both_fused | 0.524 | 0.535 | 0.506 | 240 | 0.163 |
| ds_known+sample_f | 1.000 | 0.414 | 0.560 | 78 | 0.500 |

**news20b** — TF-IDF 4-topic; known atheism/religion → sparse med/hardware

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.125 | 0.692 | 0.286 | 147 | 0.218 |
| ds_known | 0.125 | 0.662 | 0.304 | 64 | 0.500 |
| ds_random | 0.125 | 0.738 | 0.246 | 64 | 0.234 |
| ds_mem | 0.000 | 0.646 | 0.258 | 64 | 0.297 |
| sample_w_fused | 0.375 | 0.538 | 0.376 | 147 | 0.218 |
| subset_w_fused | 0.250 | 0.554 | 0.402 | 147 | 0.218 |
| both_fused | 0.250 | 0.508 | 0.406 | 147 | 0.218 |
| ds_known+sample_f | 0.375 | 0.585 | 0.436 | 64 | 0.500 |

**isolet** — ISOLET speech; known letters A–M → sparse novel N–Z

| scheme | novel | known | forgetting | n_train | novel_frac |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.350 | 0.750 | 0.200 | 234 | 0.214 |
| ds_known | 0.350 | 0.729 | 0.277 | 100 | 0.500 |
| ds_random | 0.300 | 0.750 | 0.231 | 100 | 0.210 |
| ds_mem | 0.300 | 0.771 | 0.187 | 100 | 0.350 |
| sample_w_fused | 0.250 | 0.729 | 0.256 | 234 | 0.214 |
| subset_w_fused | 0.250 | 0.625 | 0.321 | 234 | 0.214 |
| both_fused | 0.300 | 0.646 | 0.310 | 234 | 0.214 |
| ds_known+sample_f | 0.350 | 0.573 | 0.366 | 100 | 0.500 |

```bash
python3 scripts/run_downsample_vs_weight_ablation.py --quick
```
