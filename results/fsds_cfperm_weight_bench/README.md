# FSDS/CFPerm weight benchmark

FSDS RF-domain VIMP → membership → soft sample_w + hard observation subset_w on amazon/digits/news20/synth_cgcd. Pure membership is diagnostic; fused (membership ∪ private-live) is the cold-start weight channel.

## Aggregate (vs uniform)

| channel | mean Δ novel | mean Δ acc |
| --- | --- | --- |
| sample_w (membership) | -0.043 | — |
| subset_w (hard obs) | -0.037 | — |
| sample × subset | 0.020 | -0.017 |
| sample_w fused | — | 0.005 |
| subset_w fused | — | 0.007 |
| both fused | -0.018 | -0.012 |

## Per dataset

| dataset | n_ref/n_live | subset_frac | vimp_top | Δacc both | Δacc both_fused | Δnovel both_fused |
| --- | --- | --- | --- | --- | --- | --- |
| digits | 400/192 | 0.401 | 0.220 | -0.078 | -0.109 | -0.071 |
| news20 | 500/216 | 0.398 | 0.275 | 0.000 | 0.083 | 0.000 |
| synth_cgcd | 240/254 | 0.402 | 0.384 | 0.012 | 0.060 | 0.000 |
| amazon | 120/150 | 0.400 | 0.428 | 0.000 | -0.080 | 0.000 |

### Scheme detail (acc / novel / forgetting)

**digits** — known 0-4 REF; LIVE sparse inject 5-9 (~1:5)

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.938 | 0.786 | 0.980 | 0.047 | 0.000 |
| sample_w | 0.875 | 0.714 | 0.920 | 0.078 | 0.235 |
| subset_w | 0.859 | 0.571 | 0.940 | 0.078 | 0.268 |
| sample_x_subset | 0.859 | 0.714 | 0.900 | 0.095 | 0.448 |
| sample_w_fused | 0.891 | 0.786 | 0.920 | 0.075 | 0.228 |
| subset_w_fused | 0.875 | 0.714 | 0.920 | 0.073 | 0.268 |
| both_fused | 0.828 | 0.714 | 0.860 | 0.108 | 0.437 |

**news20** — TF-IDF; known sci/rec → sparse inject comp.graphics (~1:5)

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.764 | 0.800 | 0.758 | 0.204 | 0.000 |
| sample_w | 0.764 | 0.700 | 0.774 | 0.220 | 0.292 |
| subset_w | 0.778 | 0.700 | 0.790 | 0.214 | 0.267 |
| sample_x_subset | 0.764 | 0.700 | 0.774 | 0.256 | 0.509 |
| sample_w_fused | 0.847 | 0.800 | 0.855 | 0.226 | 0.301 |
| subset_w_fused | 0.819 | 0.800 | 0.823 | 0.182 | 0.265 |
| both_fused | 0.847 | 0.800 | 0.855 | 0.236 | 0.490 |

**synth_cgcd** — synthetic sparse novel C-GCD

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.714 | 1.000 | 0.711 | 0.275 | 0.000 |
| sample_w | 0.786 | 1.000 | 0.783 | 0.246 | 0.247 |
| subset_w | 0.798 | 1.000 | 0.795 | 0.225 | 0.267 |
| sample_x_subset | 0.726 | 1.000 | 0.723 | 0.304 | 0.452 |
| sample_w_fused | 0.738 | 1.000 | 0.735 | 0.283 | 0.246 |
| subset_w_fused | 0.810 | 1.000 | 0.807 | 0.271 | 0.267 |
| both_fused | 0.774 | 1.000 | 0.771 | 0.308 | 0.447 |

**amazon** — category hop: REF cats{0,1}; LIVE recirculate + sparse cats{2,3}

| scheme | acc | novel | known | forgetting | weight_gini |
| --- | --- | --- | --- | --- | --- |
| uniform | 0.680 | 0.250 | 0.816 | 0.108 | 0.000 |
| sample_w | 0.640 | 0.250 | 0.763 | 0.225 | 0.280 |
| subset_w | 0.680 | 0.417 | 0.763 | 0.158 | 0.262 |
| sample_x_subset | 0.680 | 0.500 | 0.737 | 0.217 | 0.456 |
| sample_w_fused | 0.640 | 0.167 | 0.789 | 0.158 | 0.256 |
| subset_w_fused | 0.620 | 0.083 | 0.789 | 0.192 | 0.261 |
| both_fused | 0.600 | 0.250 | 0.711 | 0.300 | 0.427 |

Weights: `sample_w = 1+λ·membership`; `subset_w` boosts observations with membership ≥ q.

```bash
python3 scripts/run_fsds_cfperm_weight_benchmark.py --quick
```
