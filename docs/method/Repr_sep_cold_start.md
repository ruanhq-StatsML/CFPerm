# Streaming C-GCD proxy — performance

Protocol: C-GCD streaming proxy (quick=False). Seeds [2026, 2027, 2028].
KPI: **novel_acc** primary; known retention / forgetting; efficiency n*.

| mode | known | novel | forgetting | Δnovel vs uniform | n* (mean) |
| --- | --- | --- | --- | --- | --- |
| repr_sep | 0.743±0.014 | 0.518±0.055 | 0.243±0.009 | +0.114 | 22.7 |
| weighted | 0.771±0.004 | 0.471±0.021 | 0.214±0.022 | +0.067 | 12.0 |
| membership | 0.799±0.012 | 0.273±0.056 | 0.188±0.012 | -0.131 | 24.0 |
| uniform | 0.783±0.028 | 0.404±0.098 | 0.224±0.019 | +0.000 | 25.3 |

## Per-stage novel_acc (mean over seeds)

| mode | s0 | s1 | s2 | s3 |
| --- | --- | --- | --- | --- |
| repr_sep | 0.700 | 0.532 | 0.477 | 0.363 |
| weighted | 0.392 | 0.539 | 0.425 | 0.529 |
| membership | 0.071 | 0.370 | 0.260 | 0.391 |
| uniform | 0.467 | 0.390 | 0.338 | 0.422 |

# Streaming cold-start performance

Streaming KPI = novel_acc (+ known/forgetting). repr_sep = opportunity-gated attribute → sep on S → expand.

## Summary (this seed)

| mode | known | novel | forgetting | feat_hit |
| --- | --- | --- | --- | --- |
| repr_sep | 0.760 | 0.457 | 0.246 | 0.250 |
| weighted | 0.771 | 0.495 | 0.244 | 0.250 |
| uniform | 0.804 | 0.516 | 0.248 | — |

Δ novel (repr − uniform): **-0.059**
Δ novel (repr − weighted): **-0.039**

## Streaming per-stage (repr_sep)

| stage | novel | known | forget | phase | O_local | gap_lift |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.444 | 0.758 | 0.269 | expand | 0.720 | 0.395 |
| 1 | 0.571 | 0.724 | 0.229 | expand | 0.745 | 0.314 |
| 2 | 0.455 | 0.723 | 0.287 | expand | 0.741 | 0.313 |
| 3 | 0.357 | 0.833 | 0.200 | expand | 0.757 | 0.325 |
