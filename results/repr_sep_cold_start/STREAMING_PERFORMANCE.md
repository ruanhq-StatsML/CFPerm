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
