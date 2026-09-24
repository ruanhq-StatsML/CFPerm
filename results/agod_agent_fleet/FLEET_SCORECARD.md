# Agent fleet CL scorecard (N≈180)

**Verdict:** even-force lower headcount HHI (0.139<0.430); even-force lower WIP Gini; even-force mean commercial score ≥ pile-on; antifraud useful lift>0 under specialty-10; CUPED var ratio=0.61<1

Observation: `per-cohort feedback features (fraud useful, CUPED pre-period X, WIP)` → Y=`Score`; intermediate=`fleet Thought = {cohort mix, gini, rehome action} — not Ŷ`

## Even force vs pile-on

| mode | n | cohorts | Gini WIP | HHI WIP | HHI headcount | mean Score | fraud lift | CUPED ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `even_force` | 180 | 10 | 0.243 | 0.129 | 0.139 | 0.778 | 0.175 | 0.613 |
| `pile_on` | 180 | 10 | 0.634 | 0.488 | 0.43 | 0.767 | -0.025 | 0.613 |

## Default cohort mix

| cohort | n |
|---|---:|
| `antifraud` | 10 |
| `cuped` | 10 |
| `review_accel` | 20 |
| `localize_fsds` | 20 |
| `board_stream` | 15 |
| `eta_gate` | 10 |
| `gray_rollback` | 15 |
| `data_contract` | 15 |
| `roi_ops` | 15 |
| `flex_reserve` | 50 |

## Reflection (even tick)

- action: `keep_specialty_mix`
- reason: fraud lift≥0, CUPED ratio≤1, WIP gini ok — average force holds

## Multi-tick CL reflection

- actions: `['rehome_flex_to_underserved', 'rehome_flex_to_underserved', 'keep_specialty_mix', 'keep_specialty_mix', 'keep_specialty_mix']`
- fraud lifts: `[-0.1, -0.2, 0.317, 0.217, 0.267]`
- CUPED ratios: `[0.966, 0.998, 0.657, 0.756, 0.711]`
- reading: early weak feedback → rehome_flex_to_underserved; later fraud lift↑ / CUPED ratio↓ → keep_specialty_mix

## Reading

- Soft-cap specialty WIP; Score only counts wired/used/roi.
- Antifraud×10: optimize human useful-rate, not AUC.
- CUPED×10: optimize Var ratio < 1, not another ranker.
- Rehome only from `flex_reserve` → underserved cohorts.

