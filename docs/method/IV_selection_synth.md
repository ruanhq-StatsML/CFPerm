# Instrumental Variable selection — synthetic argument

Hidden confounder `U → (T, Y)`. Naive OLS / partialling is biased.
2SLS with **valid** instruments recovers `τ`; invalid / weak Z do not.

Selection rule (`iv_auto`): first-stage partial **F ≥ 10** and
exclusion excess screen.

## Absolute error `|ATÊ − τ|` (lower better)

| method | abs_err | ATÊ | n_Z | valid / invalid picked |
| --- | --- | --- | --- | --- |
| naive_ols | 0.221 (0.050) | 1.721 | 0.0 | 0.0 / 0.0 |
| naive_partial | 0.221 (0.050) | 1.721 | 0.0 | 0.0 / 0.0 |
| iv_oracle | 0.020 (0.008) | 1.512 | 3.0 | 3.0 / 0.0 |
| iv_auto | 0.020 (0.008) | 1.512 | 3.0 | 3.0 / 0.0 |
| iv_strong | 0.053 (0.065) | 1.494 | 4.2 | 3.0 / 1.2 |
| iv_all | 0.073 (0.059) | 1.486 | 8.0 | 3.0 / 3.0 |
| iv_invalid | 2.597 (4.065) | 2.497 | 3.0 | 0.0 / 3.0 |
| iv_weak | 0.439 (0.491) | 1.730 | 2.0 | 0.0 / 0.0 |

## IV selection quality (`iv_auto` vs true valid set)

| precision | recall | F1 |
| --- | --- | --- |
| 1.000 | 1.000 | 1.000 |

## Confounding sweep

| α | naive | iv_oracle | iv_auto | iv_invalid | sel F1 |
| --- | --- | --- | --- | --- | --- |
| 0.0 | 0.043 | 0.009 | 0.026 | 4.217 | 0.971 |
| 0.5 | 0.060 | 0.015 | 0.037 | 4.028 | 0.931 |
| 1.0 | 0.218 | 0.024 | 0.024 | 3.839 | 1.000 |
| 1.5 | 0.390 | 0.033 | 0.051 | 3.654 | 0.971 |
| 2.0 | 0.532 | 0.042 | 0.082 | 3.476 | 0.903 |

## Sensitivity on the naive association

| metric | mean (sd) |
| --- | --- |
| E-value (point) | 5.399 (0.524) |
| E-value (CI) | 5.287 (0.510) |
| approx RR | 2.975 (0.259) |
| Rosenbaum Γ* (p-bound > 0.05) | 16.88 (3.72) |
| Γ* → ρ amplification | 0.883 (0.030) |
| OVB robustness value RV | 0.879 (0.013) |
| OVB RV (halve effect) | 0.800 (0.018) |
| OVB × vs max observed partial R² | 8.5 (6.6) |
| Conley γ* (IV crosses 0) | 3.218 (0.634) |
| naive ATÊ | 1.721 (0.050) |
| iv ATÊ | 1.512 (0.018) |

E-value: RR-scale confounder strength to kill the naive effect.
Rosenbaum Γ*: matched high/low-T sign-test odds bias; ρ=(Γ−1)/(Γ+1).
OVB RV (Cinelli–Hazlett): partial-R² both sides needed to nullify OLS.
Conley γ*: how large a direct Z→Y path can be before 2SLS crosses 0.

## Subset localization — IV strength heterogeneity

Segments = quartiles of `X0` (designed stronger first stage in Q4).

| segment | n | stage-1 F | R² | ATÊ | abs_err | scale |
| --- | --- | --- | --- | --- | --- | --- |
| Q1_low | 500 | 1.2 (1.1) | 0.007 | 2.142 | 0.642 | 0.03 |
| Q2 | 500 | 13.6 (7.1) | 0.076 | 1.563 | 0.149 | 0.22 |
| Q3 | 500 | 212.9 (35.6) | 0.563 | 1.508 | 0.041 | 1.01 |
| Q4_high | 500 | 1111.9 (213.1) | 0.869 | 1.510 | 0.010 | 2.27 |

**F heterogeneity ratio** (max/min across segments): **929.89×**.

## After localization — strategy adjustment

Map local first-stage F → tiers, then change aggregation / where you act.

| tier | share | rule |
| --- | --- | --- |
| NO_ACT | 31.2% | F < 10 → hold; do not quote IV ATE |
| ESTIMATE_ONLY | 15.6% | 10 ≤ F < 20 → publish + WATCH |
| ACT | 53.1% | F ≥ 20 → eligible for F-/precision weights |

ID coverage=69% · ACT coverage=53% · policy conflict=31%.

| aggregator | |ATÊ − τ| |
| --- | --- |
| global_iv | 0.020 (0.008) |
| naive_ols | 0.221 (0.050) |
| equal_pool_all | 0.181 (0.120) |
| equal_pool_id | 0.044 (0.055) |
| f_weighted | 0.012 (0.007) |
| precision_wt | 0.012 (0.007) |
| adaptive_act | 0.011 (0.007) |

Insights:

- best_aggregator=adaptive_act (|err|=0.011) vs global_iv=0.020 — weight by local F/precision, don't flat-pool NO_ACT
- tier_mix: NO_ACT=31% ESTIMATE_ONLY=16% ACT=53%
- coverage mean: identify=69% act=53%; policy_conflict=31% of mass would get a global IV quote while local F<10
- equal_pool_all |err|=0.181 vs f_weighted=0.012 — including weak segments pollutes the pooled ATE
- strength_gap: Q1_low F=3.4 (NO_ACT) vs Q4_high F=1246.8 (ACT)
- coverage: identify=75%  act=50%  hold_out=25%

## Run

```bash
python3 scripts/run_iv_selection_synth.py
python3 scripts/run_iv_selection_synth.py --quick
```

Code: `Python/src/iv_selection_synth.py`
