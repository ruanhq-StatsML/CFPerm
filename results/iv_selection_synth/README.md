# Instrumental Variable selection — synthetic argument

Hidden confounder `U → (T, Y)`. Naive OLS / partialling is biased.
2SLS with **valid** instruments recovers `τ`; invalid / weak Z do not.

Selection rule (`iv_auto`): first-stage partial **F ≥ 10** and
exclusion excess screen.

## Absolute error `|ATÊ − τ|` (lower better)

| method | abs_err | ATÊ | n_Z | valid / invalid picked |
| --- | --- | --- | --- | --- |
| naive_ols | 0.157 (0.041) | 1.657 | 0.0 | 0.0 / 0.0 |
| naive_partial | 0.157 (0.041) | 1.657 | 0.0 | 0.0 / 0.0 |
| iv_oracle | 0.011 (0.004) | 1.506 | 3.0 | 3.0 / 0.0 |
| iv_auto | 0.011 (0.004) | 1.506 | 3.0 | 3.0 / 0.0 |
| iv_strong | 0.022 (0.020) | 1.498 | 4.1 | 3.0 / 1.1 |
| iv_all | 0.025 (0.020) | 1.501 | 8.0 | 3.0 / 3.0 |
| iv_invalid | 2.553 (2.977) | 2.067 | 3.0 | 0.0 / 3.0 |
| iv_weak | 0.452 (0.384) | 1.461 | 2.0 | 0.0 / 0.0 |

## IV selection quality (`iv_auto` vs true valid set)

| precision | recall | F1 |
| --- | --- | --- |
| 1.000 | 1.000 | 1.000 |

## Confounding sweep

| α | naive | iv_oracle | iv_auto | iv_invalid | sel F1 |
| --- | --- | --- | --- | --- | --- |
| 0.0 | 0.027 | 0.005 | 0.010 | 3.412 | 0.971 |
| 0.5 | 0.045 | 0.008 | 0.011 | 3.412 | 0.960 |
| 1.0 | 0.155 | 0.013 | 0.013 | 3.371 | 1.000 |
| 1.5 | 0.296 | 0.018 | 0.033 | 3.301 | 0.863 |
| 2.0 | 0.428 | 0.023 | 0.065 | 3.211 | 0.745 |

## Sensitivity on the naive association

| metric | mean (sd) |
| --- | --- |
| E-value (point) | 5.067 (0.442) |
| E-value (CI) | 4.977 (0.431) |
| approx RR | 2.811 (0.218) |
| Rosenbaum Γ* (p-bound > 0.05) | 20.00 (6.55) |
| Γ* → ρ amplification | 0.898 (0.027) |
| OVB robustness value RV | 0.900 (0.011) |
| OVB RV (halve effect) | 0.831 (0.016) |
| OVB × vs max observed partial R² | 8.6 (6.4) |
| Conley γ* (IV crosses 0) | 6.663 (0.380) |
| naive ATÊ | 1.657 (0.041) |
| iv ATÊ | 1.506 (0.010) |

E-value: RR-scale confounder strength to kill the naive effect.
Rosenbaum Γ*: matched high/low-T sign-test odds bias; ρ=(Γ−1)/(Γ+1).
OVB RV (Cinelli–Hazlett): partial-R² both sides needed to nullify OLS.
Conley γ*: how large a direct Z→Y path can be before 2SLS crosses 0.

## Subset localization — IV strength heterogeneity

Segments = quartiles of `X0` (designed stronger first stage in Q4).

| segment | n | stage-1 F | R² | ATÊ | abs_err | scale |
| --- | --- | --- | --- | --- | --- | --- |
| Q1_low | 500 | 374.9 (80.1) | 0.691 | 1.501 | 0.041 | 1.27 |
| Q2 | 500 | 551.7 (144.9) | 0.764 | 1.504 | 0.019 | 1.49 |
| Q3 | 500 | 695.4 (119.0) | 0.807 | 1.505 | 0.024 | 1.69 |
| Q4_high | 500 | 940.4 (183.6) | 0.848 | 1.512 | 0.012 | 1.92 |

**F heterogeneity ratio** (max/min across segments): **2.51×**.

## Run

```bash
python3 scripts/run_iv_selection_synth.py
python3 scripts/run_iv_selection_synth.py --quick
```

Code: `Python/src/iv_selection_synth.py`
