# Instrumental Variable selection — synthetic argument

Hidden confounder `U → (T, Y)`. Naive OLS / partialling is biased.
2SLS with **valid** instruments recovers `τ`; invalid / weak Z do not.

Selection rule (`iv_auto`): first-stage partial **F ≥ 10** and
exclusion residual score ≤ 0.12.

## Absolute error `|ATÊ − τ|` (lower better)

| method | abs_err | ATÊ | n_Z | valid / invalid picked |
| --- | --- | --- | --- | --- |
| naive_ols | 0.301 (0.074) | 1.801 | 0.0 | 0.0 / 0.0 |
| naive_partial | 0.301 (0.074) | 1.801 | 0.0 | 0.0 / 0.0 |
| iv_oracle | 0.017 (0.007) | 1.510 | 3.0 | 3.0 / 0.0 |
| iv_auto | 0.030 (0.030) | 1.523 | 3.1 | 2.9 / 0.2 |
| iv_strong | 0.052 (0.046) | 1.509 | 4.5 | 3.0 / 1.5 |
| iv_all | 0.058 (0.046) | 1.501 | 8.0 | 3.0 / 3.0 |
| iv_invalid | 2.290 (2.966) | 2.308 | 3.0 | 0.0 / 3.0 |
| iv_weak | 0.373 (0.291) | 1.518 | 2.0 | 0.0 / 0.0 |

## IV selection quality (`iv_auto` vs true valid set)

| precision | recall | F1 |
| --- | --- | --- |
| 0.938 | 0.958 | 0.939 |

## Confounding sweep

| α | naive | iv_oracle | iv_auto | iv_invalid | sel F1 |
| --- | --- | --- | --- | --- | --- |
| 0.0 | 0.064 | 0.008 | 0.021 | 3.393 | 0.971 |
| 0.5 | 0.097 | 0.013 | 0.027 | 3.377 | 0.971 |
| 1.0 | 0.293 | 0.021 | 0.035 | 3.328 | 0.931 |
| 1.5 | 0.488 | 0.029 | 0.077 | 3.254 | 0.825 |
| 2.0 | 0.629 | 0.037 | 0.082 | 3.164 | 0.798 |

## Run

```bash
python3 scripts/run_iv_selection_synth.py
python3 scripts/run_iv_selection_synth.py --quick
```

Code: `Python/src/iv_selection_synth.py`
