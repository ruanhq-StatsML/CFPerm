# Attribution-guided online stacking

## Method

Online stacking learns fusion weights `w = softmax(psi)` over modality logits.

Weight socket: attribution alpha (MSG-B3 / B5 hybrid / ImgTxt hybrid / Affec VIMP, EMA) plugs in as

```
L = CE(stack_w · logits, y) + lambda * KL(stack_w || alpha)
```

Classical online stacking + an FSDS/MSG attribution prior — not a new fusion architecture.

## Variants

| variant | fusion | alpha used? |
|---|---|---|
| `mean_ce` | mean-pool | no |
| `stack_ce` | learned stack_w | logged only |
| `stack_alpha` | stack_w + KL(w||alpha) | yes |

## Smoke

| dataset | variant | MSE drop | Acc lift | Acc post | |w-alpha| | mean stack_w |
|---|---|---:|---:|---:|---:|---|
| `msrvtt` | `mean_ce` | -0.0080 | +0.000 | 0.651 | 0.000 | video:0.33 / text:0.33 / audio:0.33 |
| `msrvtt` | `stack_ce` | +0.0076 | +0.008 | 0.655 | 0.124 | video:0.37 / text:0.29 / audio:0.34 |
| `msrvtt` | `stack_alpha` | +0.0097 | +0.008 | 0.655 | 0.096 | video:0.39 / text:0.33 / audio:0.27 |
| `amazon` | `mean_ce` | +0.0339 | +0.049 | 0.543 | 0.000 | text:0.50 / image:0.50 |
| `amazon` | `stack_ce` | +0.0357 | +0.033 | 0.539 | 0.224 | text:0.47 / image:0.53 |
| `amazon` | `stack_alpha` | +0.0383 | +0.069 | 0.571 | 0.155 | text:0.39 / image:0.61 |
| `coco` | `mean_ce` | +0.0236 | +0.034 | 0.860 | 0.000 | image:0.50 / text:0.50 |
| `coco` | `stack_ce` | +0.0251 | +0.033 | 0.865 | 0.135 | image:0.39 / text:0.61 |
| `coco` | `stack_alpha` | +0.0240 | +0.036 | 0.860 | 0.060 | image:0.47 / text:0.53 |
| `coco` | `stack_uniform` | +0.0247 | +0.036 | 0.862 | 0.067 | image:0.46 / text:0.54 |
| `coco` | `stack_fixed` | +0.0213 | +0.025 | 0.857 | 0.044 | image:0.50 / text:0.50 |
| `coco` | `stack_temp` | +0.0251 | +0.037 | 0.872 | 0.200 | image:0.33 / text:0.67 |
| `coco` | `stack_erank` | +0.0246 | +0.037 | 0.854 | 0.100 | image:0.43 / text:0.57 |
| `affec` | `mean_ce` | -0.0501 | -0.055 | 0.560 | 0.000 | eye_tracking:0.20 / pupil:0.20 / cursor:0.20 / gsr_eda:0.20 / eeg:0.20 |
| `affec` | `stack_ce` | +0.0180 | +0.019 | 0.629 | 0.084 | eye_tracking:0.21 / pupil:0.21 / cursor:0.21 / gsr_eda:0.15 / eeg:0.21 |
| `affec` | `stack_alpha` | +0.0369 | +0.039 | 0.622 | 0.084 | eye_tracking:0.21 / pupil:0.21 / cursor:0.21 / gsr_eda:0.15 / eeg:0.21 |
| `affec` | `stack_uniform` | +0.0225 | +0.024 | 0.652 | 0.084 | eye_tracking:0.21 / pupil:0.21 / cursor:0.21 / gsr_eda:0.15 / eeg:0.21 |
| `affec` | `stack_fixed` | +0.0093 | +0.009 | 0.634 | 0.078 | eye_tracking:0.20 / pupil:0.20 / cursor:0.20 / gsr_eda:0.20 / eeg:0.20 |
| `affec` | `stack_temp` | +0.0742 | +0.071 | 0.625 | 0.088 | eye_tracking:0.22 / pupil:0.21 / cursor:0.23 / gsr_eda:0.12 / eeg:0.23 |
| `affec` | `stack_erank` | +0.0379 | +0.037 | 0.626 | 0.084 | eye_tracking:0.21 / pupil:0.21 / cursor:0.21 / gsr_eda:0.15 / eeg:0.21 |
| `fashion_iq` | `mean_ce` | +0.0050 | +0.012 | 0.888 | 0.000 | image:0.50 / text:0.50 |
| `fashion_iq` | `stack_ce` | +0.0025 | +0.010 | 0.888 | 0.421 | image:0.37 / text:0.63 |
| `fashion_iq` | `stack_alpha` | +0.0023 | +0.012 | 0.887 | 0.247 | image:0.57 / text:0.43 |
| `fashion_iq` | `stack_temp` | +0.0020 | +0.010 | 0.891 | 0.508 | image:0.29 / text:0.71 |
| `food101` | `mean_ce` | +0.0243 | +0.040 | 0.958 | 0.000 | image:0.50 / text:0.50 |
| `food101` | `stack_ce` | +0.0289 | +0.037 | 0.975 | 0.438 | image:0.33 / text:0.67 |
| `food101` | `stack_alpha` | +0.0239 | +0.036 | 0.961 | 0.270 | image:0.54 / text:0.46 |
| `food101` | `stack_temp` | +0.0295 | +0.037 | 0.981 | 0.542 | image:0.22 / text:0.78 |

Best (MSE drop, Acc lift): **`affec/stack_temp`**
(MSE drop=+0.0742, Acc lift=+0.071).

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py
PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py --datasets coco affec
```

## Food-101 (HF download) + Fashion-IQ

Downloaded `ethz/food101` validation via `hf download`, CLIP-packed
(`scripts/pack_food101_clip.py`). Label socket for Food-101: even/odd class id
(mode-vs-rest is too rare with 101 classes). Fashion-IQ uses existing ImgTxt pack.

| dataset | variant | MSE drop | Acc lift | Acc post |
|---|---|---:|---:|---:|
| food101 | `mean_ce` | +0.0243 | +0.040 | 0.958 |
| food101 | `stack_ce` | +0.0289 | +0.037 | 0.975 |
| food101 | `stack_alpha` | +0.0239 | +0.036 | 0.961 |
| food101 | `stack_temp` | +0.0295 | +0.037 | 0.981 |
| fashion_iq | `mean_ce` | +0.0050 | +0.012 | 0.888 |
| fashion_iq | `stack_ce` | +0.0025 | +0.010 | 0.888 |
| fashion_iq | `stack_alpha` | +0.0023 | +0.012 | 0.887 |
| fashion_iq | `stack_temp` | +0.0020 | +0.010 | 0.891 |

Reading: both streams are already strong under mean fusion (Food-101 Acc post $\approx 0.96$).
Stacking still helps MSE on Food-101 (`stack_temp` / `stack_ce`); Fashion-IQ gains are small —
same pattern as COCO (saturated baseline).
