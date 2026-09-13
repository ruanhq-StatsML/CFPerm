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
| `affec` | `mean_ce` | -0.0501 | -0.055 | 0.560 | 0.000 | eye_tracking:0.20 / pupil:0.20 / cursor:0.20 / gsr_eda:0.20 / eeg:0.20 |
| `affec` | `stack_ce` | +0.0180 | +0.019 | 0.629 | 0.084 | eye_tracking:0.21 / pupil:0.21 / cursor:0.21 / gsr_eda:0.15 / eeg:0.21 |
| `affec` | `stack_alpha` | +0.0369 | +0.039 | 0.622 | 0.084 | eye_tracking:0.21 / pupil:0.21 / cursor:0.21 / gsr_eda:0.15 / eeg:0.21 |

Best (MSE drop, Acc lift): **`amazon/stack_alpha`**
(MSE drop=+0.0383, Acc lift=+0.069).

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py
PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py --datasets coco affec
```
