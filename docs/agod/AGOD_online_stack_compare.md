# Attribution-guided online stacking (Amazon + MSR-VTT)

## Method

Online stacking learns fusion weights `w = softmax(psi)` over modality logits.

Weight socket: attribution alpha (Amazon MSG-B3 / MSR-VTT B5 hybrid, EMA) plugs in as

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

Best (MSE drop, Acc lift): **`amazon/stack_alpha`**
(MSE drop=+0.0383, Acc lift=+0.069).

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stack_compare.py
```
