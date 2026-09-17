# Grad-OnlineRFPerm — supplementary experiments

Single-stream `g_t=||∇_{θ_U} L||_2` (no per-layer multiple testing).

## A) Null / grace FPR (stationary synthetic)

| grace | P(any reject) | P(early≤5) | mean t_grad | mean Grad duty | mean MSE duty |
|---:|---:|---:|---:|---:|---:|
| 0 | 100% | 100% | 8.0 | 0.233 | 0.225 |
| 2 | 100% | 100% | 10.7 | 0.211 | 0.228 |
| 4 | 100% | 67% | 15.3 | 0.204 | 0.222 |

## B) Alpha sensitivity (mean lead g−mse)

| dataset | α=0.01 | α=0.05 | α=0.10 |
|---|---:|---:|---:|
| `synthetic` | -6.33 | -3.67 | -3.67 |
| `electricity` | -7.00 | -7.00 | -2.33 |

## C) Freeze closed-loop (post-reject MSE / always_adapt, FLOPs proxy)

### `synthetic`

| policy | mse_post | FLOPs proxy |
|---|---:|---:|
| `always_adapt` | 0.9529 (1.00×) | 1024320 (1.00×) |
| `freeze_early` | 1.6714 (1.75×) | 679061 (0.66×) |
| `freeze_low_share` | 1.0980 (1.15×) | 582613 (0.57×) |
| `no_adapt` | 13.5582 (14.23×) | 8536 (0.01×) |

### `electricity`

| policy | mse_post | FLOPs proxy |
|---|---:|---:|
| `always_adapt` | 0.3403 (1.00×) | 860480 (1.00×) |
| `freeze_early` | 0.2922 (0.86×) | 676160 (0.79×) |
| `freeze_low_share` | 0.2955 (0.87×) | 573653 (0.67×) |
| `no_adapt` | 1.7335 (5.09×) | 0 (0.00×) |

## D) Extra stream packs (mean lead g−mse)

| pack | mean lead | P(earlier) | P(≤0) |
|---|---:|---:|---:|
| `metro_interstate` | +2.00 | 0% | 0% |
| `beijing_pm25` | +0.00 | 33% | 67% |
| `stocks_MSFT` | -1.67 | 67% | 67% |
| `stocks_IWM` | -4.00 | 100% | 100% |
| `waymo_proxy` | +0.00 | 0% | 100% |

## Takeaways

- **Null:** over a long horizon alpha-investing still rejects eventually; **grace** delays first reject (mean t: 8→15) and cuts early(≤5) FPR (100%→67% at grace=4).
- **Alpha:** lead(g−mse) stays negative for α∈{0.01,0.05,0.10} on synthetic / electricity.
- **Freeze loop:** on electricity, `freeze_early` / `freeze_low_share` beat `always_adapt` on post-reject MSE (~0.86×) at ~0.7× FLOPs; on synthetic, `freeze_low_share` ≈1.15× MSE at 0.57× FLOPs. `no_adapt` collapses.
- **Extra packs:** stocks_IWM lead −4.0 (100% earlier); stocks_MSFT −1.7; waymo / beijing ≈0; metro Grad later (+2).

