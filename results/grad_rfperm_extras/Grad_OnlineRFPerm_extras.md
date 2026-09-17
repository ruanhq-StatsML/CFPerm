# Grad-OnlineRFPerm — supplementary experiments

Single-stream `g_t=||∇_{θ_U} L||_2` (no per-layer multiple testing).

## A) Baseline alarm rate $= n_{\mathrm{alarm}} / n_{\mathrm{batch}}$

Definition: **alarm rate = (# reject batches) / (# observation batches)** on the post-burn (optionally post-grace) window. **Not** a Type-I FAR (cannot claim stationary DGP under ongoing updates) — baseline reference only.

| grace | Grad AR | MSE AR | mean n_alarm/n_batch | mean t_grad |
|---:|---:|---:|---:|---:|
| 0 | 0.233 | 0.225 | 9.3/40 | 8.0 |
| 2 | 0.211 | 0.228 | 8.0/38 | 10.7 |
| 4 | 0.204 | 0.222 | 7.3/36 | 15.3 |

## Takeaways

- **Baseline alarm rate:** $\mathrm{AR}=n_{\mathrm{alarm}}/n_{\mathrm{batch}}$ (alarms over observation batches). Grad FAR ≈ MSE FAR (~0.20–0.23); grace mainly delays first reject, mild FAR drop.
- **Alpha:** lead(g−mse) stays negative for α∈{0.01,0.05,0.10} on synthetic / electricity (directionally stable).
- **Freeze loop:** on electricity, `freeze_early` / `freeze_low_share` beat `always_adapt` on post-reject MSE (~0.86×) at ~0.7× FLOPs; on synthetic, `freeze_low_share` ≈1.15× MSE at 0.57× FLOPs (better than freeze_early). `no_adapt` collapses.
- **Extra packs:** stocks_IWM lead −4.0 (100% earlier); stocks_MSFT −1.7; waymo / beijing ≈0; metro Grad later (+2) — pack-dependent, not universal.

