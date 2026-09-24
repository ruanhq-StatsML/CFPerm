# Soft IPTW burn decisions

burn soft IPTW on 1/6 packs; decisions={'BURN_SQRT': 1, 'SOFTEN_ONLY': 5}. Default KEEP_UNIFORM / SOFTEN_ONLY — α is not a FLOPs knob.

1) Do not burn soft weights to save FLOPs (they don't). 2) Burn only if gated_α beats uniform on sig-MSE. 3) If gate is mandatory and hurts, SOFTEN_ONLY with ∛. 4) Keep the real bill on duty × refit (算力账).

| dataset | duty | √ rel | ∛ rel | decision | burn? | α | reason |
|---|---:|---:|---:|---|:---:|---:|---|
| `metro_interstate` | 0.179 | 0.927 | 0.934 | `BURN_SQRT` | Y | 0.500 | both beat uniform; √ has lower rel MSE |
| `beijing_pm25` | 0.128 | 1.012 | 1.068 | `SOFTEN_ONLY` | N | 0.500 | gated IPTW does not beat uniform; √ hurts less than ∛ here |
| `stocks_AAPL` | 0.179 | 1.145 | 1.103 | `SOFTEN_ONLY` | N | 0.333 | gated IPTW does not beat uniform; if reject path is mandatory, use ∛ to limit damage (same FLOPs as √) |
| `waymo_proxy` | 0.564 | 1.083 | 1.050 | `SOFTEN_ONLY` | N | 0.333 | gated IPTW does not beat uniform; if reject path is mandatory, use ∛ to limit damage (same FLOPs as √) |
| `stocks_MSFT` | 0.231 | 1.115 | 1.026 | `SOFTEN_ONLY` | N | 0.333 | gated IPTW does not beat uniform; if reject path is mandatory, use ∛ to limit damage (same FLOPs as √) |
| `stocks_IWM` | 0.179 | 1.090 | 1.039 | `SOFTEN_ONLY` | N | 0.333 | gated IPTW does not beat uniform; if reject path is mandatory, use ∛ to limit damage (same FLOPs as √) |

## 算力账

- Adaptation FLOPs = duty × refit（真账单）
- Weighting FLOPs ≈ O(n)（α 不改总账）
- 因此：该不该烧看 MSE 风险，不看「省不算力」

