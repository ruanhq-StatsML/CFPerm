# Grad-OnlineRFPerm (`param.grad.norm`) — Sep17 MVP datasets

**Single stream:** `g_t = ||∇_{θ_U} L||_2` over all unfrozen params,
one OnlineRFPerm (no per-layer / per-param multiple testing).
Layer shares = diagnostics only after global reject.

Lead: `t_grad − t_mse` (negative ⇒ Grad earlier). seeds=[0, 1, 2, 3, 4], batch=128, n_batches=48, n_burn=8, alpha=0.05

| dataset | mean lead(g−mse) | median | P(earlier) | P(≤0) |
|---|---:|---:|---:|---:|
| `synthetic` | -3.00 | -2.0 | 100% | 100% |
| `covertype` | -1.20 | -1.0 | 60% | 80% |
| `bank` | -1.20 | +0.0 | 40% | 100% |
| `electricity` | -6.20 | -4.0 | 80% | 100% |
| `eeg` | -4.60 | -6.0 | 80% | 100% |

**Overall** (25 runs): mean lead(g−mse)=-3.24, P(Grad earlier)=72%, P(≤0)=96%.

## Method (engineering口径)

```
θ_U = unfrozen params (requires_grad=True)
g_t = ||∇_{θ_U} L||_2 = sqrt(Σ_i ||∇θ_i||²)   # ONE scalar — not mean of norms
OnlineRFPerm once on T_t = g_t - e_ref
shares_ℓ = ||g_ℓ|| / g_t   # diagnostic ranking only, no FDR per layer
```

