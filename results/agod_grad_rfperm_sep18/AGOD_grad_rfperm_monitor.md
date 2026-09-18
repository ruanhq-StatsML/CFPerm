# Grad-OnlineRFPerm (`param.grad.norm`) — Sep17 MVP datasets

**Single stream:** `g_t = ||∇_{θ_U} L||_2` over all unfrozen params,
one OnlineRFPerm (no per-layer / per-param multiple testing).
Layer shares = diagnostics only after global reject.

Lead: `t_grad − t_mse` (negative ⇒ Grad earlier). seeds=[0, 1, 2, 3, 4], batch=128, n_batches=48, n_burn=8, alpha=0.05

| dataset | mean lead(g−mse) | median | P(earlier) | P(≤0) |
|---|---:|---:|---:|---:|
| `adult` | -3.20 | -3.0 | 80% | 100% |
| `nomao` | -4.00 | -4.0 | 80% | 100% |

**Overall** (10 runs): mean lead(g−mse)=-3.60, P(Grad earlier)=80%, P(≤0)=100%.

## Method (engineering口径)

```
θ_U = unfrozen params (requires_grad=True)
g_t = ||∇_{θ_U} L||_2 = sqrt(Σ_i ||∇θ_i||²)   # ONE scalar — not mean of norms
OnlineRFPerm once on T_t = g_t - e_ref
shares_ℓ = ||g_ℓ|| / g_t   # diagnostic ranking only, no FDR per layer
```

