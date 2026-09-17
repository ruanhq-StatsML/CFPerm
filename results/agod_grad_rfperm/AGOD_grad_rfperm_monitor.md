# Grad-OnlineRFPerm (`param.grad.norm`) — Sep17 MVP datasets

Frozen `f_ref` MLP pretrained on burn-in; OnlineRFPerm on per-layer
`||∇_θ L||` (+ relative shares) vs MSE-OnlineRFPerm / MMD² / PO.

Lead: `t_grad − t_mse` (negative ⇒ Grad earlier). seeds=[0, 1, 2, 3, 4], batch=128, n_batches=48, n_burn=8, alpha=0.05

| dataset | mean lead(g−mse) | median | P(earlier) | mean lead(share−mse) | P(share earlier) |
|---|---:|---:|---:|---:|---:|
| `synthetic` | -3.40 | -3.0 | 100% | -2.80 | 80% |
| `covertype` | -1.20 | -1.0 | 60% | -1.00 | 60% |
| `bank` | -1.20 | +0.0 | 40% | +0.00 | 40% |
| `electricity` | -6.20 | -4.0 | 80% | -6.00 | 80% |
| `eeg` | -4.60 | -6.0 | 80% | -4.00 | 80% |

**Overall** (25 runs): mean lead(g−mse)=-3.32, P(Grad earlier)=72%, P(≤0)=96%.
Relative-share: mean=-2.76, P(earlier)=68%.

## Method

```
pretrain frozen f_ref MLP on burn-in (never stepped afterward)
for each stream batch t:
  T_grad = ||∇_θ L(batch; f_ref)|| - e_ref
  p = rank/EWMA; alpha-investing FDR → reject
  compare to MSE-OnlineRFPerm / MMD / PO
```

Freeze hint: earliest rejecting layer guides back-prop depth.

