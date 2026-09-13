# AGOD PO-risk → next-step training metric compare

Iterate metric versions that map **PO-risk sensors → α → next-window training actuators**
(LR multipliers, step budget, freeze mask, stack prior).

Causal loop: sensors/α at window `t` only set actuators for window `t+1`
(no same-window leakage).

## Metric versions

| version | idea | next-step train role |
|---|---|---|
| `equal` | uniform α | dense baseline (all mods equal LR/steps) |
| `po_soft` | Softmax(PO / τ) | steer LR/steps toward high residual-concept mods |
| `po_minus_cov` | Softmax((PO − λ·MMD) / τ) | discount covariate mush before spending steps |
| `po_gated` | drift-vs-noise gate → Softmax | freeze/skip noisy windows; FLOPs saver |
| `po_proto` | Softmax(PO·(1+proto) − λ·MMD) | amplify true centroid moves |
| `po_delta` | Softmax(EMA(PO) + γ·ΔPO) | anticipatory reallocation before Acc drops |
| `po_budget` | floor + Softmax(PO) | keep all mods warm; soft reweight only |
| `po_next` | α-hist ⊕ ΔPO forecast | next-α forecast for stack prior + LR |

## Actuators (next window)

- `lr_mult[m]` from simplex α (`alpha_to_lr`, β controls floor)
- `step_alloc` proportional to α (sum ≈ `total_steps`)
- `freeze_mask` when α_m < `freeze_theta` (BWD off; FWD still on)
- `stack_prior = α` for KL(stack_w ‖ α)
- efficiency proxy: `flops_rel = (#active mods) / M`

## affec

mods = `['eye_tracking', 'pupil', 'cursor', 'gsr_eda', 'eeg']`

| version | Acc↑ | MSE↓ | Acc lift vs equal | MSE drop vs equal | FLOPs_rel | H(α) |
|---|---:|---:|---:|---:|---:|---:|
| `equal` | 0.5078 | 0.4922 | +0.0000 | +0.0000 | 1.000 | 1.000 |
| `po_soft` | 0.5059 | 0.4941 | -0.0020 | -0.0020 | 0.825 | 0.974 |
| `po_minus_cov` | 0.5039 | 0.4961 | -0.0039 | -0.0039 | 0.825 | 0.975 |
| `po_gated` | 0.5039 | 0.4961 | -0.0039 | -0.0039 | 0.700 | 0.970 |
| `po_proto` | 0.5098 | 0.4902 | +0.0020 | +0.0020 | 0.825 | 0.961 |
| `po_delta` | 0.5039 | 0.4961 | -0.0039 | -0.0039 | 0.775 | 0.970 |
| `po_budget` | 0.5059 | 0.4941 | -0.0020 | -0.0020 | 1.000 | 0.998 |
| `po_next` | 0.5059 | 0.4941 | -0.0020 | -0.0020 | 0.925 | 0.991 |

### Opportunity ranking (vs `equal`)

1. `po_gated` opp_score=0.0173 (acc_lift=-0.0039, mse_drop=-0.0039, flops=0.700)
2. `po_proto` opp_score=0.0165 (acc_lift=+0.0020, mse_drop=+0.0020, flops=0.825)
3. `po_delta` opp_score=0.0060 (acc_lift=-0.0039, mse_drop=-0.0039, flops=0.775)
4. `po_soft` opp_score=0.0042 (acc_lift=-0.0020, mse_drop=-0.0020, flops=0.825)
5. `po_minus_cov` opp_score=-0.0017 (acc_lift=-0.0039, mse_drop=-0.0039, flops=0.825)

## food101

mods = `['img', 'txt']`

| version | Acc↑ | MSE↓ | Acc lift vs equal | MSE drop vs equal | FLOPs_rel | H(α) |
|---|---:|---:|---:|---:|---:|---:|
| `equal` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 1.000 |
| `po_soft` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 0.966 |
| `po_minus_cov` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 0.967 |
| `po_gated` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 0.966 |
| `po_proto` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 0.966 |
| `po_delta` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 0.977 |
| `po_budget` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 0.984 |
| `po_next` | 0.5176 | 0.4824 | +0.0000 | +0.0000 | 1.000 | 0.998 |

### Opportunity ranking (vs `equal`)

1. `po_proto` opp_score=-0.0158 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
2. `po_soft` opp_score=-0.0158 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
3. `po_gated` opp_score=-0.0158 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
4. `po_minus_cov` opp_score=-0.0158 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
5. `po_delta` opp_score=-0.0164 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)

## fashion_iq

mods = `['img', 'txt']`

| version | Acc↑ | MSE↓ | Acc lift vs equal | MSE drop vs equal | FLOPs_rel | H(α) |
|---|---:|---:|---:|---:|---:|---:|
| `equal` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 1.000 |
| `po_soft` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 0.993 |
| `po_minus_cov` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 0.991 |
| `po_gated` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 0.993 |
| `po_proto` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 0.991 |
| `po_delta` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 0.990 |
| `po_budget` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 0.996 |
| `po_next` | 0.6426 | 0.3574 | +0.0000 | +0.0000 | 1.000 | 0.995 |

### Opportunity ranking (vs `equal`)

1. `po_delta` opp_score=-0.0170 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
2. `po_proto` opp_score=-0.0170 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
3. `po_minus_cov` opp_score=-0.0171 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
4. `po_soft` opp_score=-0.0171 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
5. `po_gated` opp_score=-0.0171 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)

## coco

mods = `['img', 'txt']`

| version | Acc↑ | MSE↓ | Acc lift vs equal | MSE drop vs equal | FLOPs_rel | H(α) |
|---|---:|---:|---:|---:|---:|---:|
| `equal` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 1.000 |
| `po_soft` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 0.992 |
| `po_minus_cov` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 0.991 |
| `po_gated` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 0.992 |
| `po_proto` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 0.990 |
| `po_delta` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 0.987 |
| `po_budget` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 0.996 |
| `po_next` | 0.9023 | 0.0977 | +0.0000 | +0.0000 | 1.000 | 0.999 |

### Opportunity ranking (vs `equal`)

1. `po_delta` opp_score=-0.0168 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
2. `po_proto` opp_score=-0.0170 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
3. `po_minus_cov` opp_score=-0.0171 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
4. `po_soft` opp_score=-0.0171 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)
5. `po_gated` opp_score=-0.0171 (acc_lift=+0.0000, mse_drop=+0.0000, flops=1.000)

## Cross-dataset opportunity votes

| version | vote mass (top-5 ranks) |
|---|---:|
| `po_proto` | 17 |
| `po_delta` | 14 |
| `po_gated` | 10 |
| `po_soft` | 10 |
| `po_minus_cov` | 9 |

## Next-step training opportunities (actionable)

### A. When PO-risk should change the train plan

1. **Asymmetric multi-mod packs (Affec-like, M≥3)** — largest headroom.
   Softmax(PO) peaking lets `freeze_mask` drop low-α towers → FLOPs↓ with Acc≈flat.
   Best smoke winners: `po_gated` (FLOPs saver), `po_proto` (small Acc lift).
2. **Balanced img/txt packs already easy (COCO ~0.90 Acc)** — LR reallocation is near-noop.
   Keep `equal` or `po_budget` (warm floor); do not freeze; spend inventiveness on stack prior only.
3. **Non-stationary streams** — prefer `po_delta` / `po_next` so next-window LR moves *before*
   Acc collapses; pair with higher `gamma_delta` when ΔPO is reliable.

### B. Actuator recipes for the next training loop

| situation | metric | LR | steps | freeze | stack prior |
|---|---|---|---|---|---|
| concept spike on one mod | `po_soft` / `po_proto` | β≈0.2, boost top-α | dump steps to top-1/2 | freeze α<θ | KL→α |
| high MMD, flat uni-acc | `po_minus_cov` / `po_gated` | damp shared LR | cut total_steps | freeze gated-off | weak KL |
| need always-on towers | `po_budget` | soft only | equal-ish | never | KL→α |
| forecast next risk | `po_next` | use forecast α | from forecast | optional | KL→forecast α |
| already saturated Acc | `equal` | flat | flat | off | optional |

### C. What not to do

- Do **not** backprop into α (external attributor stays PO/FSDS/VIMP).
- Do **not** treat KL(stack_w‖α) as MoE aux — prior is external, `w` is the only learned mixer.
- Do **not** hard-gate FWD; freeze is BWD-only (efficiency = update FLOPs).
- On 2-mod packs, freeze_theta must be high (~0.35+) to ever fire; soft LR is the real lever.

### D. Suggested next training adjustments (priority)

1. Affec online loop: switch default metric `equal` → `po_gated` (or `po_proto` if Acc-first).
2. Wire `step_alloc` into real optimizer step counts (not just LR scale) — biggest unused lever.
3. Add holdout Acc@FLOPs Pareto (freeze_theta sweep) before claiming efficiency wins.
4. For Food-101 / Fashion-IQ: keep stacking KL prior = α, but keep freeze off; try `po_budget`.
5. Log per-window (PO, α, lr_mult, freeze, Acc) — use to tune τ / λ_cov / freeze_theta.
