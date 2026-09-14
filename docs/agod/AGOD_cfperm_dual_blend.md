# CFPerm-gated dual/blend (v8)

Gate = **CFPerm subset**, not OnlineRFPerm.

## How to estimate / evaluate (by layer)

| Layer | Estimate | Evaluate |
|---|---|---|
| **L0 Gate** | `DRPerm`: recent=`W=0`, current=`W=1`; PO-risk + permute-W → `(T, p, reject)` | Synthetic **size** (null≈α) / **power** (alt↑); stream **duty** (selective) |
| **L1 Intensity** | `intensity = 0.55·po_gap + 0.30·p_strength + 0.15·T_strength` → `(λ, beijing?)` | beijing_frac; λ̄; calm rejects stay near uniform |
| **L2 Shape** | `hard_support` / soft CV `PO^p` / `blend(mix)` | hard-rank (Spearman, P@20%); shape ≠ gate |
| **L3 Policy** | `dual` (mild→hard, BJ→soft CV) vs `blend_50` | **sig-only** hard-subset next-MSE (primary); pack MSE only under beijing |

`blend mix` = L2 convex combo under shared λ. `dual` = L3 regime switch — not a substitute for mix.

### Critical implementation notes

1. **Permutation floor**: need `1/(n_perm+1) < α` (use `n_perm≥39` for α=0.05).
2. **Stream propensity**: default `e_mode='known'` (e≡n₁/n). Fitted `e(X)` absorbs covariate-shift on adjacent packs → p≈1, duty=0.
3. Cache L0 once per timestep; reuse across L2/L3 modes.

## CFPerm subset inventory

| Piece | Role for AGOD stream |
|---|---|
| **DRPerm** | **Primary L0** — batch shift via PO-risk + permute-W |
| **RRPerm** | Optional L0 via R-risk (`--risk rr`) |
| **CFPerm-VIMP** (permuCATE / LOCO / GRF) | **Post-hoc** feature attribution after reject — not the stream gate |

Optional ablation only: Jaccard(CFPerm reject set, OnlineRFPerm reject set).

## Smoke (16×96, n_perm=39, known e)

| dataset | duty | intensitȳ | beijing_frac | hard-subset best | notes |
|---|---:|---:|---:|---|---|
| `metro_interstate` | 0.20 | 0.18 | 0.08 | uni | mild rejects; pack MSE not yet for hard |
| `beijing_pm25` | 0.13 | 0.17 | 0.00 | **hard_m** | hard claim direction OK |
| `stocks_AAPL` | 0.00 | 0.05 | 0.00 | uni | calm — gate correctly quiet |

Synthetic L0: **size=0.0, power=0.8**.

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_cfperm_dual_blend.py \
  --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT \
  --n-batches 40 --batch-size 256 --n-perm 39
```

API: `agod/cfperm_gate.py` (`cfperm_batch_test`, `cfperm_intensity_to_temper`, `synthetic_shift_trial`).

Legacy RFPerm dual/blend numbers: `docs/agod/AGOD_obs_po_weights.md` — secondary once CFPerm gate is primary.
