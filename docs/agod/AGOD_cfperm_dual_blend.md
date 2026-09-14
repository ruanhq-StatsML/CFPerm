# CFPerm-gated dual/blend (v8.1 multi-seed)

Gate = **CFPerm DRPerm** (`e_mode=known`), not OnlineRFPerm.

## This hour’s idea

**Multi-seed mean (seeds 0/1/2) + Jaccard(CFPerm, RFPerm) reject sets** on 6 packs (24×128, `n_perm=39`).

## How to estimate / evaluate (by layer)

| Layer | Estimate | Evaluate |
|---|---|---|
| **L0 Gate** | DRPerm: recent=`W=0`, current=`W=1`; PO-risk + permute-W → `(T,p,reject)` | Synth size/power; stream **duty** |
| **L1 Intensity** | `0.55·po_gap + 0.30·p_str + 0.15·T_str` → `(λ, beijing?)` | beijing_frac; calm ≈ uniform |
| **L2 Shape** | hard_support / soft CV / blend(mix) | hard-rank; shape ≠ gate |
| **L3 Policy** | dual vs blend_50 | **sig-only hard-subset** primary; pack MSE under beijing |
| **Ablation** | — | **Jaccard(CFPerm, RFPerm)** reject sets |

### Implementation locks

1. `n_perm ≥ 39` for α=0.05 (else min p ≥ α).
2. Stream propensity **`e_mode=known`** (fitted e collapses duty on adjacent packs).
3. Cache L0 once per timestep across modes.
4. Multi-seed: average metrics; report mean duty / Jaccard.

## CFPerm subset inventory

| Piece | Role |
|---|---|
| **DRPerm** | Primary L0 |
| **RRPerm** | Optional L0 (`--risk rr`) |
| **CFPerm-VIMP** | Post-hoc only |

## Synthetic L0 (from smoke)

- **size = 0.0**, **power = 0.8** (`n_perm=39`, concept-drift DGP)

## 6-pack multi-seed (24×128, seeds 0/1/2)

### Duty / Jaccard

| dataset | CFPerm duty | intensitȳ | beijing_frac | Jaccard vs RFPerm | RFPerm duty |
|---|---:|---:|---:|---:|---:|
| `metro_interstate` | 0.42 | 0.25 | 0.07 | **0.22** | 0.23 |
| `beijing_pm25` | 0.13 | 0.17 | 0.05 | **0.18** | 0.14 |
| `stocks_AAPL` | 0.01 | 0.10 | 0.02 | **0.00** | 0.16 |
| `stocks_MSFT` | 0.03 | 0.11 | 0.03 | **0.07** | 0.20 |
| `stocks_IWM` | 0.04 | 0.14 | 0.03 | **0.07** | 0.20 |
| `waymo_proxy` | 0.04 | 0.11 | 0.00 | **0.08** | 0.35 |

**Takeaway:** Jaccard ≪ 1 → CFPerm and OnlineRFPerm open on **different** batches; RFPerm tables are not a drop-in substitute.

### Sig-only hard top-20% next-MSE (↓)

| dataset | uni | hard_m | dual | b50 | best |
|---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | **3.217e+06** | 3.241e+06 | 3.231e+06 | 3.242e+06 | uni |
| `beijing_pm25` | **6829** | 7472 | 7476 | 7443 | uni |
| `stocks_AAPL` | **0.001686** | 0.001709 | 0.001709 | 0.001709 | uni |
| `stocks_MSFT` | 0.0009212 | **0.0008966** | 0.0008966 | 0.0008966 | hard_m |
| `stocks_IWM` | 0.000429 | 0.0004164 | 0.0004206 | **0.0004156** | b50 |
| `waymo_proxy` | 0.02101 | **0.01958** | 0.01958 | 0.01958 | hard_m |

Hard wins: hard_m=2, b50=1, uni=3 (this shorter multi-seed bench).

### Rel. pack MSE vs uniform (CFPerm-sig)

| dataset | hard_m | dual | b50 | duty |
|---|---:|---:|---:|---:|
| `metro_interstate` | +6.5% | +6.0% | +6.2% | 0.42 |
| `beijing_pm25` | +21.8% | +19.9% | +18.5% | 0.13 |
| `stocks_AAPL` | **−4.0%** | −4.0% | −4.0% | 0.01 |
| `stocks_MSFT` | **−2.9%** | −2.9% | −2.9% | 0.03 |
| `stocks_IWM` | +1.6% | +0.3% | +0.2% | 0.04 |
| `waymo_proxy` | +2.2% | +2.2% | +2.2% | 0.04 |

On low-duty stock packs, hard_m helps pack MSE; on high-duty metro/beijing, reweighting hurts pack MSE (beijing soft path barely fires — beijing_frac≤0.07).

## Locked recipe (partial, CFPerm era)

- **L0** = CFPerm DRPerm, `e_mode=known`, `n_perm≥39`
- **Hard claim path** still `hard_m` when duty is sparse (stocks)
- **dual / blend** not yet re-winning beijing pack under CFPerm reject sets (beijing_frac too low on this bench) → next: intensity→beijing threshold / temper schedule / full 40×256
- Do **not** treat OnlineRFPerm v6/v7 numbers as primary once CFPerm gate is on

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_cfperm_dual_blend.py \
  --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT stocks_IWM waymo_proxy \
  --n-batches 40 --batch-size 256 --n-perm 39 --seeds 0 1 2 --jaccard-rfperm
```

API: `agod/cfperm_gate.py`. Legacy RFPerm dual/blend: `docs/agod/AGOD_obs_po_weights.md`.
