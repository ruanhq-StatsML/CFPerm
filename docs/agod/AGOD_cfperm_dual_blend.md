# CFPerm-gated dual/blend (v8.2 beijing_gate retune)

Gate = **CFPerm DRPerm** (`e_mode=known`). This hour: **retune `beijing_gate`** so dual soft path fires.

## Problem (v8.1)

CFPerm reject intensities typically sit in **~0.15–0.35**. Old `beijing_gate=0.45` → `beijing_frac≲0.07` → dual almost always = hard_support → pack MSE hurt on metro/beijing.

## Fix

- Remap `(λ, is_beijing)` from **cached intensity** (`remap_gate_temper`) — no re-run of CFPerm.
- Scan `beijing_gate ∈ {0.15,0.20,0.25,0.30,0.35,0.45}` on metro / beijing / AAPL / MSFT (16×96, seeds 0/1).
- Lock default **`beijing_gate=0.25`**.

## Scan summary (dual pack MSE vs uniform, CFPerm-sig)

| beijing_gate | bj_frac@reject̄ | dual rel pack vs unī |
|---|---:|---:|
| 0.15 | 1.00 | **+0.6%** |
| 0.20 | 1.00 | **+0.6%** |
| **0.25** | **1.00** | **+0.6%** |
| 0.30 | 0.33 | +5.8% |
| 0.35 | 0.33 | +5.8% |
| 0.45 (old) | 0.20 | **+9.0%** |

At ≤0.25: reject intensities open soft CV (`cv_PO^0` / `cv_PO^0.25`); pack damage collapses vs hard-only. At ≥0.30: mostly hard_support again.

Per-pack at `beijing_gate=0.25` (mean seeds):

| dataset | duty | bj@reject | dual fam | dual pack vs uni (sig) |
|---|---:|---:|---|---:|
| `metro_interstate` | 0.17 | 1.00 | `cv_PO^0` | ~+1% (vs ~+13% hard) |
| `beijing_pm25` | 0.17 | 1.00 | `cv_PO^0/0.25` | ~+1% (vs ~+8% hard) |
| `stocks_MSFT` | 0.03 | 1.00 | `cv_PO^0.25` | **−0.7%** |
| `stocks_AAPL` | 0.00 | — | — | calm |

## Locked recipe update

| Item | Choice |
|---|---|
| L0 | CFPerm DRPerm, `e_mode=known`, `n_perm≥39` |
| L1 beijing | **`beijing_gate=0.25`** (was 0.45) |
| L1 temper | gate=0.20, lam_max=0.75 |
| Hard claim | `hard_m` still the hard-only path |
| Dual | mild→hard_support; intensity>0.25→soft CV |

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_cfperm_dual_blend.py \
  --beijing-gate 0.25 \
  --scan-beijing-gates 0.15 0.20 0.25 0.30 0.35 0.45 \
  --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT \
  --n-batches 16 --batch-size 96 --seeds 0 1 --skip-synth
```

Next: full 40×256 multi-seed under locked 0.25; optional dual+blend_50 hybrid only on BJ rejects.
