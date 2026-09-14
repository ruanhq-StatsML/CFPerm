# CFPerm-gated dual/blend (v8.3 full multi-seed)

Gate = **CFPerm DRPerm** (`e_mode=known`), **`beijing_gate=0.25`**.

## This hour

**Full 40×256 × 6 packs × seeds {0,1,2}** under locked v8.2 threshold.

## Synthetic L0

- **size = 0.0**, **power = 0.8** (`n_perm=39`)

## Stream packs (mean over 3 seeds)

| dataset | duty | intensitȳ | bj@reject | dual fam |
|---|---:|---:|---:|---|
| `metro_interstate` | 0.44 | 0.24 | **1.00** | `cv_PO^0.25` |
| `beijing_pm25` | 0.19 | 0.17 | **1.00** | `cv_PO^0` |
| `stocks_AAPL` | 0.03 | 0.11 | **1.00** | `cv_PO^0` |
| `stocks_MSFT` | 0.09 | 0.15 | **1.00** | `cv_PO^0` |
| `stocks_IWM` | 0.01 | 0.11 | **1.00** | — |
| `waymo_proxy` | 0.09 | 0.10 | **1.00** | `cv_PO^0` |

Soft path fires on essentially all rejects (`bj@reject=1`).

## Rel. pack MSE vs uniform (CFPerm-sig)

| dataset | hard_m | dual | b50 | duty |
|---|---:|---:|---:|---:|
| `metro_interstate` | +5.2% | **+1.8%** | +4.0% | 0.44 |
| `beijing_pm25` | +8.0% | **+0.4%** | +8.6% | 0.19 |
| `stocks_AAPL` | −1.1% | +0.3% | **−1.2%** | 0.03 |
| `stocks_MSFT` | +1.3% | **−0.8%** | −0.2% | 0.09 |
| `stocks_IWM` | +9.9% | **+1.4%** | +8.7% | 0.01 |
| `waymo_proxy` | +5.3% | **−0.2%** | +6.0% | 0.09 |

**Takeaway:** with `beijing_gate=0.25`, **dual ≪ hard_m** on high-duty packs (metro/beijing), and dual is near-uniform or better on most packs. Soft CV is doing the work v8.1 lacked.

## Sig-only hard top-20% next-MSE (↓)

| dataset | uni | hard_m | dual | b50 | best |
|---|---:---|---:---|---:---|---:|---|
| `metro_interstate` | **2.029e+06** | 2.061e+06 | 2.047e+06 | 2.060e+06 | uni |
| `beijing_pm25` | 3152 | 3369 | **3146** | 3349 | dual |
| `stocks_AAPL` | 0.001854 | 0.001888 | 0.001853 | **0.001827** | b50 |
| `stocks_MSFT` | 0.002170 | 0.002156 | **0.002123** | 0.002144 | dual |
| `stocks_IWM` | 0.000509 | 0.000527 | **0.000498** | 0.000532 | dual |
| `waymo_proxy` | 0.02224 | **0.02198** | 0.02218 | 0.02201 | hard_m |

Hard-subset wins: dual=3, hard_m=1, b50=1, uni=1.

## Locked recipe (v8.3)

| Item | Choice |
|---|---|
| L0 | CFPerm DRPerm, `e_mode=known`, `n_perm≥39` |
| L1 beijing | **`beijing_gate=0.25`** (validated 40×256) |
| Default policy | **`dual`** under CFPerm (soft on BJ intensity) |
| Hard-only claim | still report `hard_m` separately |
| Blend | mix≤0.5; useful on sparse AAPL pack |

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_cfperm_dual_blend.py \
  --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT stocks_IWM waymo_proxy \
  --n-batches 40 --batch-size 256 --n-perm 39 --seeds 0 1 2 \
  --beijing-gate 0.25 \
  --out results/agod_cfperm_dual_blend_v83
```

Next: dual+blend_50 hybrid on BJ-only; or post-reject FSDS/VIMP feature subset.
