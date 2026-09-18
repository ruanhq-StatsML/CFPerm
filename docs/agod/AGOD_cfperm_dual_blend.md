# CFPerm-gated dual/blend (v8.4 dual_b50 hybrid)

Gate = **CFPerm DRPerm** (`e_mode=known`), **`beijing_gate=0.25`**.

## This hour

**Idea:** `dual_b50` — mild rejects → hard_support; Beijing intensity → fixed **blend_50** (hard/qrt mix=0.5) instead of soft CV.

Compare vs locked soft-CV **dual** on metro/beijing packMSE. Keep `hard_m` for hard claim.

## Locked carry-forward (v8.3)

| Item | Choice |
|---|---|
| L0 | CFPerm DRPerm, `e_mode=known`, `n_perm≥39` |
| L1 beijing | **`beijing_gate=0.25`** |
| Default policy (pending v8.4) | **`dual`** unless dual_b50 wins metro/beijing |
| Hard-only claim | still report `hard_m` |

## Modes

| mode | mild reject | beijing reject |
|---|---|---|
| `hard_m` | hard_support | hard_support |
| `dual` | hard_support | soft CV (PO^p) |
| `dual_b50` | hard_support | blend hard/qrt mix=0.5 |
| `blend_50` | hard_support | blend hard/qrt mix=0.5 |

Note: under `bj@reject≈1`, `dual_b50` ≈ `blend_50` on packMSE; difference is mild-path framing only when mild rejects exist.

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_cfperm_dual_blend.py \
  --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT stocks_IWM waymo_proxy \
  --n-batches 40 --batch-size 256 --n-perm 39 --seeds 0 1 2 \
  --beijing-gate 0.25 \
  --out results/agod_cfperm_dual_blend_v84
```

Results pending bench.
