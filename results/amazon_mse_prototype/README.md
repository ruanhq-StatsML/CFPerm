# Amazon MSE prototype

The only number is **online rating MSE**. Predict the arriving category, then update.

## Shapes (live, n_per=80 quick / 240 full)

Amazon:
- `X` `(N, P)` TF-IDF, vocab frozen on batch 0. `P ≤ 128` (Gift-card head often ~109).
- `y` `(N,)` stars in `[1, 5]`
- `batch` `(N,)` category id `0..8`
- `N = 9 * n_per`
- SDC prototypes `mu` `(5, P)`
- GPM bases `M` `(P, k)`

MSR-VTT is not scored here: `s` is `(n, 2049) = 768+512+768+1`.

## Methods

- `plateau` — RidgeSGD, current Amazon winner
- `bank` — Wu instance memory: cosine-weighted ratings from stored `(x,y)`, then enqueue
- `sdc` — Yu star prototypes (5 means). Kept as the prototype ablation; it does **not** cut MSE here
- `gpm_typed` — plateau RidgeSGD, project `dw` iff heatmap `ĉ` is loud and `δ̂` is quiet

```
python3 scripts/run_amazon_mse_prototype.py
python3 scripts/run_amazon_mse_prototype.py --quick
```

## Live numbers (4 seeds, n_per=240, X is `(2160, 128)`)

| method | online MSE |
| --- | --- |
| bank | 1.566 (0.063) |
| plateau | 1.590 (0.068) |
| gpm_typed | 1.593 (0.070) |
| sdc (5 star means) | 2.289 (0.022) |

The instance bank is the one that cuts MSE. Five star means do not. Typed GPM is a wash vs plateau. Do not oversell.
