# Weighted cold-start under category shift

Weighted cold-start (upweight localized subset + selected features) vs uniform partial_fit on C-GCD proxy; Office-Home PDA checks incomparable rate rises as label overlap falls.

## CIFAR-100 C-GCD proxy (known → staged novel influx)

| method | mean known acc | mean novel acc | mean forgetting | feat hit |
| --- | --- | --- | --- | --- |
| weighted | 0.816 | 0.956 | 0.462 | 0.542 |
| uniform | 0.799 | 0.943 | 0.447 | — |

Δ novel acc (weighted − uniform): **0.012**

Δ forgetting (uniform − weighted; higher = weighted forgets less): **-0.015**

## Office-Home PDA proxy (partial label overlap)

| overlap | shared acc | incomparable rate | adapt rate |
| --- | --- | --- | --- |
| 0.1 | 0.767 | 1.000 | 0.000 |
| 0.4 | 0.939 | 0.000 | 0.333 |
| 0.8 | 0.804 | 0.000 | 0.667 |

Weights: sample `w=1+λ·membership` in top-K VIMP subspace; soft-freeze non-selected features in `partial_fit`.

Boundary: Jaccard(label_ref, label_live) < 0.2 → **incomparable** (skip directed adapt).

```bash
python3 scripts/run_weighted_cold_start.py --quick
```
