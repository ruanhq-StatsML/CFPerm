# First-category performance (per dataset)

For each novel category, take the **earliest** stream indices (frac=0.25, n_first=6), then record accuracy on that cold-start window. Compare weights vs no-weights.

| dataset | uniform first_acc | weighted | both_fused | Δ weighted | Δ fused |
| --- | --- | --- | --- | --- | --- |
| har | 0.667 | 0.833 | 0.833 | 0.167 | 0.167 |
| mice | 0.500 | 0.250 | 0.250 | -0.250 | -0.250 |
| cnae9 | 0.125 | 0.250 | 0.250 | 0.125 | 0.125 |
| jvowel | 0.500 | 0.375 | 0.375 | -0.125 | -0.125 |
| wine | 1.000 | 1.000 | 1.000 | 0.000 | 0.000 |
| cancer | 1.000 | 1.000 | 1.000 | 0.000 | 0.000 |
| news20b | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| covtype | 0.500 | 0.500 | 0.500 | 0.000 | 0.000 |
| cgcd_stream | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

mean Δ weighted (tabular boards): **-0.010**
mean Δ both_fused (tabular boards): **-0.010**

## Per-category (both_fused) detail

### har

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 3 | 2 | 0.500 | 1.000 |
| 4 | 2 | 1.000 | 0.000 |
| 5 | 2 | 1.000 | 1.000 |

### mice

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 4 | 2 | 0.000 | 1.000 |
| 5 | 2 | 0.500 | 0.500 |
| 6 | 2 | 0.500 | 0.500 |
| 7 | 2 | 0.000 | 0.000 |

### cnae9

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 5 | 2 | 0.000 | 0.000 |
| 6 | 2 | 0.500 | 0.000 |
| 7 | 2 | 0.500 | 0.500 |
| 8 | 2 | 0.000 | 0.000 |

### jvowel

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 5 | 2 | 1.000 | 1.000 |
| 6 | 2 | 0.500 | 0.500 |
| 7 | 2 | 0.000 | 0.000 |
| 8 | 2 | 0.000 | 0.500 |

### wine

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 2 | 2 | 1.000 | 1.000 |

### cancer

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 1 | 2 | 1.000 | 1.000 |

### news20b

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 2 | 2 | 0.000 | 0.000 |
| 3 | 2 | 0.000 | 0.000 |

### covtype

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 4 | 2 | 1.000 | 1.000 |
| 5 | 2 | 0.000 | 0.000 |

### cgcd_stream

| category | n_first_eval | acc (fused) | acc (uniform) |
| --- | --- | --- | --- |
| 12 | 1 | 0.000 | 0.000 |
| 13 | 1 | 0.000 | 0.000 |
| 14 | 1 | 0.000 | 0.000 |
| 15 | 1 | 0.000 | 0.000 |
| 16 | 1 | 0.000 | 0.000 |
| 17 | 1 | 0.000 | 0.000 |
| 18 | 1 | 0.000 | 0.000 |
| 19 | 1 | 0.000 | 0.000 |
| 20 | 1 | 0.000 | 0.000 |
| 21 | 1 | 0.000 | 0.000 |
| 22 | 1 | 0.000 | 0.000 |
| 23 | 1 | 0.000 | 0.000 |

```bash
python3 scripts/run_first_category_performance.py --quick
```
