# Gap-guided reading budget (not TMLE, not causal)

Two-sample localization of W on X. **Same expert**: one m̂=argmax π for every query.
**Context packing**: one reader on concatenated columns of E=E(π). Not an opinion pool.

## Same-expert specialist AUC (ê_m, one block, all rows)

| setting | oracle | π | VIMP | instance s(x) | same random |
|---|---:|---:|---:|---:|---:|
| synthetic GT=valence | 0.959 | 0.959 | 0.959 | 0.931 | 0.466 |
| inject valence | 0.747 | 0.747 | 0.600 | 0.698 | 0.664 |
| diffuse equal shift | nan | 0.716 | 0.698 | 0.694 | 0.727 |
| shuffle valence (neg.) | 0.505 | 0.474 | 0.470 | 0.501 | 0.498 |

## Context packing AUC (one RF on concat E)

| setting | pack GT | pack π | pack VIMP | concat top-2(π) | pack all |
|---|---:|---:|---:|---:|---:|
| synthetic GT=valence | 0.960 | 0.960 | 0.960 | 0.952 | 0.947 |
| inject valence | 0.753 | 0.753 | 0.600 | 0.755 | 0.767 |
| diffuse equal shift | nan | 0.710 | 0.698 | 0.780 | 0.876 |
| shuffle valence (neg.) | 0.505 | 0.459 | 0.475 | 0.474 | 0.496 |

## Stage B: query update of E (π, f frozen)

| setting | prior AUC | shrink AUC | query-only AUC | switch rate | Δ AUC |
|---|---:|---:|---:|---:|---:|
| synthetic GT=valence | 0.959 | 0.957 | 0.931 | 0.101 | -0.002 |
| inject valence | 0.747 | 0.745 | 0.698 | 0.414 | -0.002 |
| diffuse equal shift | 0.716 | 0.725 | 0.694 | 0.286 | +0.009 |
| shuffle valence (neg.) | 0.474 | 0.481 | 0.501 | 0.641 | +0.007 |

```bash
python3 scripts/run_gap_guided_inference.py
```
