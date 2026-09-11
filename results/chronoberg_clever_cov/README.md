# Chronoberg clever-covariate prototype

Temporal batch **W = 1750 vs 1950** on Hugging Face `spaul25/Chronoberg` test books.
Modalities: hashed text n-grams, valence / arousal / dominance (pooled lexicons).

## Consensus relative contribution π

| modality | π |
|---|---:|
| text | 0.282 |
| valence | 0.343 |
| arousal | 0.109 |
| dominance | 0.266 |

## Raw X vs clever-Z vs stack X+Z

| setting | AUC raw | AUC Z | AUC X+Z | Δ stack |
|---|---:|---:|---:|---:|
| observational 1750/1950 | 0.739 | 0.692 | 0.690 | -0.049 |
| CD: Y=valence, X=text+A+D | 0.710 | 0.637 | 0.638 | -0.072 |
| inject valence (α=0.95) | 0.859 | 0.816 | 0.875 | +0.015 |
| inject text (α=0.75) | 0.999 | 1.000 | 1.000 | +0.001 |
| text-shuffled + inject valence | 0.888 | 0.816 | 0.913 | +0.025 |
| synthetic GT=valence | 0.770 | 0.879 | 0.944 | +0.175 |

Clever-Z is the compact detector: instance-level relative contributions `π_m(x)`
and `logit ê_m(X_m)` (2 × n_modalities columns). TMLE `H_m` is used only in PO-risk targeting.

```bash
python3 scripts/run_chronoberg_clever_cov.py
```
