# Chronoberg clever-covariate prototype

Temporal batch **W = 1750 vs 1950** on Hugging Face `spaul25/Chronoberg` test books.
Modalities: hashed text n-grams, valence / arousal / dominance (pooled lexicons).

## Consensus relative contribution π

| modality | π |
|---|---:|
| text | 0.311 |
| valence | 0.307 |
| arousal | 0.131 |
| dominance | 0.250 |

## GT board (mean over seeds, 5-fold CV AUC)

| setting | AUC raw | AUC Z | AUC X+Z | Δ stack | π on GT |
|---|---:|---:|---:|---:|---:|
| inject valence (α=0.70) | 0.826 | 0.864 | 0.892 | +0.066 | 0.373 |
| text⊥ + inject valence | 0.855 | 0.880 | 0.911 | +0.055 | 0.466 |
| synthetic GT=valence | 0.880 | 0.883 | 0.902 | +0.021 | 0.734 |

Observational 1750/1950 (no GT, diffuse shift): raw AUC 0.720 vs clever-Z 0.694.

Clever-Z is `Z_m = π_m · logit ê_m(X_m)` (n_modalities columns), estimated OOF.
Detection uses 5-fold stratified CV. GT rows average 3 subsample seeds.

```bash
python3 scripts/run_chronoberg_clever_cov.py
```
