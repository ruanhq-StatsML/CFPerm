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

Human feature-blocks → Stage-1 consensus π → Stage-2 training guided by π (π-weighted RF via column replicates; opinion pool of block RFs; X+Z stacking as the feature view of the same weights).

| setting | AUC raw | π-RF | pool | X+Z | Δ π-RF | π on GT | π-RF on GT |
|---|---:|---:|---:|---:|---:|---:|---:|
| inject valence (α=0.70) | 0.826 | 0.916 | 0.837 | 0.892 | +0.091 | 0.373 | 0.482 |
| text⊥ + inject valence | 0.855 | 0.941 | 0.865 | 0.911 | +0.085 | 0.466 | 0.572 |
| synthetic GT=valence | 0.880 | 0.897 | 0.901 | 0.902 | +0.016 | 0.734 | 0.828 |

Observational 1750/1950 (no GT, diffuse shift): raw AUC 0.720 vs π-RF 0.714 vs π-pool 0.703.

Clever-Z is `Z_m = π_m · logit ê_m(X_m)` (n_modalities columns), estimated OOF.
Stage-2 uses the same π to *train*: π-weighted RF (column replicates so sklearn max_features samples p_j ∝ π_m/|B_m|), opinion pool of block RFs, and X+Z stacking. Random-subspace BAWF is reported as an architecture ablation (it is not comparable to RF AUC). Adaptive logit can saturate on mean-shift inject.
Detection uses 5-fold stratified CV. GT rows average 3 subsample seeds.
Observational Chronoberg is a diffuse multi-modality shift; wins are on concentrated GT.

```bash
python3 scripts/run_chronoberg_clever_cov.py
```
