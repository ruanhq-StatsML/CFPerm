# DiffusionDB temporal attribution (prompt tokens → image_nsfw)

## Definition
- **X** = lightweight prompt token embedding (`tfidf`, d=512)
- **Y** = `image_nsfw` (image attribute from metadata; no PNG download)
- **T** = 2 equal-count time windows (early T=0 → late T=1)

Span: `2022-08-06 22:02:00+00:00` → `2022-08-20 13:30:00+00:00` (2M gallery; short calendar span).

## Signed direction (outcome)
- Ȳ_early=0.1671, Ȳ_late=0.2178, ΔȲ=0.0507 → **positive/正向**

## Methods
1. **Covariate drift**: RF Domain VIMP (X→T)
2. **Tip-cmean**: signed μ_late−μ_early per token feature
3. **FSDS**: Scaler→Var→SelectKBest→HGB/LR predicting high-Y (threshold from early quantile); train early / test late
4. **Blend rank**: mean of percentile ranks (f_score, vimp_cov, |Δ|)

- FSDS HGB AUC(late)=0.630 AP=0.432
- FSDS LogReg AUC(late)=0.630 AP=0.441

## Top blend tokens

| rank | token | f_score | vimp_cov | Δ | sign |
|---:|---|---:|---:|---:|---:|
| 1 | `intricate` | 36.18 | 0.0144 | 0.0064 | 1 |
| 2 | `sharp` | 42.98 | 0.0090 | 0.0043 | 1 |
| 3 | `focus` | 31.90 | 0.0117 | 0.0062 | 1 |
| 4 | `greg` | 21.57 | 0.0159 | 0.0071 | 1 |
| 5 | `sharp focus` | 34.14 | 0.0074 | 0.0043 | 1 |
| 6 | `greg rutkowski` | 14.27 | 0.0187 | 0.0070 | 1 |
| 7 | `art` | 25.50 | 0.0066 | 0.0057 | 1 |
| 8 | `the` | 15.30 | 0.0105 | -0.0097 | -1 |
| 9 | `rutkowski` | 13.95 | 0.0126 | 0.0075 | 1 |
| 10 | `art by` | 27.45 | 0.0043 | 0.0042 | 1 |
| 11 | `and` | 17.66 | 0.0074 | 0.0040 | 1 |
| 12 | `by` | 6.22 | 0.0151 | 0.0087 | 1 |
| 13 | `portrait` | 22.13 | 0.0050 | 0.0027 | 1 |
| 14 | `highly` | 18.11 | 0.0059 | 0.0026 | 1 |
| 15 | `fantasy` | 17.36 | 0.0047 | 0.0030 | 1 |

## Note
Full CLIP token grid (77×768) can replace the light embedding later; protocol (T windows, FSDS, cmean sign, domain VIMP) stays the same.

Artifacts under `results/diffusiondb_temporal_fsds/`.
