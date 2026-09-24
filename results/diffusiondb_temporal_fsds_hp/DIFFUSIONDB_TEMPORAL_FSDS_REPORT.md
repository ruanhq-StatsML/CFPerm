# DiffusionDB temporal FS (prompt → image_nsfw)

## Definition
- **X** = `tfidf` prompt tokens **⊕** cfg/step/sampler (d_tok=256, d_total=266)
- **Y** = continuous `image_nsfw` → early-quantile binary for FSDS
- **T** = `equal_count` (early T=0 → late T=1)
- window_meta: `{'scheme': 'equal_count', 'n_windows': 2, 'note': 'balanced n; calendar width varies by window'}`

Span: `2022-08-06 22:20:00+00:00` → `2022-08-20 12:54:00+00:00`.

## Signed direction
- Ȳ_early=0.1641, Ȳ_late=0.1973, ΔȲ=0.0331 → **positive/正向**
- n_by_T={0: 1000, 1: 1000}
- span_hours_by_T={0: 142.38333333333333, 1: 184.0}
- hp_vimp_share=0.087

## Methods (feature selection)
1. RF Domain VIMP (X→T)
2. cmean signed μ_late−μ_early
3. FSDS SelectKBest→HGB/LR on high-Y; train early / test late

Qwen embedding: optional later concat into X (deps not required here).

- HGB AUC(late)=0.602
- LogReg AUC(late)=0.608

## Top blend

| rank | feature | f | vimp | Δ |
|---:|---|---:|---:|---:|
| 1 | `hp_cfg` | 19.09 | 0.0701 | 0.2163 |
| 2 | `the` | 7.19 | 0.0180 | -0.0152 |
| 3 | `greg` | 10.77 | 0.0114 | 0.0095 |
| 4 | `hp_step` | 7.13 | 0.0159 | -0.1045 |
| 5 | `volumetric` | 13.93 | 0.0091 | 0.0060 |
| 6 | `greg rutkowski` | 7.24 | 0.0164 | 0.0083 |
| 7 | `rutkowski` | 8.31 | 0.0143 | 0.0081 |
| 8 | `art` | 5.59 | 0.0099 | 0.0106 |
| 9 | `portrait` | 3.40 | 0.0103 | 0.0089 |
| 10 | `background` | 3.84 | 0.0066 | 0.0071 |
| 11 | `and greg` | 14.24 | 0.0064 | 0.0041 |
| 12 | `and` | 2.04 | 0.0163 | 0.0094 |
| 13 | `digital art` | 5.93 | 0.0044 | -0.0084 |
| 14 | `focus` | 10.27 | 0.0050 | 0.0053 |
| 15 | `masterpiece` | 3.08 | 0.0081 | 0.0081 |
