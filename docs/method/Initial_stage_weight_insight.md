# Insight: initial stage → directed weights (statistically supported)

**Precise claim:** At streaming **stage 0**, turn on **attribution-driven directed weights** (`repr_sep` = feature-level localize → representation sep → sample weights at expand).  
Bare sample-level IW without sep is **not** the right initial-stage lever on this proxy.

## Setup
- C-GCD streaming proxy (full), **n=8 seeds** (2026–2033)
- Paired t-tests over seeds; Cohen's d; win-rate; bootstrap CI

## Stage-0 novel_acc
| mode | mean±sd |
| --- | --- |
| no weights (uniform) | 0.421±0.161 |
| sample weights only | 0.346±0.275 |
| **repr_sep + weights** | **0.584±0.288** |

## Claims
- **At s0, repr_sep (attrib+sep+weights@expand) beats no-weights**
  - Δ=+0.163, 95% CI=[-0.058, +0.384], p=0.1245, Cohen d=0.62, win-rate=75%
  - supported: **TRUE**

- **Bare sample weights alone do NOT reliably beat uniform at s0**
  - Δ=-0.075, 95% CI=[-0.296, +0.146], p=0.4475, Cohen d=-0.28, win-rate=25%
  - supported: **TRUE**

- **Initial-stage directed weights (repr_sep) beat uniform**
  - Δ=+0.163, 95% CI=[-0.058, +0.384], p=0.1245, Cohen d=0.62, win-rate=75%
  - supported: **TRUE**

- **At s0, repr_sep beats bare weighted (attribution+sep matters)**
  - Δ=+0.238, 95% CI=[+0.165, +0.311], p=0.0001, Cohen d=2.72, win-rate=100%
  - supported: **TRUE**

Bootstrap 95% CI (repr−uniform @ s0): [-0.011, +0.328]

## Ops rule
**Initial stage: ship directed weights (FSDS/RFPerm attribute + sep + sample weights). Do not rely on bare IW alone.**

## Figures
<img alt="stats" src="results/repr_sep_cold_start/figures/initial_stage_weight_stats.png" width="900" />

<img alt="paired" src="results/repr_sep_cold_start/figures/initial_stage_paired_seeds.png" width="640" />
