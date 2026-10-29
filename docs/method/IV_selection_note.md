# Instrumental Variable selection (synthetic argument)

Not “important-variable” VIMP — **instrumental variables** `Z` for
endogenous treatment under hidden confounding.

## Claim

When `U → (T, Y)` is unobserved, naive uplift / OLS for `τ = ∂Y/∂T` is
biased. A **valid** IV (`Z → T`, exclusion from `Y`) + 2SLS recovers `τ`.
Selection must pass **strength** (partial F) and an **exclusion screen**.

## Synthetic DGP

```
T = Xβ_t + Z_valid γ_v + Z_weak γ_w + Z_invalid γ_i + α U + ε_t
Y = τ T + Xβ_y + Z_invalid δ + α U + ε_y
```

| Z type | first stage | exclusion |
| --- | --- | --- |
| valid | strong | holds |
| invalid | weak leak | broken (direct `→ Y`) |
| weak | tiny | holds but useless |

## Selection (`iv_auto`)

1. partial F of each `Z_j` after residualizing on `X` — keep if `F ≥ 10`
2. exclusion residual `|corr(Z_j, Y | X, T̂, Z_{-j})|` — drop if large
3. 2SLS on the kept set

## What the suite shows (8 seeds, n=2000, α=1, τ=1.5)

| method | \|ATÊ−τ\| | ATÊ | notes |
| --- | --- | --- | --- |
| `naive_ols` | **0.301** | 1.80 | confounded uplift |
| `iv_oracle` | 0.017 | 1.51 | true valid Z |
| `iv_auto` | **0.030** | 1.52 | F≥10 + exclusion excess |
| `iv_invalid` | 2.29 | 2.31 | exclusion broken |
| `iv_weak` | 0.37 | 1.52 | weak first stage |

Selection quality (`iv_auto` vs true valid set): **P=0.94, R=0.96, F1=0.94**.

Confounding sweep: as `α` goes 0→2, naive error 0.06→0.63; `iv_auto` stays ~0.02–0.08.

## Sensitivity + subset localization (8 seeds, n=2000, α=1)

| tool | question | synthetic result |
| --- | --- | --- |
| **E-value** | how strong must U be (RR scale) to explain away the *naive* association? | point ≈ **5.1**, CI ≈ **5.0** |
| **Rosenbaum Γ\*** | odds-ratio bias at which matched high-vs-low contrast loses significance | Γ\* ≈ **10** (fairly robust) |
| **Subset localization** | does first-stage F vary across segments (X0 quartiles)? | Q1 F≈375 → Q4 F≈940 (**2.5×**) |

Designed heterogeneity: valid first stage scaled by `1 + het·σ(X0)` → Q4 stronger than Q1.
ATÊ stays near τ=1.5 in every quartile; **strength** is what moves, not the point estimate.

```bash
python3 scripts/run_iv_selection_synth.py
```

Code: `Python/src/iv_selection_synth.py`
