# Confirm + tip-cmean formulation (prototype)

## Joint fire
`Fire = 1{Δ_tip ≥ ε_Δ} · 1{confirm}` with
`confirm = 1{collapse ≥ thr_col} · 1{PO ≥ thr_PO}`.

## How magnitude thresholds are set
- recipe: **`burn_mean_plus_k_sd`** with `k=5.0`
- `ε_Δ = mean(Δ_tip_burn) + k·sd = 0.1322 + 5.0·0.0485 = **0.3750**`
- `thr_collapse = 0.0090 + k·0.0162 = **0.0901**`
- `thr_PO = 0.2307 + k·0.0347 = **0.4044**`

Burn-only calibration; engineering gates (not Type-I α).

## Result
- true tips [0, 1]: t*=20 delay=0
- wrong tips (6,7): t*=None delay=None
- signed tip shifts at fire: `{0: 2.439173048460785, 1: 2.567120775359228}`
