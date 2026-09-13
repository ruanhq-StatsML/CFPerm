# How to tune (after embedding `benchmark_feature_selection`)

Code path is fixed. Tuning = 4 knobs + one swap hook.

## Pipeline (locked)

```
hop t-2 → t-1
   └─ benchmark_feature_selection(X0,X1)  → vimp
         └─ modality_mass(vimp)           → π_m
               └─ modality_hop_weights(π, γ) → sample weights for Ridge
PO / same vimp top-k                      → mask (zero|noise|missing) → ΔMSE
```

Amazon stays **uniform** `hop_ridge` (one modality ⇒ π collapses).

## Knobs (tune in this order)

| # | knob | where | start | what it does |
| --- | --- | --- | --- | --- |
| 1 | **γ** | `modality_hop_weights(..., gamma=)` | `4` | sharper IW → more weight on π-aligned past batches. Try `{2,4,8}` |
| 2 | **ewma_pi** | `modality_pi_shares(..., ewma=)` / `run_msrvtt_hop_po(ewma_pi=)` | `0.3` | smooth π across hops. `0` = raw each hop; `0.5–0.7` if π jumps |
| 3 | **α** | Ridge `alpha` | `3` | only if under/overfit on reconstruction MSE |
| 4 | **k_mask** | `run_mask_ablation(k_tokens=)` / `k_per_mod` | `16` / `32` | attribution resolution, not training MSE |

Do **not** tune feature reweighting by VIMP — already lost to `hop_ridge`.

## Selector swap (the real lever)

Default body of `benchmark_feature_selection` = RF-Domain VIMP.

```python
# Python/src/msrvtt_multimodal_attribution.py
def benchmark_feature_selection(X0, X1, seed=SEED, n_estimators=40):
    return rf_domain(X0, X1, seed=seed, n_estimators=n_estimators)
```

Swap candidates (same return `(vimp, meta)`):

1. **RF-Domain** (current) — covariate hop, fast
2. **PO-risk `po_tau_vimp`** — if you care about treatment/outcome-linked shift
3. **MMD / group LOCO** (`group_mmd_loco` block scores → expand to coord VIMP) — when RF is noisy
4. Pass `selector=` into `modality_pi_shares(...)` without editing the default

Judge by: (a) online cross-modal MSE lift vs uniform hop, (b) π_m path stability, (c) mask ΔMSE that matches intuition (video vs audio).

## Suggested sweep (MSR-VTT)

```text
for γ in {2,4,8}:
  for ewma_pi in {0.0, 0.3, 0.6}:
      run_msrvtt_hop_po(use_mod=True, gamma=γ, ewma_pi=ewma_pi)
# pick by online MSE; then freeze and report mask ΔMSE once
```

Amazon: leave γ=4 on `hop_ridge`; only re-check if TF-IDF dim / n_per changes.

## What “good” looks like

- Multimodal: hop+π ≤ hop (uniform), π_m not collapsing to 1/3 forever
- Mask: top-π modality’s zero/noise ΔMSE > 0 and larger than low-π modality
- Amazon: hop_ridge still ≈ 1.38; no regression from selector experiments
