# How to tune (after embedding `benchmark_feature_selection`)

Code path is fixed. Tuning = knobs below + selector swap.

## Pipeline (locked)

```
hop t-2 → t-1
   ├─ (ĉ, δ̂) from heatmap / concept intensity
   ├─ benchmark_feature_selection → vimp → modality_mass → π_m
   │     └─ attribution_ewma_rate(ĉ, δ̂, Δπ) → λ_new
   │           └─ EWMA persistence = 1−λ_new  (next-batch ensemble memory)
   └─ modality_hop_weights(π, γ) → sample weights for Ridge
PO / same vimp top-k → mask (zero|noise|missing) → ΔMSE
```

**Attribution → next-batch ensemble EWMA** (the useful bit):

| signal | λ_new (weight on *new* hop) | ensemble effect |
| --- | --- | --- |
| ĉ large / δ̂ quiet (covariate) | ↓ | remember past bank⊕Ridge mix |
| δ̂ large (concept) | ↑ | forget fast, track Ridge |
| π_m jumps | ↑ bump | modality mix moved → refresh |

Amazon `attr_adapter` and MSR-VTT `ewma_pi="auto"` both use `attribution_ewma_rate`.

Amazon training MSE still prefers uniform **`hop_ridge`**; `attr_adapter` is the ensemble path that *uses* this guide.

## Knobs (tune in this order)

| # | knob | where | start | what it does |
| --- | --- | --- | --- | --- |
| 1 | **γ** | `modality_hop_weights(..., gamma=)` | `4` | IW sharpness. Try `{2,4,8}` |
| 2 | **ewma_pi** | `run_msrvtt_hop_po(ewma_pi=)` | `"auto"` | `"auto"` = attribution-guided; or fix `{0,0.3,0.6}` |
| 3 | **α** | Ridge `alpha` | `3` | only if under/overfit |
| 4 | **k_mask** | mask top-k | `16` / `32` | attribution resolution |

Do **not** tune feature reweighting by VIMP — already lost to `hop_ridge`.

## Selector swap (the real lever)

Default body of `benchmark_feature_selection` = RF-Domain VIMP.

```python
# Python/src/msrvtt_multimodal_attribution.py
def benchmark_feature_selection(X0, X1, seed=SEED, n_estimators=40):
    return rf_domain(X0, X1, seed=seed, n_estimators=n_estimators)
```

Swap candidates (same return `(vimp, meta)`): RF-Domain / PO-τ / MMD-LOCO.
Or pass `selector=` into `modality_pi_shares(...)`.

Judge by: (a) online MSE vs uniform hop, (b) π path + λ_new stability, (c) mask ΔMSE.

## Suggested sweep (MSR-VTT)

```text
for γ in {2,4,8}:
  run_msrvtt_hop_po(use_mod=True, gamma=γ, ewma_pi="auto")
# optional ablate auto vs fixed:
for ewma_pi in {0.0, 0.3, 0.6, "auto"}: ...
```

Amazon: keep `hop_ridge` for scoreboard; inspect `attr_adapter` history `ewma_rate` to see the guide working.
