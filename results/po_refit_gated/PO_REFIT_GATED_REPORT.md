# Rolling PO-learner refit, gated re-adjustment

Uniform stays the default. When a new batch arrives, **refit** the
PO-learner: `T=0` = most recent control (`B_{t-2}`), `T=1` = 上一批 ∪ 这一批.
√PO weights apply **only if** mean PO-risk on T=1 / T=0 ≥ `gate=1.25`.

## Board (next-batch Ridge MSE)

| scene | uniform_pair | gated_pair | switch | adaptive | resid | resid_po | oracle | uniform_hop | gated_hop | dre_hop |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `similar` | 0.275 | 0.277 | 0.276 | 0.288 | 0.272 | 0.275 | 0.272 | 0.294 | 0.307 | 0.306 |
| `covariate` | 0.336 | 0.338 | 0.341 | 0.358 | 0.391 | 0.410 | 0.331 | 0.428 | 0.445 | 0.846 |
| `concept` | 1.674 | 1.460 | 1.635 | 1.628 | 1.372 | 1.372 | 1.378 | 1.376 | 1.368 | 1.366 |
| `mixed` | 1.963 | 1.748 | 1.843 | 1.852 | 1.616 | 1.628 | 1.624 | 1.616 | 1.625 | 1.782 |

Fire rates (PO-ratio γ=1.25 vs residual γ=2.0):

| scene | gated_pair | gated_hop | switch | resid | resid_po | soft_pair |
|---|---:|---:|---:|---:|---:|---:|
| `similar` | 0.02 | 0.15 | 0.15 | 0.02 | 0.02 | 0.02 |
| `covariate` | 0.02 | 0.19 | 0.19 | 0.21 | 0.21 | 0.02 |
| `concept` | 0.17 | 0.10 | 0.10 | 0.19 | 0.19 | 0.17 |
| `mixed` | 0.19 | 0.17 | 0.17 | 0.27 | 0.27 | 0.19 |

## Concept stream, MSE by hop (predict `B_{t+1}`)

Cut is `concept_at=4`. Predicting `B_4` uses hop `t=3` (both
pre-cut) — that first post-shift batch is structurally
unpredictable from labels. Gate can fire from `t=4`.

| t (train) | test | uniform_pair | switch | resid | resid_po | oracle | gated_hop |
|---|---|---:|---:|---:|---:|---:|---:|
| 1 | `B_2` | 0.283 | 0.266 | 0.262 | 0.262 | 0.262 | 0.294 |
| 2 | `B_3` | 0.284 | 0.284 | 0.284 | 0.284 | 0.284 | 0.298 |
| 3 | `B_4` ← cut | 6.833 | 6.807 | 6.834 | 6.813 | 6.833 | 6.715 |
| 4 | `B_5` | 2.104 | 1.911 | 0.311 | 0.332 | 0.311 | 0.316 |
| 5 | `B_6` | 0.279 | 0.284 | 0.279 | 0.279 | 0.290 | 0.300 |
| 6 | `B_7` | 0.259 | 0.259 | 0.259 | 0.259 | 0.286 | 0.286 |

- **similar / covariate**: `P(Y|X)` stable → uniform slightly better;
  always-on √PO hurts. PO-ratio gate fire-rate should stay low.
- **PO-ratio vs residual gate**: a *global* concept flip raises φ² on
  both T=0 and T=1, so ρ=r̄1/r̄0 stays near 1 and switch/adaptive
  keep pooling. Residual gate = MSE(fit 上一批 → 这一批) / train MSE;
  that is the 'batches clearly different' detector.
- **switch vs resid**: same train-set idea; differ only in the gate.
- **resid vs resid_po**: drop old batch, then optional √PO on the new one.
- **oracle**: knows `concept_at`; train-set upper bound, not deployable.
- **dre_hop**: X-only density ratio. Misses label shift; overreacts to X-hop.
- pair = last two batches; hop = new batch only.

11 methods, Ridge next-batch MSE, gate γ=1.25.
