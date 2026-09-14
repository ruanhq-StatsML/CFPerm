# Rolling PO-learner refit, gated re-adjustment

Uniform stays the default. When a new batch arrives, **refit** the
PO-learner: `T=0` = most recent control, `T=1` = 上一批 ∪ 这一批.
√PO weights apply **only if** mean PO-risk on T=1 / T=0 ≥ `gate=1.25`.

| scene | uniform_pair | gated_pair | always_pair | adaptive | fire_adp | uniform_hop | gated_hop |
|---|---:|---:|---:|---:|---:|---:|---:|
| `similar` | 0.273 | 0.277 | 0.333 | 0.293 | 0.17 | 0.289 | 0.306 |
| `covariate` | 0.370 | 0.373 | 0.449 | 0.406 | 0.21 | 0.517 | 0.538 |
| `concept` | 1.695 | 1.515 | 1.409 | 1.609 | 0.12 | 1.410 | 1.399 |
| `mixed` | 2.176 | 2.003 | 1.884 | 1.983 | 0.21 | 1.844 | 1.885 |

- **similar**: batches share P(Y|X) → uniform should be slightly better; gated fire-rate should stay low. Always-on √PO is worse.
- **concept / mixed**: a clear hop → drop the old batch (hop / adaptive). Pair-uniform mixes pre- and post-shift and is worst on concept.
- pair = train on last two; hop = new batch only; adaptive = pair+uniform when quiet, hop+√PO when the hop contrast fires.

4 seeds, 8 batches × 100, Ridge next-batch MSE. Gate γ=1.25.
