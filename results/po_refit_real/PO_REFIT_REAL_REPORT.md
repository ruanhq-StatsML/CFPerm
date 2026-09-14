# Real-data consecutive-batch PO-refit (no K-fold)

Each hop trains on the latest consecutive batches and scores the
**next** batch. Residual gate (no K-fold):

`ρ = err(fit B_{t-1} → B_t) / err(fit B_{t-2} → B_{t-1})`

Streams are time- or space-ordered: Metro Interstate traffic,
NYC green taxi, California latitude, diabetes readmission
source→target. Continuous tasks report **RMSE** (↓); discrete report **Acc** (↑).

## `rf`

| dataset | task | n_batches | uniform_pair | gated_pair | switch | resid | resid_po | uniform_hop | dre_hop | fire_resid |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `interstate` | RMSE | 12 | 742.5653 | 753.6823 | 749.8380 | 748.6585 | 794.8598 | 831.6054 | 919.6351 | 0.20 |
| `nyc_taxi` | RMSE | 12 | 3.0610 | 3.4956 | 2.9863 | 2.9796 | 6.4155 | 3.0234 | 3.1062 | 0.10 |
| `diabetes_readmit` | acc | 12 | 0.6590 | 0.6590 | 0.6520 | 0.6520 | 0.6520 | 0.6550 | 0.6550 | 0.00 |
| `california` | RMSE | 12 | 0.8222 | 0.8362 | 0.8307 | 0.8328 | 0.8557 | 0.8523 | 0.9210 | 0.10 |
| `diabetes` | RMSE | 8 | 62.3650 | 63.6214 | 64.1201 | 62.3974 | 62.3974 | 65.0094 | 72.6396 | 0.00 |
| `wine_red` | RMSE | 7 | 0.7196 | 0.7196 | 0.7108 | 0.7159 | 0.7159 | 0.7354 | 0.7703 | 0.00 |

## `xgb`

| dataset | task | n_batches | uniform_pair | gated_pair | switch | resid | resid_po | uniform_hop | dre_hop | fire_resid |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `interstate` | RMSE | 12 | 714.8927 | 727.5632 | 791.4930 | 729.3593 | 727.2737 | 875.4076 | 837.0926 | 0.20 |
| `nyc_taxi` | RMSE | 12 | 3.1940 | 3.9710 | 3.4359 | 3.4343 | 3.5016 | 3.3324 | 3.8698 | 0.20 |
| `diabetes_readmit` | acc | 12 | 0.6430 | 0.6430 | 0.6460 | 0.6460 | 0.6460 | 0.6245 | 0.6215 | 0.00 |
| `california` | RMSE | 12 | 0.7666 | 0.7673 | 0.7871 | 0.7579 | 0.7579 | 0.8307 | 0.8459 | 0.00 |
| `diabetes` | RMSE | 8 | 69.4372 | 69.5872 | 70.8420 | 68.9944 | 68.9944 | 72.4591 | 74.2521 | 0.00 |
| `wine_red` | RMSE | 7 | 0.7441 | 0.7537 | 0.7620 | 0.7413 | 0.7413 | 0.7759 | 0.7612 | 0.00 |

- **resid** drops the old batch only when the consecutive hop error jumps.
- **resid_po** = resid train-set + √PO on the new batch.
- **gated_pair / switch** use PO-ratio, which misses a global label map flip.
- No K-fold; no shuffle. Row order is the stream clock.
