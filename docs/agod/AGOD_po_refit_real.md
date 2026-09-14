# Real consecutive-batch PO-risk adaptation

Overview: `docs/agod/AGOD_overview.md`.
LaTeX: `docs/agod/AGOD_po_refit_real_tables.tex`,
`docs/agod/AGOD_performance_tables.tex`.

Batch size **200**, row order is the stream clock (no shuffle,
no K-fold). Gate γ=1.5: fire only when consecutive OOS probe
error jumps. Always-on DRE / always-on √PO are off this board.

- **uniform_pair**: last two batches, w=1.
- **rfperm**: same rows; on fire, T=1 gets `w=√po_risk0`.
- **resid**: residual hop-gate (drop old batch) as a reference.

Continuous tasks report **RMSE** (↓); discrete report **Acc** (↑).

## `rf`

| dataset | clock | task | n_batches | uniform_pair | rfperm | resid | fire_rfperm | fire_resid |
|---|---|---|---:|---:|---:|---:|---:|---:|
| `interstate` | time | RMSE | 24 | 714.8018 | 714.8018 | 698.0895 | 0.00 | 0.14 |
| `nyc_taxi` | time | RMSE | 24 | 2.5859 | 2.5851 | 2.5836 | 0.05 | 0.09 |
| `electricity` | time | acc | 24 | 0.8036 | 0.8039 | 0.8014 | 0.23 | 0.14 |
| `airlines` | time | acc | 24 | 0.6727 | 0.6727 | 0.6727 | 0.00 | 0.00 |
| `bike_hour` | time | RMSE | 24 | 57.0618 | 56.9280 | 57.8255 | 0.05 | 0.18 |
| `beijing_pm25` | time | RMSE | 24 | 78.5127 | 82.6325 | 81.1645 | 0.14 | 0.27 |
| `occupancy` | time | acc | 24 | 0.8589 | 0.8307 | 0.8586 | 0.27 | 0.27 |
| `diabetes_readmit` | shift | acc | 24 | 0.6382 | 0.6382 | 0.6382 | 0.00 | 0.00 |
| `california` | spatial | RMSE | 24 | 0.6929 | 0.6982 | 0.6973 | 0.05 | 0.09 |

## `xgb`

| dataset | clock | task | n_batches | uniform_pair | rfperm | resid | fire_rfperm | fire_resid |
|---|---|---|---:|---:|---:|---:|---:|---:|
| `interstate` | time | RMSE | 24 | 683.0381 | 683.0381 | 677.2343 | 0.00 | 0.14 |
| `nyc_taxi` | time | RMSE | 24 | 2.7259 | 2.7158 | 2.8246 | 0.05 | 0.09 |
| `electricity` | time | acc | 24 | 0.8100 | 0.8130 | 0.8068 | 0.23 | 0.18 |
| `airlines` | time | acc | 24 | 0.6482 | 0.6482 | 0.6482 | 0.00 | 0.00 |
| `bike_hour` | time | RMSE | 24 | 49.2405 | 49.6221 | 50.7515 | 0.05 | 0.18 |
| `beijing_pm25` | time | RMSE | 24 | 76.6374 | 79.3827 | 82.7780 | 0.14 | 0.27 |
| `occupancy` | time | acc | 24 | 0.8432 | 0.8457 | 0.8509 | 0.27 | 0.23 |
| `diabetes_readmit` | shift | acc | 24 | 0.6273 | 0.6273 | 0.6273 | 0.00 | 0.00 |
| `california` | spatial | RMSE | 24 | 0.6500 | 0.6494 | 0.6486 | 0.05 | 0.05 |

Quiet streams should match uniform. A real P(Y|X) hop should fire
rfperm; RMSE/Acc then shows whether reweighting helped the
**next** batch. Subset localization is off this board — too thin
at batch size 200.
