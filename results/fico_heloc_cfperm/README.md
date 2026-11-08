# FICO HELOC CFPerm — VIMP power + subset localization

- Clean n = `9861` (raw `10459`).
- Split features (top-3 LR |coef|): `NumInqLast6M`, `NumInqLast6Mexcl7days`, `ExternalRiskEstimate`.
- Settings: `B=150` trees, `n_perm=200`, `n_perm_cate=50`, `max_n_per_arm=1500`.

## Power

| scope | reject_rate | recover_rate | top1_hit | top3_hit | n_cells |
|---|---:|---:|---:|---:|---:|
| overall | 0.222 | 0.111 | — | — | 9 |
| loco | 0.000 | 0.000 | 0.000 | 0.333 | 3 |
| permucate | 0.333 | 0.000 | 0.000 | 0.000 | 3 |
| grf | 0.333 | 0.333 | 0.333 | 0.333 | 3 |

## Cells

| split | vimp | reject | recovered | top features |
|---|---|---:|---:|---|
| `NumInqLast6M` | loco | 0 | 0 | `MSinceOldestTradeOpen`, `PercentInstallTrades`, `MSinceMostRecentInqexcl7days` |
| `NumInqLast6M` | permucate | 0 | 0 | `MSinceMostRecentInqexcl7days`, `MSinceOldestTradeOpen`, `PercentInstallTrades` |
| `NumInqLast6M` | grf | 0 | 0 | `MSinceMostRecentInqexcl7days`, `ExternalRiskEstimate`, `AverageMInFile` |
| `NumInqLast6Mexcl7days` | loco | 0 | 0 | `PercentTradesWBalance`, `NumInqLast6Mexcl7days`, `ExternalRiskEstimate` |
| `NumInqLast6Mexcl7days` | permucate | 0 | 0 | `MSinceMostRecentInqexcl7days`, `ExternalRiskEstimate`, `AverageMInFile` |
| `NumInqLast6Mexcl7days` | grf | 0 | 0 | `MSinceMostRecentInqexcl7days`, `NetFractionRevolvingBurden`, `MSinceOldestTradeOpen` |
| `ExternalRiskEstimate` | loco | 0 | 0 | `NumTrades90Ever2DerogPubRec`, `NumTrades60Ever2DerogPubRec`, `NetFractionRevolvingBurden` |
| `ExternalRiskEstimate` | permucate | 1 | 0 | `NumSatisfactoryTrades`, `PercentTradesWBalance`, `NetFractionRevolvingBurden` |
| `ExternalRiskEstimate` | grf | 1 | 1 | `ExternalRiskEstimate`, `NetFractionRevolvingBurden`, `MSinceOldestTradeOpen` |

## Subset localization (CFPerm bins)

Quantile bins of rejected / top VIMP features; `delta_y = mean(Y|W=1) - mean(Y|W=0)` localizes where the batch gap concentrates. See `subset_localization.csv`.

## Post-hoc localization (PFI degree → driver cut)

See [`POSTHOC.md`](POSTHOC.md). Pipeline: Ridge-PFI degree ranks drivers → `subset_on_driver` median cut with bootstrap ΔVIMP CI + ΔMSE → z-cut SMD with PSI / TVD / Hellinger.

| split | driver | concentrated | contrast | ΔMSE low | ΔMSE high |
|---|---|---|---:|---:|---:|
| `NumInqLast6M` | `PercentInstallTrades` | low | -0.0001 | 0.0280 | 0.0182 |
| `NumInqLast6Mexcl7days` | `NumInqLast6M` | — | — | — | 0.4259 |
| `ExternalRiskEstimate` | `NumInqLast6M` | high | 0.0277 | 0.0194 | 0.0579 |

Tables: `posthoc_feature_rank.csv`, `posthoc_subset_driver.csv`, `posthoc_zcut_dist.csv`, `posthoc_summary.csv`.

### Takeaway

Under B=150 / n_perm=200, **GRF** is the only VIMP that both rejects and recovers the continuous split `ExternalRiskEstimate`. Inquiry median-splits are weak / discrete alternatives (often no reject). Post-hoc localization still surfaces where ΔMSE / ΔVIMP concentrates (high arm on `ExternalRiskEstimate` split).

