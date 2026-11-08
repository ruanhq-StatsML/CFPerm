# Post-hoc subset localization (FICO HELOC)

After CFPerm / median-split batches: PFI degree → driver `j★`, then `subset_on_driver` (low/high median cut) with bootstrap ΔVIMP CI + ΔMSE; z-cuts report SMD + PSI / TVD / Hellinger.

## Summary

| split | driver | concentrated | contrast | ΔMSE low | ΔMSE high | hit@k |
|---|---|---|---:|---:|---:|---:|
| `NumInqLast6M` | `PercentInstallTrades` | low | -0.0001 | 0.0280 | 0.0182 | 0.20 |
| `NumInqLast6Mexcl7days` | `NumInqLast6M` |  | nan | nan | 0.4259 | 0.20 |
| `ExternalRiskEstimate` | `NumInqLast6M` | high | 0.0277 | 0.0194 | 0.0579 | 0.20 |

## Files

- `posthoc_feature_rank.csv` — |ΔPFI| ranking per split
- `posthoc_subset_driver.csv` — low/high ΔVIMP (±1.96 SE) and ΔMSE
- `posthoc_zcut_dist.csv` — group-wise SMD + PSI/TVD/Hellinger
- `posthoc_summary.csv` / `.json`

