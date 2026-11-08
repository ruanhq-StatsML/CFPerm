# FICO subset localization insights (auto)

Package: `feature_localize` → `subset_on_driver` → `subset_localize` (SMD + PSI / TVD / Hellinger).

Settings: max_n_per_arm=`1500`, B_boot=`30`, n_perm_pfi=`25`, weak-cd draws=`5`, intensity=`2.5`.

## Comparison snapshot

| design | mean |contrast| | mean ΔMSE low | mean ΔMSE high | mean hit@k |
|---|---:|---:|---:|---:|
| median_split_covshift | 0.0139 | 0.0237 | 0.1673 | 0.20 |
| weak5_nonlinear_cd | 0.0150 | 0.1462 | 0.1360 | 0.36 |

## Insights

### median_split_covshift

- Driver-cut concentration mostly on **`low`** (1/3 cells).
- |contrast| mean=`0.0139` (max=`0.0277`) — larger ⇒ stronger low/high asymmetry.
- Mean ΔMSE: low=`0.0237`, high=`0.1673` (high/low ratio≈`7.06`).
- PFI hit@k mean=`0.20` (AUC mean=`0.92`).
- Top Hellinger localizations (z-cut × feature):
  - `ExternalRiskEstimate` grp=4 / `NumInqLast6Mexcl7days`: H=1.000, PSI=54.609, SMD=-1.528 (n0=322, n1=125)
  - `ExternalRiskEstimate` grp=4 / `NumInqLast6M`: H=1.000, PSI=54.570, SMD=-1.787 (n0=300, n1=108)
  - `MSinceOldestTradeOpen` grp=4 / `NumInqLast6Mexcl7days`: H=1.000, PSI=54.281, SMD=-2.175 (n0=301, n1=200)
  - `ExternalRiskEstimate` grp=3 / `NumInqLast6M`: H=1.000, PSI=54.246, SMD=-2.067 (n0=315, n1=208)
  - `ExternalRiskEstimate` grp=3 / `NumInqLast6Mexcl7days`: H=1.000, PSI=54.229, SMD=-1.940 (n0=309, n1=232)
- Largest |SMD| pockets:
  - `MSinceMostRecentTradeOpen` grp=3 → `NetFractionRevolvingBurden` SMD=35.781
  - `MSinceOldestTradeOpen` grp=0 → `NetFractionRevolvingBurden` SMD=35.151
  - `MSinceOldestTradeOpen` grp=3 → `NetFractionRevolvingBurden` SMD=33.242
- True drivers appearing in PFI top-5: `3` row-hits across draws (unique features: ['ExternalRiskEstimate', 'NumInqLast6M', 'NumInqLast6Mexcl7days']).

### weak5_nonlinear_cd

- Driver-cut concentration mostly on **`high`** (5/5 cells).
- |contrast| mean=`0.0150` (max=`0.0300`) — larger ⇒ stronger low/high asymmetry.
- Mean ΔMSE: low=`0.1462`, high=`0.1360` (high/low ratio≈`0.93`).
- PFI hit@k mean=`0.36` (AUC mean=`0.63`).
- Top Hellinger localizations (z-cut × feature):
  - `MSinceMostRecentDelq` grp=0 / `ExternalRiskEstimate`: H=0.594, PSI=17.001, SMD=3.258 (n0=26, n1=35)
  - `MSinceMostRecentDelq` grp=0 / `AverageMInFile`: H=0.573, PSI=14.273, SMD=3.800 (n0=26, n1=35)
  - `MSinceMostRecentDelq` grp=0 / `NetFractionRevolvingBurden`: H=0.560, PSI=13.810, SMD=-4.492 (n0=26, n1=35)
  - `MSinceMostRecentDelq` grp=0 / `MSinceMostRecentInqexcl7days`: H=0.357, PSI=6.137, SMD=1.660 (n0=26, n1=35)
  - `MSinceMostRecentInqexcl7days` grp=0 / `NetFractionRevolvingBurden`: H=0.340, PSI=2.545, SMD=-6.286 (n0=65, n1=66)
- Largest |SMD| pockets:
  - `NumBank2NatlTradesWHighUtilization` grp=0 → `NetFractionRevolvingBurden` SMD=11.588
  - `MSinceMostRecentInqexcl7days` grp=0 → `NetFractionRevolvingBurden` SMD=-6.286
  - `MSinceMostRecentInqexcl7days` grp=0 → `AverageMInFile` SMD=6.266
- True drivers appearing in PFI top-5: `9` row-hits across draws (unique features: ['AverageMInFile', 'ExternalRiskEstimate', 'MSinceMostRecentInqexcl7days', 'MaxDelq2PublicRecLast12M', 'NetFractionRevolvingBurden', 'NumInqLast6M', 'NumSatisfactoryTrades', 'PercentTradesNeverDelq']).

## Cross-design justification

- **Median-split** localizes a *single strong axis* of covariate shift; ΔVIMP/ΔMSE often sit on one arm of that axis, and Hellinger peaks along the split feature's own bins.
- **Weak-5 nonlinear CD** keeps `P(X)` shared, so localization is about *where in X-space the concept map bites*: `subset_on_driver` finds the PFI-top cut with larger |ΔVIMP| / ΔMSE, and z-cuts with high Hellinger/PSI mark the pockets where selected features diverge across batches under the new labeling.
- Operational takeaway: after CFPerm flags a shift, run this package to **name the subset (low/high of j★)** that concentrates the damage — that is the actionable post-hoc localization, not a CS/CD decomposition.

Tables: `*_summary.csv`, `*_subset_driver.csv`, `*_zcut_dist.csv`, `subset_driver_comparison.csv`.

