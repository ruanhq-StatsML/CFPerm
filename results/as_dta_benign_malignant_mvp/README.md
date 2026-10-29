# Benign vs Malignant OOD MVP (PO mean + variance)

## Justify

Magnitude (ΔMean PO) first; if quiet → hold. If Mean↑, check Variance/Peak/CVaR/top-share: flat → benign OOD → conformal; spiked → malignant → freeze-lexical.

- **Same:** Both can raise mean PO / OOD score
- **Diff:** Benign = uniform; Malignant = heavy-tail / high subenv variance

## Separation panel

| batch | mean PO | subStd | peak | CVaR | verdict | route |
| --- | --- | --- | --- | --- | --- | --- |
| ref | 0.1332 | 0.0067 | 0.0147 | 0.1499 | quiet | hold |
| benign_ood | 0.2110 | 0.0176 | 0.0292 | 0.2433 | benign_ood | conformal_calibrate |
| malignant_ood | 0.3899 | 0.0255 | 0.0551 | 0.4416 | malignant_ood | freeze_lexical_adapt_syntactic |

## Anti-patterns

Do **not** hard-filter high-PO rows or pseudo-label flips from PO.
See dashboard right panel for selection-bias / label-noise numbers.

![](benign_vs_malignant_ood.png)

```bash
python3 scripts/run_as_dta_benign_malignant_mvp.py
```
