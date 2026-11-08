# FICO HELOC — updated result comparison

Settings: cleaned HELOC ($n=9861$), `B=150`, `n_perm=200`, `max_n_per_arm=1500`.
PermuCATE uses RandomForest outcome models.

## CFPerm power

| design | VIMP | reject | driver hit | recall@5 |
|---|---|---:|---:|---:|
| median-split CS | overall | 0.222 | 0.111 recover | — |
| median-split CS | loco | 0.000 | 0.000 | — |
| median-split CS | permucate | 0.333 | 0.000 | — |
| median-split CS | grf | 0.333 | 0.333 | — |
| **weak-5 nonlinear CD** | **overall** | **0.867** | **1.000 @top5** | **0.493** |
| weak-5 nonlinear CD | loco | 0.600 | 1.000 | 0.520 |
| weak-5 nonlinear CD | permucate | 1.000 | 1.000 | 0.480 |
| weak-5 nonlinear CD | grf | 1.000 | 1.000 | 0.480 |

Median-split cells: inquiry splits never reject; only `ExternalRiskEstimate` rejects (permucate+grf), and only GRF recovers the axis.

## Subset localization

| design | \|contrast\| | ΔMSE low | ΔMSE high | PFI hit@k |
|---|---:|---:|---:|---:|
| median-split CS | 0.014 | 0.024 | 0.167 | 0.20 |
| weak-5 nonlinear CD | 0.015 | 0.146 | 0.136 | 0.36 |

- Median-split: Hellinger/PSI near max on the split axis; ΔMSE high/low ≈ 7×.
- Weak-5 CD: shared `P(X)`; all 5 draws concentrate on **high** of PFI-top; more true-driver hits in PFI top-5.

## Justification

HELOC is a **multi-weak-feature** book. A single median axis is a low-power foil (reject 0.22). Weak-5 + highly nonlinear concept drift matches that regime and yields reject ≈ 0.87 with any-driver@top5 = 1.0.

## LaTeX

- `docs/method/FICO_HELOC_CFPerm_supplement.tex` — full note
- `docs/method/FICO_HELOC_CFPerm_tables_only.tex` — tables only
