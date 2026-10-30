# PO-risk Var(u) estimators

PO-risk Var(u) estimators on PermuCATE φ: Causal Forest leaf-size-debiased (+ leaf merge), BCF posterior, bootstrap fallback.

**Recommend:** Default: Causal Forest leaf-size-debiased + leaf-wise merge (cf_Var_u_tree / cf_V_debiased). BCF when Bayesian pipeline exists. Bootstrap only as fallback.

BCF backend: `python_bagged_bcf`

## Separation (malignant vs benign)

| estimator | benign | malignant | gap | ratio | OK? |
| --- | --- | --- | --- | --- | --- |
| subenv_s | 0.8724 | 4.896 | 4.023 | 5.61 | yes |
| emp_var_u | 62.45 | 1480 | 1418 | 23.7 | yes |
| cf_V_naive | 0.002822 | 16.88 | 16.87 | 5.98e+03 | yes |
| cf_V_debiased | 0 | 16.37 | 16.37 | inf | yes |
| cf_Var_u_tree | 3.992 | 8.101e+06 | 8.101e+06 | 2.03e+06 | yes |
| cf_Var_u_delta | 0 | 4.521e+04 | 4.521e+04 | inf | yes |
| bcf_Var_u | 32 | 938.1 | 906.1 | 29.3 | yes |
| boot_Var_u | 24.42 | 364.7 | 340.3 | 14.9 | yes |

hit-rate: **1.00**

```bash
python3 scripts/run_as_dta_po_risk_var.py --quick
```

LaTeX: `docs/method/AS_DTA_po_risk_var_estimators.tex`
