# PO-risk Var(u) — synthetic + real-world justify

Product landing justify: CF-debiased Var(u) + Mag×Var triage on synthetic DGPs and real-world text boards.

**Recommend:** Ship CF leaf-size-debiased + leaf merge as default Var(u); gate with Mag×Var + Y-flip; BCF optional; bootstrap fallback.

- Synthetic hit-rate: **1.00** (5/5)
- Real-world pair hit-rate: **1.00** (4/4)
- Meets expectation: **True**

## Product logic (next step)

1. fit φ, τ̂
2. u=(φ−τ̂)²
3. Mag×Var (μ, s)
4. Var(u) CF-debiased (default)
5. Y-flip if text integrity
6. actuator via triage()
7. honesty: bias_ratio

## Synthetic

| board | oracle | pred | V_debiased | bias_ratio | OK? |
| --- | --- | --- | --- | --- | --- |
| quiet | quiet | quiet | 0 | 1.00 | yes |
| benign_covariate | benign | benign | 0 | 1.00 | yes |
| malignant_sparse | malignant | malignant | 0 | 1.00 | yes |
| malignant_dense | malignant | malignant | 0 | 1.00 | yes |
| label_delay | malignant | malignant | 0 | 1.00 | yes |

## Real-world

| board | ben pred | mal pred | V_gap | pair OK? |
| --- | --- | --- | --- | --- |
| hc3_finance | quiet | malignant | 0 | yes |
| news20_sci | quiet | malignant | 0 | yes |
| ag_news | quiet | malignant | 0 | yes |
| imdb | quiet | malignant | 0 | yes |

```bash
pip install -e Python/po_risk_var
po-risk-var justify --quick
# or: python3 -m po_risk_var.cli justify --quick
```
