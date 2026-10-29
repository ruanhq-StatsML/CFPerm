# OOD vs intentional injection — with labels Y

With Y: PO_i=1-p̂(Y_i|x_i). Uniform ΔMean + flat concentration → natural OOD; rising top-20% PO share / CVaR + error/FN → intentional injection. Effectiveness = Δerror/ΔFN/Δrecall/Δcost after actuator.

## Verdicts

| case | truth | verdict | match | ΔMean | ΔCVaR | ΔTop20 | error | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ref_clean | quiet | quiet | True | 0.0000 | -0.0033 | 0.0000 | 0.143 | 0.018 |
| ood_uniform | ood_uniform | ood_uniform | True | 0.0043 | -0.0012 | -0.0100 | 0.137 | 0.022 |
| intentional_injection | intentional_injection | intentional_injection | True | 0.0498 | 0.0421 | -0.0478 | 0.201 | 0.131 |
| cross_domain_ood | quiet | quiet | True | -0.0849 | -0.0906 | 0.1010 | 0.077 | 0.024 |

Verdict accuracy: **100%**

## Effectiveness with Y (after actuator)

| case | worked | Δrecall | ΔFN | Δcost% | Δerror |
| --- | --- | --- | --- | --- | --- |
| intentional_injection | True | +0.100 | -100 | -56.0 | -0.120 |
| ood_uniform | True | -0.019 | +19 | -36.4 | -0.042 |
| cross_domain_ood | True | -0.001 | +1 | -32.8 | -0.023 |

## Why Y is informative

- `PO_i = 1 - p̂(Y_i|x_i)` is true-label risk — not just unlabeled OOD energy.
- Intentional injection **concentrates** PO (top-20% share / CVaR) and lifts FN.
- Natural OOD lifts Mean more **uniformly**; concentration stays flatter.
- After the actuator, Y gives a direct scorecard: ΔFN / Δrecall / Δcost.

```bash
python3 scripts/run_as_dta_ood_vs_intentional.py
```
