# 客群 / region Mag×Var (PermuCATE PO)

客群/region Mag×Var on PermuCATE unit PO — triage scoped per segment (finance / medicine / topic), not one global alarm.

## Estimator

- Unit $u=(\varphi-\hat\tau)^2$; $Y$=response; $W$=REF/LIVE
- Severity $=J + Y$-flip; triage **per segment**

## Table

| segment | region | sev(benign) | sev(malignant) | gap | separates? |
| --- | --- | --- | --- | --- | --- |
| general | hc3_all | 0.039 | 0.413 | 0.375 | yes |
| finance | hc3_finance | 0.087 | 0.448 | 0.361 | yes |
| medicine | hc3_medicine | 0.012 | 0.479 | 0.467 | yes |
| sci | news20_sci | -0.006 | 0.458 | 0.464 | yes |

segment hit-rate: **1.00**

Same Mag×Var geometry, applied per 客群/region. Finance and medicine can need different actuators even under one global detector — segment-local J / severity is the product table.

```bash
python3 scripts/run_as_dta_segment_magvar.py --quick
```
