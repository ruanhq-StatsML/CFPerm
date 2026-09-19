# LOCO PO-risk tip timing v2 sweep

shift_at=20, seeds=[0, 1, 2, 3, 4]

| mode | true delay mean±std | hit±1 | miss | early/run | wrong hit±1 | wrong miss |
|---|---:|---:|---:|---:|---:|---:|
| `sum` | 0.00±0.00 | 100% | 0% | 0.00 | 100% | 0% |
| `abs_dev` | -2.20±4.40 | 80% | 0% | 0.40 | 40% | 0% |
| `po` | 0.00±0.00 | 100% | 0% | 0.00 | 100% | 0% |
| `collapse` | -0.80±1.60 | 80% | 0% | 0.20 | 100% | 0% |
| `cusum` | 0.00±0.00 | 100% | 0% | 0.00 | 100% | 0% |
| `confirm` | 0.00±0.00 | 100% | 0% | 0.00 | 0% | 100% |

**Preferred default:** `confirm` (true hit±1=100%, wrong miss=100%, early/run=0.00).

confirm seed0 tip t*=20 delay=0.
