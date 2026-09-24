# RSI Iteration Log

| tick | UTC | idea | justify (1 line) | change | tests |
|---|---|---|---|---|---|
| 1 | 2026-09-24 | label-perm null + excess_auc + probe_eff on adjacent board | raw AUC unidentifiable as skill; excess vs chance / FLOPs proxy is cross-pack comparable | `agod/transfer_null.py` + board enrich | `test_transfer_null` + board helpers |
| 2 | 2026-09-24 | equal-width ECE on transfer HGB/LogReg probs | ranking≠calibration; high excess+high ECE means transferable order but untrustworthy probs | `expected_calibration_error` + board fields | `test_ece_*` + board |
