# RSI Iteration Log

| tick | UTC | idea | justify (1 line) | change | tests |
|---|---|---|---|---|---|
| 1 | 2026-09-24 | label-perm null + excess_auc + probe_eff on adjacent board | raw AUC unidentifiable as skill; excess vs chance / FLOPs proxy is cross-pack comparable | `agod/transfer_null.py` + board enrich | `test_transfer_null` + board helpers |
| 2 | 2026-09-24 | equal-width ECE on transfer HGB/LogReg probs | ranking≠calibration; high excess+high ECE means transferable order but untrustworthy probs | `expected_calibration_error` + board fields | `test_ece_*` + board |
| 3 | 2026-09-24 | blocked bootstrap CI90 on mean excess / ECE | adjacent pairs share endpoints → iid bootstrap understates var; block=2 restores conservative CI | `blocked_bootstrap_ci` + pack summary | `test_blocked_bootstrap_ci_*` |
| 4 | 2026-09-24 | reservoir pair subsample under max_pairs | fixed-k unbiased SSRS for pack means; linspace≠simple random | `select_pair_indices` default reservoir | `test_reservoir_sample_*` |
