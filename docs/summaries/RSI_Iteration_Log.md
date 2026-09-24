# RSI Iteration Log

| tick | UTC | idea | justify (1 line) | change | tests |
|---|---|---|---|---|---|
| 1 | 2026-09-24 | label-perm null + excess_auc + probe_eff on adjacent board | raw AUC unidentifiable as skill; excess vs chance / FLOPs proxy is cross-pack comparable | `agod/transfer_null.py` + board enrich | `test_transfer_null` + board helpers |
| 2 | 2026-09-24 | equal-width ECE on transfer HGB/LogReg probs | ranking≠calibration; high excess+high ECE means transferable order but untrustworthy probs | `expected_calibration_error` + board fields | `test_ece_*` + board |
| 3 | 2026-09-24 | blocked bootstrap CI90 on mean excess / ECE | adjacent pairs share endpoints → iid bootstrap understates var; block=2 restores conservative CI | `blocked_bootstrap_ci` + pack summary | `test_blocked_bootstrap_ci_*` |
| 4 | 2026-09-24 | reservoir pair subsample under max_pairs | fixed-k unbiased SSRS for pack means; linspace≠simple random | `select_pair_indices` default reservoir | `test_reservoir_sample_*` |
| 5 | 2026-09-24 | **pivot** PO ref/probe/refit efficiency scorecard | ranking≠MSE; probe always-on FLOPs vs reject-only refit → rank_eff/mse_eff | `agod/po_eff.py` + scorecard script | `test_po_eff` |
| 6 | 2026-09-24 | gate duty → E[adapt FLOPs] + rank_eff_E | planning needs E[cost]=duty·n·fit; probe duty-invariant; budget_ratio≈duty | `expected_adapt_flops` | `test_expected_flops_*` |
| 7 | 2026-09-24 | √ vs ∛ IPTW power softness scorecard + **3min cadence** | same FLOPs ⇒ α is bias–variance only; ∛≤√ often but uniform still wins | `agod/po_power_eff.py` + tomorrow demo board | `test_po_power_eff` |
| 8 | 2026-09-24 | tomorrow demo PNG + close P3 image-OOD | demo needs one figure; PO image-OOD is negative — don't spend FLOPs there | `plot_rsi_tomorrow_demo.py` | `test_plot_rsi_tomorrow_demo` |
| 9 | 2026-09-24 | Grad-RFPerm freeze MSE–FLOPs scorecard | freeze_eff=(1−mse_rel)/flops_rel; Pareto when both &lt;1 vs always_adapt | `agod/freeze_eff.py` | `test_freeze_eff` |
| 10 | 2026-09-24 | 3-panel tomorrow demo (rank + soft + freeze Pareto) | one figure for review; green quadrant = dominate always_adapt | `plot_rsi_tomorrow_demo.py` | plot smoke |
| 11 | 2026-09-24 | duty breakeven + TOMORROW_INDEX | prefer refit on budget iff duty &lt; 1/n_control; one-click review index | `duty_breakeven_*` + index md | `test_po_eff` |
| 12 | 2026-09-24 | n_control sensitivity in scorecard MD | larger control windows erase refit budget edge | `PO_EFF_SCORECARD` §n_control | regen |
| 13 | 2026-09-24 | walkthrough note + refresh demo artifact | 12 tests green; one figure for human review tomorrow | `RSI_Walkthrough.md` + artifacts png | 12 passed |
| 14 | 2026-09-24 | Δmse_eff(∛−√) + strict freeze JSON | soft gap quantifies how much ∛ saves vs √ per FLOP; NaN→null | `po_power_eff` + freeze dump | power+freeze tests |
