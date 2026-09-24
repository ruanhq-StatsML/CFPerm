# Tomorrow review index — RSI PO/AGOD efficiency

Branch: `cursor/rsi-model-stats-eff-abce` · cadence ~3 min

## Open this first
1. **Figure:** [`rsi_tomorrow_demo.png`](./rsi_tomorrow_demo.png) (3 panels)
2. **Narrative:** [`../../docs/summaries/RSI_Tomorrow_Demo_Board.md`](../../docs/summaries/RSI_Tomorrow_Demo_Board.md)
3. **Log:** [`../../docs/summaries/RSI_Iteration_Log.md`](../../docs/summaries/RSI_Iteration_Log.md)

## Claims (quiz yourself)
| # | Claim | Where |
|---|---|---|
| 1 | refit wins rank_eff 6/6; mse_eff often &lt;0 | `PO_EFF_SCORECARD.md` |
| 2 | duty breakeven = 1/n_control; below ⇒ refit cheaper than probe | scorecard headline |
| 3 | ∛≤√ on ~67% packs; uniform still best on 5/6 | `../agod_po_power_eff/PO_POWER_EFF.md` |
| 4 | electricity freeze Pareto (~0.86× MSE @ 0.7× FLOPs) | `../agod_freeze_eff/FREEZE_EFF.md` |
| 5 | Do not spend PO FLOPs on image-OOD | `docs/agod/AGOD_image_ood_bench.md` |

## Rebuild in 30s
```bash
PYTHONPATH=. python3 scripts/run_po_eff_scorecard.py
PYTHONPATH=. python3 scripts/run_po_power_eff_scorecard.py
PYTHONPATH=. python3 scripts/run_freeze_eff_scorecard.py
PYTHONPATH=. python3 scripts/plot_rsi_tomorrow_demo.py
```
