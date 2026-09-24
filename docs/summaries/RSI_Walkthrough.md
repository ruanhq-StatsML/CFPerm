# RSI Walkthrough — what to open tomorrow

Evidence figure (3 panels):

<img alt="RSI PO efficiency demo" src="/opt/cursor/artifacts/rsi_tomorrow_demo.png" />

Also in-repo: `results/agod_po_eff/rsi_tomorrow_demo.png`

## 60-second script
1. Open [`TOMORROW_INDEX.md`](../../results/agod_po_eff/TOMORROW_INDEX.md)
2. Read panel **A**: refit bars taller → rank per FLOP
3. Panel **B**: uniform dominates IPTW; ∛ softer when you must gate
   - **median** Δmse_eff(∛−√)≈+2e-4 (mean is metro-skewed — ignore mean)
4. Panel **C**: green quadrant = freeze beats always_adapt on MSE *and* FLOPs (electricity)

## Reproducibility
```bash
PYTHONPATH=. python3 -m pytest tests/test_po_eff.py tests/test_po_power_eff.py tests/test_freeze_eff.py tests/test_plot_rsi_tomorrow_demo.py -q
# 12 passed
```

Branch: `cursor/rsi-model-stats-eff-abce`
