# Streaming performance curves — weights vs no weights

Seeds [2026, 2027, 2028], quick=False.

## Novel / known / forgetting by stage

<img alt="stage curves" src="/workspace/results/repr_sep_cold_start/figures/streaming_stage_curves.png" width="900" />

## One-pager

<img alt="combo" src="/workspace/results/repr_sep_cold_start/figures/streaming_performance_curves.png" width="900" />

## Efficiency (budget → novel_acc)

<img alt="efficiency" src="/workspace/results/repr_sep_cold_start/figures/streaming_efficiency_curves.png" width="560" />

## Δ novel vs no weights

<img alt="delta" src="/workspace/results/repr_sep_cold_start/figures/streaming_delta_vs_noweight.png" width="560" />

```bash
python3 scripts/plot_streaming_performance.py
```
