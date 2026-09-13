# AGOD adapter layers + dynamic L3 gate

## How many layers / how to characterize

| Layer | Role | Adjustable knobs |
|---|---|---|
| L0 sensor | cov/concept scores | MMD / PO / MSG |
| L1 state | EMA Softmax α | ema, τ |
| L2 soft LR | continuous next-stage step sizes | β, lr0 |
| L3 hard gate | sparse adapt (proj BWD only) | θ / quantile / hyst / rand |

FWD inference always paid. Efficiency claim is **adapt FLOPs**, not latency.

## Dynamic L3 policies

- `none` — soft LR only
- `fixed` — α ≥ θ
- `quantile` — drop bottom-α mass
- `ema_theta` — θ tracks EMA(mean α)
- `hysteresis` — on/off thresholds (less chatter)
- `random` — non-attribution sparsity control

## Smoke board (Amazon + MSR-VTT)

| Dataset | L3 gate | Acc lift | flops_rel | util | switches/win |
|---|---|---:|---:|---:|---:|
| amazon | none | +0.114 | 1.000 | +0.114 | 0.00 |
| amazon | fixed | +0.114 | 0.778 | +0.158 | 0.20 |
| amazon | quantile | +0.114 | 0.733 | +0.156 | 0.00 |
| amazon | ema_theta | +0.114 | 0.778 | +0.158 | 0.20 |
| amazon | hysteresis | +0.114 | 0.733 | +0.156 | 0.00 |
| amazon | random | +0.109 | 0.911 | +0.127 | 0.60 |
| msrvtt | none | +0.016 | 1.000 | +0.016 | 0.00 |
| msrvtt | fixed | +0.019 | 0.841 | +0.020 | 1.40 |
| msrvtt | quantile | +0.024 | 0.810 | +0.029 | 1.20 |
| msrvtt | ema_theta | +0.009 | 0.873 | +0.007 | 1.00 |
| msrvtt | hysteresis | +0.009 | 0.841 | +0.005 | 1.00 |
| msrvtt | random | +0.045 | 0.841 | +0.057 | 1.20 |

### Readout

- **Amazon:** attribution gates (`fixed/quantile/hysteresis/ema_theta`) keep Acc lift ≈ `none` (+0.114) while cutting adapt FLOPs to ~0.73–0.78; `random` is weaker on util and chatters more.
- **MSR-VTT:** `quantile` is the best attributed gate on Acc/util; `hysteresis/ema_theta` are stabler (fewer switches) but milder lift; `random` can look strong on Acc in short smoke but is not attribution-guided.

```bash
PYTHONPATH=. python3 scripts/run_agod_adapter_gate_layers.py
```
