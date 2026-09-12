# MSR-VTT AGOD portable prototype

Same control plane as Amazon (`agod/`): sensor -> EMA state -> param-group LR actuator -> optional hard-gate adaptor.

| Policy | Acc lift | flops_rel | cost_utility |
|---|---:|---:|---:|
| B1 | +0.000 | 1.000 | +0.000 |
| B2 | +0.000 | 1.000 | +0.000 |
| B5 | +0.000 | 1.000 | +0.000 |
| B5g | +0.000 | 0.714 | +0.000 |
| B5r | +0.000 | 0.841 | +0.000 |

- B5-B1 Acc lift = **+0.000**
- B5g flops_rel = **0.714** (vs B1 1.000)

```bash
PYTHONPATH=. python3 scripts/run_agod_msrvtt_mmd_lr.py
```
