# MSR-VTT AGOD portable prototype

Same control plane as Amazon (`agod/`): **sensor → EMA state → param-group LR actuator → optional hard-gate adaptor**.

## Setup

- Data: `data/msrvtt/packed/{video,text,audio}_feat.npy` + `labelsmsr.npy` (N=5000, binary)
- Stream: 1 ref + 6 windows with **induced P(Y) drift** (`p1` from 0.25 → 0.75)
- FLOPs units are **projection-head only** (features pre-extracted); adaptor savings are visible

## Smoke summary

| Policy | Acc lift | Acc post | flops_rel | cost_utility |
|---|---:|---:|---:|---:|
| B1 equal | +0.070 | 0.671 | 1.000 | +0.070 |
| B2 RF con−cov | +0.064 | 0.667 | 1.000 | +0.064 |
| B5 MMD-cov+PO | +0.063 | 0.675 | 1.000 | +0.063 |
| B5g B5+gate | +0.048 | 0.620 | **0.746** | +0.053 |

- B5−B1 Acc lift ≈ **−0.007** (≈ on par; not Amazon-style +0.03 on this smoke)
- B5g saves **~25%** proj adapt FLOPs (`flops_rel=0.746`) with Acc lift trade-off
- Under B5, text often gets highest α / LR as label ratio drifts

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_msrvtt_mmd_lr.py
```

Outputs: `results/agod_msrvtt/agod_msrvtt_mmd_lr.json`, `AGOD_MSRVTT_MMD_LR_Acc_Dashboard.png`
