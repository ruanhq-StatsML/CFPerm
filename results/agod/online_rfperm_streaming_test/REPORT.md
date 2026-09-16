# Streaming OnlineRFPerm test (continuous)

> RF shallow probe = component model on the tidy infer stream.  
> Continuous clock: freeze + rolling ratio, and sliding consecutive OOS.  
> **If a hop is present and RF does not fire → under this probe the chain looks smooth.**

## Setup

| item | value |
|------|-------|
| source streams | `results/agod/online_rfperm_multi_datasets/*/stream.jsonl` |
| faith / Y | answer-precision → binary $Y$ |
| win / gate | 20 / 1.25 |
| ref | quiet window before cut (`n_ref=100` on default multi) |

## Modes

1. **orf_freeze** — fit RF on quiet; stream per-obs error; fire when rolling mean / quiet ≥ γ  
2. **orf_slide** — slide window by 1; OnlineRFPerm consecutive OOS gate  
3. **smooth_control** — rewrite hop rows as quiet resamples; expect **no** fire

Delay = first trail observation index with fire (`0` = first obs after cut).

## Results

| dataset | y quiet→hop | freeze delay | slide delay | smooth first1 (fr/sl) | read |
|---------|-------------|-------------:|------------:|-----------------------|------|
| halueval | 0.00→1.00 | 0 | 1 | None/None | hop_caught_smooth_control_quiet — RF component sees the law change |
| squad | 0.00→1.00 | 0 | 1 | None/None | hop_caught_smooth_control_quiet — RF component sees the law change |
| hotpotqa | 0.00→1.00 | 0 | 1 | None/None | hop_caught_smooth_control_quiet — RF component sees the law change |
| truthfulqa | 0.00→1.00 | 0 | 1 | None/None | hop_caught_smooth_control_quiet — RF component sees the law change |

## How to read

- Hop injected + fire soon → RF component **caught** the quality-law change.  
- Hop injected + **no** fire → inference chain looks **丝滑 / smooth** under this probe (or labels/features hide the hop).  
- Smooth control should stay quiet; if it fires, the gate is too hot.

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_streaming_test.py --win 20 --gate 1.25
```
