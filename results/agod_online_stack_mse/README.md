# Online stacking + erank/modality balance — MSE (MSR-VTT)

## Setup

1. **Online stacking**: modality logits fused by learned `w=softmax(ψ)`
2. **Erank / modality balance** (`stack_erank_bal` only):
   - `L_erank = ReLU(erank_floor − erank(G_h))`
   - `L_align = mean max(cos_ij, 0)²`
   - `L_bal = KL(w ‖ prior)` (prior=α if erank collapsed else uniform)
3. Task = CE; **MSE = Brier** on holdout (`pre→post` drop = good)

## Smoke (honest)

| variant | MSE pre | MSE post | MSE drop↑ | frac drop>0 | Acc↑ | mean erank aux |
|---|---:|---:|---:|---:|---:|---:|
| `mean_ce` | 0.2477 | 0.2576 | -0.0099 | 0.50 | +0.010 | nan |
| `stack_ce` | 0.2618 | 0.2648 | -0.0030 | 0.50 | -0.004 | nan |
| `stack_erank_bal` | 0.2607 | 0.2718 | -0.0111 | 0.50 | -0.012 | 2.99 |

Best MSE-drop: **`stack_ce`** (-0.0030).

On low-corr MSR-VTT, hidden erank stays ≈|M| so `L_erank` barely fires;
mean MSE drop is **not** positive overall — same regime lesson as soft_decorr.
Amazon high-corr is the next place to expect MSE↓.

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stack_mse.py
```
