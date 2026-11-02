# Adaptation opportunity probe

O_local = (D L R M)^(1/4) gates sep/expand this batch; O_persist = O_local·Stab discounts carrying S/offset into unknown next batch.

| stage | O_local | O_persist | D | L | R | M | Stab | action | gap_lift | expand |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.773 | 0.270 | 1.00 | 0.41 | 0.86 | 1.00 | 0.35 | sep+expand | 0.356 | True |
| 1 | 0.668 | 0.434 | 0.55 | 0.42 | 0.86 | 1.00 | 0.65 | sep+expand | 0.441 | True |
| 2 | 0.707 | 0.542 | 0.69 | 0.42 | 0.86 | 1.00 | 0.77 | sep+expand | 0.477 | True |

control: structured O_local=**0.773** vs flat **0.001** (lift **0.772**)

Track benches: C-GCD (CIFAR-100 / IN-100 / Tiny / CUB); PDA/OSDA (Office-Home / VisDA / DomainNet).

```bash
python3 scripts/run_adaptation_opportunity_probe.py --quick
```
