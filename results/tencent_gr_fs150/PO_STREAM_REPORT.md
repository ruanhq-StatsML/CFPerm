# TencentGR 150-feat PO-risk stream

Clock = last prefix event time (`_seq_t_end`). Label = `future_cnv` (Acc ↑).
Batches = 24 × 200 (used 4800 / 6000 users, 150 feats).
Gate γ=1.5.

| method | Acc | fire |
|---|---:|---:|
| uniform last-two | 0.8298 | — |
| DRE last-two | 0.8223 | — |
| rfperm √PO | 0.8298 | 0.00 |
| resid drop-old | 0.8298 | 0.00 |

Quiet hops should match uniform. Fire only when consecutive OOS probe
error jumps with a non-vacuous denominator.
