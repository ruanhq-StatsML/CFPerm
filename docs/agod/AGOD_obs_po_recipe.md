# AGOD obs-PO recipe (CFPerm gate, v8.1)

## Thesis

认 **hard**; packMSE only under **beijing**; **L0 = CFPerm DRPerm** (not OnlineRFPerm).

## Locked

| Item | Choice |
|---|---|
| L0 | CFPerm DRPerm, `e_mode=known`, `n_perm≥39`, α=0.05 |
| L1 | intensity → λ / beijing (gate 0.20 / beijing 0.45) |
| Hard claim | `hard_m` = hard_support on reject (mild) |
| RFPerm-era dual | mild→hard, BJ→soft CV, `n_recent=1`; blend mix≤0.5 |
| Reject | `dual_r2` (n_recent=2); blend mix=0.75 |

## Open under CFPerm reject sets (v8.1)

- Multi-seed Jaccard(CFPerm, RFPerm) ≈ **0.0–0.22** → gates disagree
- Sparse-duty stocks: hard_m pack −3%…−4%
- High-duty metro/beijing: hard reweight hurts pack; beijing soft path rare (beijing_frac≲0.07)
- Next: full 40×256 multi-seed; retune beijing intensity threshold / temper; dual+blend_50 hybrid only if BJ fires

## Eval protocol

1. Synth size/power for L0
2. Stream duty + Jaccard vs RFPerm
3. Sig-only hard-subset next-MSE (primary hard claim)
4. Sig-only pack MSE under beijing (conditional)
