# AGOD obs-PO recipe (CFPerm gate, v8.2)

## Thesis

认 **hard**; packMSE only under **beijing**; **L0 = CFPerm DRPerm** (not OnlineRFPerm).

## Locked

| Item | Choice |
|---|---|
| L0 | CFPerm DRPerm, `e_mode=known`, `n_perm≥39`, α=0.05 |
| L1 | intensity → λ / beijing; **`beijing_gate=0.25`** (v8.2; was 0.45) |
| temper | gate=0.20, lam_max=0.75 |
| Hard claim | `hard_m` = hard_support on reject (mild) |
| Dual | mild→hard_support; intensity>0.25→soft CV |
| RFPerm-era extras | blend mix≤0.5; reject dual_r2 / n_recent=2 |

## Evidence

- v8.1: Jaccard(CFPerm,RFPerm)≈0.00–0.22; sparse stocks hard_m pack −3%…−4%
- v8.2 scan: beijing_gate≤0.25 → bj@reject≈1.0, dual pack rel ≈ **+0.6%**; gate=0.45 → bj@reject≈0.2, dual pack rel ≈ **+9%**

## Open

- Full 40×256 multi-seed under beijing_gate=0.25
- dual+blend_50 hybrid only when BJ fires
- Optional: post-reject FSDS / VIMP feature subset
