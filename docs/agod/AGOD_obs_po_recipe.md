# AGOD obs-PO recipe (CFPerm gate, v8.4)

## Thesis

认 **hard**; packMSE only under **beijing**; **L0 = CFPerm DRPerm**.

## Locked

| Item | Choice |
|---|---|
| L0 | CFPerm DRPerm, `e_mode=known`, `n_perm≥39`, α=0.05 |
| L1 | intensity → λ / beijing; **`beijing_gate=0.25`** |
| temper | gate=0.20, lam_max=0.75 |
| Default policy | **`dual`** (mild→hard_support; intensity>0.25→soft CV) — pending dual_b50 ablation |
| Hard claim path | `hard_m` reported separately |
| Blend | mix≤0.5; optional on sparse packs |

## Evidence (v8.3: 40×256 × 6 packs × 3 seeds)

- Synth: size=0, power=0.8
- bj@reject = **1.0** on all packs with rejects
- Pack vs uni: dual **+1.8% / +0.4%** on metro/beijing vs hard_m **+5.2% / +8.0%**
- Hard-subset: dual wins 3/6 packs

## v8.4 in progress

- **dual_b50**: BJ rejects → blend_50 instead of soft CV
- Decide: keep dual if dual_b50 loses metro/beijing packMSE

## Open (after v8.4)

- Post-reject FSDS / VIMP feature subset
- Optional RRPerm L0 ablation
- temper/lam schedule polish
