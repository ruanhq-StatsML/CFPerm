# IPTW power softness: √ vs ∛ (same FLOPs)

gated∛≤gated√ on 67% of packs; median Δmse_eff(∛−√)=+0.000227 (mean=-7.59e+03, metro-skewed); best_sig counts {'gated_sqrt': 1, 'uniform': 5}. Same FLOPs — softness is the only knob.

sqrt vs cbrt IPTW share adaptation FLOPs; compare rel MSE / soft wins, not compute. Prefer ∛ when PO ranks well but √ IPTW hurts. Δmse_eff(∛−√)>0 ⇒ ∛ loses less (or gains more) per FLOP. Prefer **median** gap over mean — one pack (metro) can dominate.

- batches=40 · batch_size=100 · datasets=6 · ∛≤√ rate=0.6666666666666666

| dataset | duty | unif | gated√ rel | gated∛ rel | gated√ mse_eff | gated∛ mse_eff | ∛≤√? | best | burn | reading |
|---|---:|---:|---:|---:|---:|---:|:---:|---|---|---|
| `metro_interstate` | 0.1795 | 8.602e+05 | 0.9274 | 0.9340 | 4.835e+05 | 4.394e+05 | N | `gated_sqrt` | `BURN_SQRT` | duty=0.18: gated√ best sig MSE; √ beats ∛ here; Δmse_eff(∛−√)=-4.41e+04 |
| `beijing_pm25` | 0.1282 | 2364 | 1.0118 | 1.0676 | -301.1272 | -1731 | N | `uniform` | `SOFTEN_ONLY` | duty=0.13: keep uniform — gated IPTW does not buy sig MSE; √ beats ∛ here; Δmse_eff(∛−√)=-1.43e+03 |
| `stocks_AAPL` | 0.1795 | 0.0007784 | 1.1450 | 1.1028 | -0.0008735 | -0.0006193 | Y | `uniform` | `SOFTEN_ONLY` | duty=0.18: keep uniform — gated IPTW does not buy sig MSE; ∛≤√ on this pack; Δmse_eff(∛−√)=+0.000254 |
| `waymo_proxy` | 0.5641 | 0.009992 | 1.0832 | 1.0496 | -0.002048 | -0.00122 | Y | `uniform` | `SOFTEN_ONLY` | duty=0.56: keep uniform — gated IPTW does not buy sig MSE; ∛≤√ on this pack; Δmse_eff(∛−√)=+0.000828 |
| `stocks_MSFT` | 0.2308 | 0.0006164 | 1.1154 | 1.0257 | -0.000428 | -9.522e-05 | Y | `uniform` | `SOFTEN_ONLY` | duty=0.23: keep uniform — gated IPTW does not buy sig MSE; ∛≤√ on this pack; Δmse_eff(∛−√)=+0.000333 |
| `stocks_IWM` | 0.1795 | 0.0005025 | 1.0904 | 1.0388 | -0.0003513 | -0.0001507 | Y | `uniform` | `SOFTEN_ONLY` | duty=0.18: keep uniform — gated IPTW does not buy sig MSE; ∛≤√ on this pack; Δmse_eff(∛−√)=+0.000201 |

## Burn policy (α ≠ FLOPs)

- See [`Soft_Weight_Burn_Logic.md`](../../docs/summaries/Soft_Weight_Burn_Logic.md).
- Default: do not burn; `SOFTEN_ONLY` → always ∛ if gate is mandatory.
- Burn only when gated_α beats uniform on sig-MSE.

## Tomorrow-demo takeaway

1. Softness ≠ free lunch: uniform still wins most packs.
2. When you *do* gate IPTW, ∛ ≤ √ on a majority of packs here.
3. FLOPs identical → choose α by MSE risk, not by compute.

