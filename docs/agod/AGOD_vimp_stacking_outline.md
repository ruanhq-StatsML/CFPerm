# VIMP / FSDS → α → stack\_w (portable online stacking)

## Statistical intuition

MoE aux pushes gates toward **internal balance**.
Here the prior is an **external attributor**:

- PO / residual concept (pseudo-outcome)
- MMD / distance (covariate)
- VIMP / Fisher (importance gate)

So `α` is not “use experts evenly”; it is “who carries concept drift mass”.
Online stacking then spends fusion mass `w` on that attributor.

## Algorithm outline

**Input:** modality features `{X_m}`, labels `y`, stream windows `t=1..T`.

**Sensors (per window):**
1. Estimate per-modality PO, MMD, VIMP (FSDS / MSG).
2. Score `g_m = f(PO_m, VIMP_m) − λ·MMD_m` (concept↑, covariate↓).
3. `α_raw = Softmax(g/τ)`, `α ← EMA(α)`.

**Online stacking:**
4. Towers → logits `z_m`; `w = Softmax(ψ)`.
5. `ẑ = Σ_m w_m z_m`.
6. Loss:
   - `mean_ce`: CE on mean-pool (no `w`)
   - `stack_ce`: CE(`ẑ`,`y`)
   - `stack_alpha`: CE(`ẑ`,`y`) + `λ_kl·KL(w ‖ α)`   ← weight socket

**Output:** holdout Brier/MSE drop, Acc lift, `‖w−α‖`.

## Portable / naive use-case

Any multimodal online learner with (i) an attribution score and (ii) a
simplex fusion weight can reuse the same socket: **attributor → α → KL to `w`**.
The actuator is ~20 LOC on top of CE stacking — no new backbone.

Paper-ready note (tables + code listing):
`docs/agod/AGOD_naive_stack_socket.tex`

Code: `agod/online_stack.py` (`alpha_stack_kl`, `WEIGHT_MODES`, `stack_weight_aux`),
`scripts/run_agod_online_stack_compare.py`.
