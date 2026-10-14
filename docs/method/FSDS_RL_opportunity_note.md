# FSDS → RL closed-loop: dataset opportunities

Locked TSS signs stay fixed: **covariate ↑ → shrink η / IW / memory**;
**concept ↑ → raise η / relearn / refresh reward heads**.

## Opportunity matrix

| dataset | FSDS channel | best-fit module | why it is a closed-loop opportunity |
| --- | --- | --- | --- |
| Meta-World ML1/10/45 | task-param cov + sparse reward concept | meta_adapt + region_lr | task descriptor = top-K FSDS features → few-shot inner loop |
| Procgen (train→test levels) | layout/style cov; gem/progress concept | reward_refresh + region_lr | OOD levels are cov; sparse-reward games flip concept mid-curriculum |
| D4RL (umaze→medium/large) | state-support cov | region_lr (shrink) + IW | offline support mismatch; naive η↑ hurts |
| CausalWorld | mass/friction cov; objective concept | all three | causal params → descriptor; safety objective flip → reward mix |
| Meta-Dataset | domain/style cov; class map concept | meta_adapt | cross-domain few-shot; descriptor conditions init |
| MEWA | collab-mode concept | reward_refresh | partner/protocol change = reward-head staleness |
| *existing* MSR-VTT modality π | modality cov/concept | region_lr (already) | reuse as RL-adjacent control prior |
| *existing* Amazon hop IW | category cov / rating concept | reward_refresh analogy | mix bank vs re-fit mirrors reward-head refresh |

## Proxy boards in this PR

Full sims are not installed in the cloud checkout. Proxies keep the **same
typed taxonomy** so the three modules are executable today:

| proxy | stands in for |
| --- | --- |
| `bandit_{cov,concept,mixed}` | Meta-World task / reward drift |
| `maze_proxy` | D4RL support mismatch |
| `causal_proxy` | CausalWorld mass/friction + objective flip |
| `procgen_proxy` | Procgen layout OOD + gem reward flip |
| `meta_vision` | Meta-Dataset domain hop + class map rotation |

## Extra opportunities (beyond the six listed)

1. **Offline→online handoff** (D4RL medium-replay → online fine-tune): FSDS on
   the first online hop decides IW vs concept re-fit before any env steps.
2. **Multi-objective safety RL**: treat safety / progress / energy as reward
   components; concept drift on safety head triggers weight refresh (CausalWorld /
   Safety-Gym pattern).
3. **Curriculum generators** (Procgen / MiniGrid): FSDS on level features gates
   when to raise vs shrink η — avoids naive “harder level ⇒ always boost LR”.
4. **Reuse MSR-VTT / Amazon boards as non-RL priors**: modality π and hop IW
   already implement module-3 / module-2 polarity; RL boards should not invent a
   new loss, only new actuators.

## What to measure next on real sims

- Adapt steps to target return (meta).
- Reward–true-objective correlation after refresh.
- OOD return gap vs `naive_boost` (must stay positive under covariate boards).
- Detection latency of mid-stream concept flips (CausalWorld / Procgen gem).

## Proxy-suite findings (this checkout)

Measured offline IW value + reward-mix MSE on 7 boards × 4 seeds
(`scripts/run_fsds_rl_closed_loop.py`):

| signal | finding |
| --- | --- |
| vs `static` | every FSDS module lifts policy value; largest lifts on `bandit_concept` / `bandit_mixed` |
| reward refresh | cuts `mix_mse` by 10–100× on bandit / maze / meta_vision when holdout accepts the empirical mix |
| `fsds_full` | strongest on `causal_proxy` (objective flip); else `meta_only` ≈ best for IW value |
| `naive_boost` | still competitive on *policy* IW — TSS polarity shows up more cleanly on **mix_mse / w_err** than on policy value in these linear proxies |
| holdout gate | required: unconstrained Ridge refresh overfit on `causal_proxy` and raised mix_mse |

**Practical closed-loop rule:** use FSDS `(c,δ)` to *choose the actuator*
(meta inner-loop vs reward-mix re-fit vs region η), not to invent a new loss.
Policy-η boost alone is an incomplete story under pure covariate boards — pair
it with IW / memory when swapping in real D4RL / Procgen loaders.
