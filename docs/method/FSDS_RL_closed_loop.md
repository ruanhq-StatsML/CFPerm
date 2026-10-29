# FSDS-empowered RL closed loop

Three modules driven by FSDS typed drift `(c, δ)` — **locked signs**:
covariate large → shrink η; concept large → raise η / refresh reward heads.

## Modules

| id | role |
| --- | --- |
| meta_adapt | few-shot policy inner-loop conditioned on drift descriptor |
| reward_refresh | mix weights over reward components from typed drift |
| region_lr | TSS η per feature group |
| meta_iv | PO-risk uplift monitor + Important-Variable selection gates δ / η |

## Boards (proxies — full sims not installed here)

| board | stands in for | drift |
| --- | --- | --- |
| bandit_cov | Meta-World task params | covariate |
| bandit_concept | reward / map change | concept |
| bandit_mixed | both | both |
| maze_proxy | D4RL support mismatch | covariate |
| causal_proxy | CausalWorld mass/friction + objective flip | both |
| procgen_proxy | Procgen layout/style OOD + gem reward flip | both |
| meta_vision | Meta-Dataset domain hop + class map rotation | both |

## Results (offline IW value ↑ better)

### `bandit_cov`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO |
| --- | --- | --- | --- | --- | --- |
| static | 2.1768 (0.2575) | 0.2160 | 0.3869 | nan | nan |
| fsds_full | 2.4298 (0.2223) | 0.0166 | 0.0088 | nan | nan |
| fsds_meta_iv | 2.4147 (0.2130) | 0.0166 | 0.0088 | 7.0 | 0.0244 |
| naive_boost | 2.4535 (0.2281) | 0.2160 | 0.3869 | nan | nan |

### `bandit_concept`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO |
| --- | --- | --- | --- | --- | --- |
| static | 2.6171 (0.1627) | 0.5753 | 2.9835 | nan | nan |
| fsds_full | 2.9832 (0.1777) | 0.0408 | 0.0276 | nan | nan |
| fsds_meta_iv | 2.9610 (0.1773) | 0.0408 | 0.0276 | 7.0 | 0.0394 |
| naive_boost | 2.9984 (0.1817) | 0.5753 | 2.9835 | nan | nan |

### `bandit_mixed`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO |
| --- | --- | --- | --- | --- | --- |
| static | 2.7624 (0.3620) | 0.5753 | 3.2165 | nan | nan |
| fsds_full | 3.1196 (0.3431) | 0.0420 | 0.0278 | nan | nan |
| fsds_meta_iv | 3.0828 (0.3383) | 0.0420 | 0.0278 | 7.0 | 0.0394 |
| naive_boost | 3.1447 (0.3545) | 0.5753 | 3.2165 | nan | nan |

### `maze_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO |
| --- | --- | --- | --- | --- | --- |
| static | 0.7813 (0.0431) | 0.4546 | 0.3003 | nan | nan |
| fsds_full | 0.9460 (0.0474) | 0.0931 | 0.0051 | nan | nan |
| fsds_meta_iv | 0.9258 (0.0358) | 0.0931 | 0.0051 | 6.0 | 0.0046 |
| naive_boost | 0.9591 (0.0461) | 0.4546 | 0.3003 | nan | nan |

### `causal_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO |
| --- | --- | --- | --- | --- | --- |
| static | -0.4389 (0.0222) | 0.3342 | 0.1217 | nan | nan |
| fsds_full | -0.3626 (0.0219) | 0.3527 | 0.2928 | nan | nan |
| fsds_meta_iv | -0.3438 (0.0228) | 0.3527 | 0.2928 | 6.0 | 0.0038 |
| naive_boost | -0.3616 (0.0210) | 0.3342 | 0.1217 | nan | nan |

### `procgen_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO |
| --- | --- | --- | --- | --- | --- |
| static | -0.0206 (0.0123) | 0.3894 | 0.0320 | nan | nan |
| fsds_full | 0.0424 (0.0197) | 0.1699 | 0.0136 | nan | nan |
| fsds_meta_iv | 0.0440 (0.0197) | 0.1699 | 0.0136 | 7.0 | 0.0011 |
| naive_boost | 0.0514 (0.0167) | 0.3894 | 0.0320 | nan | nan |

### `meta_vision`

| method | value ↑ | w_err ↓ | mix_mse ↓ | n_iv | R_PO |
| --- | --- | --- | --- | --- | --- |
| static | 3.4954 (1.0497) | 0.5435 | 7.6288 | nan | nan |
| fsds_full | 3.6339 (1.1438) | 0.0616 | 0.1113 | nan | nan |
| fsds_meta_iv | 3.6358 (1.1578) | 0.0616 | 0.1113 | 7.0 | 0.0431 |
| naive_boost | 3.6773 (1.1213) | 0.5435 | 7.6288 | nan | nan |

## Run

```bash
python3 scripts/run_fsds_rl_closed_loop.py
python3 scripts/run_fsds_rl_closed_loop.py --quick
python3 scripts/run_fsds_rl_closed_loop.py --compare-meta-iv
```

Code: `Python/src/fsds_rl_closed_loop.py`

## Opportunity notes

- Install `gymnasium` + Meta-World / Procgen / D4RL / CausalWorld to replace proxies.
- Existing FSDS boards (MSR-VTT modality π, Amazon hop IW) already implement
  module-3 style control; this package adds reward refresh + meta inner-loop.
- `fsds_meta_iv` reuses the repo meta-learner (PO-risk uplift) + IV selection
  to focus TSS / descriptor on shift-driving features.
- Do **not** use `naive_boost` in production — it raises η on covariate shift.
