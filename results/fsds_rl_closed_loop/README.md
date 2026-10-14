# FSDS-empowered RL closed loop

Three modules driven by FSDS typed drift `(c, δ)` — **locked signs**:
covariate large → shrink η; concept large → raise η / refresh reward heads.

## Modules

| id | role |
| --- | --- |
| meta_adapt | few-shot policy inner-loop conditioned on drift descriptor |
| reward_refresh | mix weights over reward components from typed drift |
| region_lr | TSS η per feature group |

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

| method | value ↑ | w_err ↓ | mix_mse ↓ |
| --- | --- | --- | --- |
| static | 2.1768 (0.2575) | 0.2160 | 0.3869 |
| lr_boost_only | 2.4092 (0.2148) | 0.2160 | 0.3869 |
| reward_only | 2.3574 (0.2256) | 0.0166 | 0.0088 |
| meta_only | 2.4421 (0.2385) | 0.2160 | 0.3869 |
| fsds_full | 2.4298 (0.2223) | 0.0166 | 0.0088 |
| naive_boost | 2.4535 (0.2281) | 0.2160 | 0.3869 |

### `bandit_concept`

| method | value ↑ | w_err ↓ | mix_mse ↓ |
| --- | --- | --- | --- |
| static | 2.6171 (0.1627) | 0.5753 | 2.9835 |
| lr_boost_only | 2.9770 (0.1761) | 0.5753 | 2.9835 |
| reward_only | 2.8521 (0.1572) | 0.0408 | 0.0276 |
| meta_only | 2.9961 (0.1789) | 0.5753 | 2.9835 |
| fsds_full | 2.9832 (0.1777) | 0.0408 | 0.0276 |
| naive_boost | 2.9984 (0.1817) | 0.5753 | 2.9835 |

### `bandit_mixed`

| method | value ↑ | w_err ↓ | mix_mse ↓ |
| --- | --- | --- | --- |
| static | 2.7624 (0.3620) | 0.5753 | 3.2165 |
| lr_boost_only | 3.1019 (0.3370) | 0.5753 | 3.2165 |
| reward_only | 3.0098 (0.3461) | 0.0420 | 0.0278 |
| meta_only | 3.1363 (0.3469) | 0.5753 | 3.2165 |
| fsds_full | 3.1196 (0.3431) | 0.0420 | 0.0278 |
| naive_boost | 3.1447 (0.3545) | 0.5753 | 3.2165 |

### `maze_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ |
| --- | --- | --- | --- |
| static | 0.7813 (0.0431) | 0.4546 | 0.3003 |
| lr_boost_only | 0.9272 (0.0379) | 0.4546 | 0.3003 |
| reward_only | 0.9035 (0.0329) | 0.0931 | 0.0051 |
| meta_only | 0.9494 (0.0454) | 0.4546 | 0.3003 |
| fsds_full | 0.9460 (0.0474) | 0.0931 | 0.0051 |
| naive_boost | 0.9591 (0.0461) | 0.4546 | 0.3003 |

### `causal_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ |
| --- | --- | --- | --- |
| static | -0.4389 (0.0222) | 0.3342 | 0.1217 |
| lr_boost_only | -0.3737 (0.0231) | 0.3342 | 0.1217 |
| reward_only | -0.3888 (0.0240) | 0.3527 | 0.2928 |
| meta_only | -0.3675 (0.0203) | 0.3342 | 0.1217 |
| fsds_full | -0.3626 (0.0219) | 0.3527 | 0.2928 |
| naive_boost | -0.3616 (0.0210) | 0.3342 | 0.1217 |

### `procgen_proxy`

| method | value ↑ | w_err ↓ | mix_mse ↓ |
| --- | --- | --- | --- |
| static | -0.0206 (0.0123) | 0.3894 | 0.0320 |
| lr_boost_only | 0.0412 (0.0122) | 0.3894 | 0.0320 |
| reward_only | 0.0150 (0.0162) | 0.1699 | 0.0136 |
| meta_only | 0.0530 (0.0157) | 0.3894 | 0.0320 |
| fsds_full | 0.0424 (0.0197) | 0.1699 | 0.0136 |
| naive_boost | 0.0514 (0.0167) | 0.3894 | 0.0320 |

### `meta_vision`

| method | value ↑ | w_err ↓ | mix_mse ↓ |
| --- | --- | --- | --- |
| static | 3.4954 (1.0497) | 0.5435 | 7.6288 |
| lr_boost_only | 3.6564 (1.1459) | 0.5435 | 7.6288 |
| reward_only | 3.5814 (1.1369) | 0.0616 | 0.1113 |
| meta_only | 3.6648 (1.1394) | 0.5435 | 7.6288 |
| fsds_full | 3.6339 (1.1438) | 0.0616 | 0.1113 |
| naive_boost | 3.6773 (1.1213) | 0.5435 | 7.6288 |

## Run

```bash
python3 scripts/run_fsds_rl_closed_loop.py
python3 scripts/run_fsds_rl_closed_loop.py --quick
```

Code: `Python/src/fsds_rl_closed_loop.py`

## Opportunity notes

- Install `gymnasium` + Meta-World / Procgen / D4RL / CausalWorld to replace proxies.
- Existing FSDS boards (MSR-VTT modality π, Amazon hop IW) already implement
  module-3 style control; this package adds reward refresh + meta inner-loop.
- Do **not** use `naive_boost` in production — it raises η on covariate shift.
