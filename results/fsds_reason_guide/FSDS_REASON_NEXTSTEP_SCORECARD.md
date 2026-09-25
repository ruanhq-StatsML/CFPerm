# FSDS Reasoning Empowerment — Adjust Next Step

> Skill: `fsds_reasoning_empowerment` · Y=1{node on oracle path} · Guide={prune,expand,reorder,backtrack}

## Verdict (simulation ledger)

- Overlap=0.045 → regime=covariate_shift (skill §5.5 router).
- FSDS P(Y=1|X) AUC=0.995 vs score-only AUC=0.822 on train nodes.
- Test success: greedy 13.3% → guided 93.3% (Δ=+80.0%).
- Tokens/task: 824 → 819 (Δ saved=+5).
- Net/task ledger: -1.434 → 2.823 (Δ=+4.257).
- Incremental NetValue≈4.453/task; ROI≈29.69 (sim ledger, not causal $).
- Y=1{node on oracle path} is the direct driver of task success → revenue; Guide maps importance→prune/expand.

## Regime & attribution

- **overlap** = 0.0453
- **regime** = `covariate_shift`
- **components** = ['MMD-LOCO', 'RF-Binary-VIMP']
- **AUC P(Y=1|X)** = 0.9945

### Top fused features

| feature | importance |
|---|---:|
| `state_entropy` | 0.0528 |
| `state_length` | 0.0481 |
| `tokens` | 0.0345 |
| `is_root` | 0.0329 |
| `score` | 0.0322 |
| `depth` | 0.0319 |
| `path_length` | 0.0318 |
| `has_keywords` | 0.0317 |

## Search effect (held-out tasks)

n_train=40 · n_test=30

| policy | success | tokens/task | cost/task | net/task | waste_off_opt |
|---|---:|---:|---:|---:|---:|
| score-greedy | 13.3% | 824.3 | 2.101 | -1.434 | 0.565 |
| FSDS-guided | 93.3% | 819.1 | 1.844 | 2.823 | 0.071 |
| **Δ** | +80.0% | +5.2 | +0.257 | +4.257 | +0.494 |

### Guide action counts (test)

```
{
  "prune": 12,
  "expand": 170,
  "reorder": 448,
  "backtrack": 0
}
```

## Economic ledger (skill §7)

```
NetValue = Σ Y_i · V_task · N_tasks − Σ (1−Y_i) · C_node
Incremental ≈ ΔRevenue + ΔCost_savings + ΔAUC·V_auc − Cost_FSDS
```

| term | value / task |
|---|---:|
| ΔRevenue (via success) | +4.0000 |
| ΔCost savings | +0.2572 |
| ΔAUC · V_auc | +0.3458 |
| Cost_FSDS | −0.1500 |
| **Net incremental** | **+4.4530** |
| **ROI** (Δ / invest) | **29.69** |

Params: `{'v_task': 5.0, 'c_token': 0.002, 'c_call': 0.01, 'v_auc': 2.0, 'cost_fsds': 0.15}`

## Justification

1. **Why adjust next step with FSDS**: search nodes carry X (struct/semantic/score/…) and a labelable Y (on-oracle-path). FSDS predicts Y under shift; importance∝P(Y=1) is exactly the signal Guide needs for expand vs prune.
2. **Why ROI is computable here**: success→V_task, off-opt nodes→C_node, tokens/calls metered; Δ vs score-greedy is an A/B-style ledger on the same trees.
3. **What this does *not* claim**: dollar ROI in production RAP/LATS, or that fused importance is a causal CATE without further identification. Numbers are **simulation scorecard** under the skill's NetValue formula.
4. **Incremental value**: positive when guided raises success and/or cuts wasted tokens enough to beat Cost_FSDS; see table above.

