# FSDS Reasoning Empowerment — Adjust Next Step

> Skill: `fsds_reasoning_empowerment` v1 · 用 FSDS 归因节点特征，Guide **下一步**（剪枝/扩展/重排序/回溯），并给出可审计的经济 ledger。

## 1. 问题

推理搜索（ToT / RAP / MCTS / Beam…）每一步都要回答：

- 哪些节点扩展？哪些剪枝？何时回溯？如何重排？

本仓答案：**FSDS 用节点特征 \(X\) 预测 \(Y=1\{\text{node on optimal path}\}\)** → `importance(i)∝P(Y_i=1)` → Guide 决策。

## 2. 核心链路

```
轨迹节点 X (struct∥semantic∥score∥search∥cost∥time∥env)
  → Overlap router (§5.5)
      Overlap<0.3 → MMD-LOCO + RF-Binary-VIMP
      Overlap≥0.3 → PO-risk-LOCO + PermuCATE (+ RF-VIMP fuse)
  → fused feature importance + RF P(Y=1|X)
  → Guide: prune / expand / reorder / backtrack
  → 对比 score-greedy baseline
  → NetValue / ROI ledger
```

## 3. 实现

| 件 | 路径 |
|---|---|
| 核心 | `agod/fsds_reason_guide.py` |
| 记分卡 | `scripts/run_fsds_reason_guide.py` |
| 测试 | `tests/test_fsds_reason_guide.py` |
| 产物 | `results/fsds_reason_guide/` |
| **LaTeX prototype (RAP)** | `docs/method/FSDS_Reasoning_Empowerment_RAP.tex` |

与 `agod/tot_eff.py`（策略 ToT：adapt×α×freeze）正交：本模块是 **推理节点级** 的 FSDS→下一步；`tot_eff` 是 **效率策略** ToT。

## 4. 经济公式（skill §7）

```
Revenue     = 1{success} · V_task
Cost        = tokens·P_token + calls·P_call + ½·waste_off_opt
NetValue    = Revenue − Cost
Incremental ≈ ΔRevenue + ΔCost_savings + ΔAUC·V_auc − Cost_FSDS
ROI         = Incremental / invest
```

**Y 是经济价值的直接驱动变量**：正确识别有效节点 → 成功↑收入↑；正确剪枝无效节点 → 浪费↓成本↓。

## 5. Justification（为何这样接）

1. **为何用 FSDS 而不是纯 score**：score 有 decoy；FSDS 融合结构/语义/搜索行为，在分布偏移下仍可学 \(X→Y\)。
2. **为何 Y=最优路径指示**：可标注、细粒度、可干预（剪/扩改变路径）、与任务成功对齐。
3. **为何 Guide 四动作**：与 skill §6 一致；优先级 prune > backtrack > expand > reorder，对应成本/失败/收入/效率。
4. **ROI 边界**：当前数字来自 **合成 RAP/ToT 记分卡**（已知 oracle 路径）。证明的是公式可算、Δ 可测、符号与机制一致；**不是**生产美元 ROI 的因果识别。

## 6. 跑法

```bash
PYTHONPATH=. python3 scripts/run_fsds_reason_guide.py --n-train 40 --n-test 30
PYTHONPATH=. python3 -m pytest tests/test_fsds_reason_guide.py -q
```

## 7. 与仓库其它件的关系

- Tabular FSDS / MMD-LOCO / PO：`scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py`
- 策略 ToT：`agod/tot_eff.py` · `docs/summaries/Stats_Empower_Tree_of_Thoughts.md`
- 本 skill：**把同一套 FSDS 归因接到「下一步」**，并显式算增量净价值。
