# 古典预测 × Agent 飞轮：机会在哪

> 接 sandbox 预测 bakeoff；**不**接 AGOD excess/PO。  
> 问：古典预测的增益，怎么变成 **agent 飞轮** 的可迭代机会？

代码：`sandbox/forecast_flywheel.py` · bakeoff：`sandbox/forecast_bakeoff.py`

---

## 1. 古典预测告诉我们的事实（已跑）

| pack | 古典赢家 | 含义 |
|---|---|---|
| metro | **HGB**（lift≈+43%） | 有可学结构 → **值得养一个预测 agent** |
| waymo | **Ridge** | 线性就够 → **便宜 agent**，别上重模型 |
| pm25 | **naive**（HGB 负 lift） | 近随机游走 → **不要浪费 HGB agent** |

机会 **不是**「再堆一个更复杂的模型」，而是：

> **按 pack 路由模型身份** —— 飞轮的第一刀是 *谁上场*，不是 *谁更深*。

---

## 2. Agent 飞轮长什么样（古典版）

```text
        ┌── observe window ──┐
        ▼                    │
   classical forecast        │
        │                    │
   residual / surprise       │
        │                    │
   decide: idle | log | retrain | fallback_naive
        │                    │
        └────── log ─────────┘
```

| 环节 | 古典预测贡献 | Agent 机会 |
|---|---|---|
| Observe | 时间序窗 | 每个 pack 一个 watcher agent |
| Forecast | naive/Ridge/HGB | **专科 agent**（按 bakeoff 选） |
| Surprise | \|resid\|/σ | 触发人审 / 换特征 / 记日志 |
| Retrain | 周期性 refit | 飞轮「自更新」——仍是古典回归 |
| Fallback | rolling lift\<0 | 自动降级 naive，防瞎烧算力 |

这和「180 人商业车道」无关：这里飞轮 = **预测误差驱动的自更新环**。

---

## 3. 机会清单（可派）

| ID | 机会 | 为什么现在能说 | 飞轮动作 |
|---|---|---|---|
| FW-metro-hgb | metro 专职 HGB agent | bakeoff lift 大 | surprise→retrain |
| FW-waymo-ridge | waymo 线性 agent | ridge 赢 HGB | 在线/定期 ridge |
| FW-pm25-naive | pm25 只跑 naive | HGB 无增益 | 只监 surprise（regime） |
| FW-router | pack→model 路由 agent | 三包三赢家 | 读 bakeoff 卡派工 |
| FW-theme | DiffDB theme 上下文 | k-means 词簇 | 给生成/审出 agent 贴 cluster id（旁路） |

**P0：** router + metro HGB + waymo ridge。  
**P1：** pm25 naive 看守（防 regime 变了还当 RW）。

---

## 4. 效果怎么谈（诚实）

飞轮跑出来的数是：

- `surprise_rate`：多大比例步触发异常  
- `n_retrain`：自更新次数  
- `fell_back_to_naive`：HGB 是否被滚动 lift 拉下马  

**能说：** 古典预测增益 → 专科 agent 分工与降级策略。  
**不能说：** 因果、PO 效率、商业 ROI（那是另一条线）。

---

## 6. DPO + 标签污染 + decision-path 回放

喂飞轮的一条具体方法（已实现 `sandbox/flywheel_dpo_replay.py`）：

1. **每步 JSON**：`{step,t,y_true,y_hat,residual,surprise,action,model,decision_path}`  
2. **区域污染**：随机抽若干 period → `indices_subset` →  
   `y[idx] = y[np.random.permutation(idx)]`  
3. **回放**：对每步 `decision_path` 中段 shuffle（保留 observe/log 端点）  
4. **DPO**：chosen=干净路径奖励，rejected=污染路径 / 打乱路径；  
   \(L=-\log\sigma(\beta(r_c-r_r))\)

跑：`PYTHONPATH=. python3 scripts/run_flywheel_dpo_replay.py`

---

## 7. 10 分钟迭代怎么用

每 tick：抽一个 FW-* 机会 → 改 flywheel 决策阈值或加一个 pack → 看 surprise/fallback 变没变 → 记一行。  
约束：**不许**把 AGOD null/PO 塞进这个飞轮。
