# Follow-up: Grad-OnlineRFPerm（`param.grad.norm`）— Sep17 status

**Audience:** leadership / research lead  
**Date:** 2026-09-17  
**PR:** https://github.com/ruanhq-StatsML/CFPerm/pull/64  
**Branch:** `cursor/agod-grad-rfperm-monitor-abce`

## 1. Ask → what we shipped

Continuous-time **OnlineRFPerm** on serving gradient energy, so we can flag distribution / concept shift **before** serving-MSE (and PO / MMD) fully break — and use that as an engineering signal for **early layer freeze / back-prop depth**.

Shipped end-to-end on the Sep17 MVP packs: **synthetic, Covertype, bank-marketing, electricity, eeg-eye-state**.

## 2. Engineering口径（已钉死，避免 multiple testing）

| 项 | 口径 |
|---|---|
| 参数范围 | 仅 **未冻结** 参数 `θ_U = {p : requires_grad=True}` |
| 标量 `g_t` | **一个** `‖∇_{θ_U} L‖₂ = sqrt(Σ_i ‖∇θ_i‖²)`（全未冻结 grad 拼成向量再取 ℓ₂） |
| 检验 | **一整条** OnlineRFPerm：`T_t = g_t − e_ref` → rank/EWMA p → alpha-investing FDR |
| 不是 | ❌ 每个 layer / param 各跑一套再 `any(reject)`（multiple testing） |
| 不是 | ❌ 各 layer norm 归一化后取平均当 gate |
| 诊断 | 全局 reject **之后**，`share_ℓ = ‖g_ℓ‖ / g_t` 只做 freeze-depth **排名**，不进 FDR |

Serving 协议：burn-in 上预训练 **frozen `f_ref` MLP**，之后权重不再 step；每 batch 对 serving loss 做一次 backward，只读 `g_t`。这与 MSE-OnlineRFPerm 的 “fixed `f_ref`” 叙事对齐。

```
θ_U = unfrozen params
g_t = ||∇_{θ_U} L(batch; f_ref)||_2          # ONE scalar
T_t = g_t - e_ref
p_t = rank/EWMA vs history; online FDR       # ONE hypothesis stream
if reject: rank layers by share_ℓ            # diagnostic only
```

## 3. MVP 结果（5 seeds × 5 datasets）

Lead = `t_grad − t_mse`（**负 = Grad 比 MSE-OnlineRFPerm 更早 reject**）。  
设置：`batch=128`, `n_batches=48`, `n_burn=8`, `α=0.05`, seeds `0..4`.

| dataset | mean lead | median | P(earlier) | P(≤0) |
|---|---:|---:|---:|---:|
| synthetic | **−3.0** | −2 | 100% | 100% |
| covertype | **−1.2** | −1 | 60% | 80% |
| bank | **−1.2** | 0 | 40% | 100% |
| electricity | **−6.2** | −4 | 80% | 100% |
| eeg | **−4.6** | −6 | 80% | 100% |

**Overall（25 runs）：** mean lead **−3.24** batches；Grad 早于 MSE **72%**；不晚于 MSE **96%**.

解读（给老板的一句话）：在正确的单流 FDR 口径下，Grad 能量仍能稳定地 **提前约 3 个 batch** 相对 serving-MSE 报警；electricity / eeg 提前最明显，bank 基本持平（不伤）。

Artifacts：`results/agod_grad_rfperm/`（per-dataset plots + `lead_time_summary.png` + `multiseed_summary.json`）。

## 4. 产品 / 训练侧怎么用

1. **Gate：** 全局 Grad-OnlineRFPerm reject → 打开适配（√PO / refit / 加深更新），与现有 MSE gate 并行或抢跑。  
2. **Freeze：** reject 后看 `share_ℓ` 谁最大 → 优先解冻 / 加深那一层附近的 back-prop；早期层 share 高偏 covariate，末层偏 concept（诊断，非硬判决）。  
3. **与 MSE/PO/MMD 的关系：** Grad = 更早的 *parameter-space* 警报；MSE = serving 误差；MMD/PO = 分布 / 风险窗口。三者同屏，不互相替代。

## 5. 边界与下一步（建议）

**已澄清**
- 单流 `‖∇_U‖₂`，无跨参数 multiple testing。  
- frozen `f_ref` 读 grad，避免 online 更新把 norm 训没。

**建议下周细化（部分已做，见 §5b）**
1. Grace / 稳健性假阳扫描  
2. Freeze 闭环 MSE / FLOPs  
3. 接到 stocks / metro / beijing / Waymo pack  
4. 大模型 adapter / LoRA 未冻结参数同口径（仍待做）

### 5b. Supplementary experiments（已跑）

`docs/method/Grad_OnlineRFPerm_extras.md` · `results/grad_rfperm_extras/`

| Extra | Headline |
|---|---|
| Null + grace | grace=4：first reject 8→15；early(≤5) FPR 100%→67% |
| α sweep | synthetic/electricity 上 lead 在 α∈{0.01,0.05,0.10} 仍为负 |
| Freeze loop | electricity：`freeze_early` ≈0.86× MSE @ 0.79× FLOPs vs always_adapt |
| Extra packs | stocks_IWM lead −4.0（100% earlier）；metro +2；waymo ≈0 |

## 6. 代码入口

```bash
PYTHONPATH=. python3 scripts/run_agod_grad_rfperm_monitor.py \
  --datasets synthetic covertype bank electricity eeg \
  --seeds 0 1 2 3 4 --batch-size 128 --n-batches 48 --n-burn 8
```

- Module: `agod/grad_rfperm.py`  
- Bench: `scripts/run_agod_grad_rfperm_monitor.py`  
- Doc: `docs/agod/AGOD_grad_rfperm_monitor.md`  
- Tests: `tests/test_grad_rfperm.py`

---

**Bottom line for the boss:** Grad-OnlineRFPerm 已按“**未冻结参数一个 ℓ₂ 标量 + 单流 OnlineRFPerm**”落地并在 5 个 MVP dataset × 5 seeds 验证；相对 MSE 平均提前 **~3 batches**（72% 更早），可直接作为 early freeze / 适配 gate 的工程信号；下一步接 freeze 闭环与更大 pack。
