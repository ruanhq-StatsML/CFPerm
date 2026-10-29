# FSDS → RL closed-loop: dataset opportunities

Locked TSS signs stay fixed: **covariate ↑ → shrink η / IW / memory**;
**concept ↑ → raise η / relearn / refresh reward heads**.

## Product case (intuitive)

FSDS does not invent a new RL loss. It reads a typed hop `(c, δ)` and
picks **which actuator** to move: meta inner-loop, reward-mix weights, or
region learning rate. The product story is three linked claims.

### 1. FSDS 增量 = 选对执行器，不是再加一个 loss

| 信号 | 直观 | 执行器 |
| --- | --- | --- |
| `c` 大、`δ` 静 | 世界长什么样变了，规则没变 | 缩 η / 重加权 / 记旧记忆 |
| `δ` 大 | 规则 / 目标变了 | 抬 η、重学映射、刷新奖励头 |
| 都静 | 别抖 | hold |

相对 `static`，七块 proxy 板 offline IW **全部正提升**（约 +4%～+21%，
`procgen` 从负变正）。增量来自“该动哪根杆”，不是新目标函数。

### 2. 奖励权重保真（合成回报）— 最硬的产品证据

真实标量回报常是合成的：`R = Σ w_k r_k`（progress / energy / safety …
或 dense progress / sparse gem）。漂移后旧 `w` 会把策略训歪。

- **保真指标**：`mix_mse = ‖Rc @ ŵ − R‖²`，`w_err = ‖ŵ − w_true‖`
- **FSDS 动作**：`δ` 大 → holdout 门控下把 `w` 朝新 batch 的经验 mix 推；
  仅 `c` 大 → 保留旧 mix，只对漂移头做轻抬权
- **数字**：bandit / maze / meta_vision 上 `mix_mse` 压 **~44×～116×**，
  `w_err` 压 **~5×～14×**

直观：**合成回报对准了，策略才在学“该学的目标”。**
policy IW 可以靠盲目加 LR 抬一点；mix 保真抬不上去 —— 这才是相对
`naive_boost` 的差异化。

### 3. 策略回报逻辑 — 适应速度 vs 目标是否正确

- **策略回报**（offline IW value）：策略在日志上拿到多少 —— 看“学得快不快”
- **合成回报保真**：训练信号是否还指向真实目标 —— 看“学得对不对”

`meta_only` 常在 policy IW 上最高（少步内环 + descriptor）；`reward_only`
在 policy 上中等，但把 mix 拉回真目标。`fsds_full` = 先校正 `w`，再用
typed η / descriptor 做内环 —— **先对目标，再加速**。

### 4. 基线保证：`static` vs `fsds_full`，合成回报与 `causal_proxy`

| 保证 | 含义 |
| --- | --- |
| vs `static` | 七板 policy 全升；有奖励组件的板上 mix 大幅更贴 `R` |
| 合成回报 | `fsds_full` 与 `reward_only` 共享刷新后的 `ŵ`；bandit/maze/meta_vision 上 mix_mse 从 O(0.3–7) 落到 O(0.01–0.1) |
| `causal_proxy` | 质量/摩擦 cov + 中段目标翻到 safety：`fsds_full` 是 **唯一 best module**（−0.439 → −0.363）；这里单靠 meta 不够，必须联用刷新 |
| `causal` 的边界 | 组件相关时经验 mix 易过拟合 → **holdout 门控**；未过门则 gate=0，退回旧 `w`（保守保证） |

一句话产品口径：**FSDS 增量 = 漂移分型选执行器；奖励保真保证合成回报对准；
策略回报保证相对无控制基线有升；`causal_proxy` 证明目标翻转时必须 full 闭环。**

### 5. Meta-learner uplift 监控 + IV selection（本轮增量）

把仓库本体的 **meta-learner OOD VIMP**（batch-as-treatment PO-risk / uplift
distance）接到闭环的 concept 通道，再做 **Important-Variable selection**：

| 步骤 | 做什么 |
| --- | --- |
| uplift monitor | hop 上 `W=batch`，`Y=R` → crossfit PO → `R_PO=E[τ̂²]`、ATE-proxy=`E[po]` |
| IV select | VIMP ≥ q-quantile，再截 top-K（CFPerm 风格：只放过漂移驱动特征） |
| 执行器门控 | 无 IV 的 feature group → η 保持 `η0`；descriptor 只编码选中 IV |

方法名：`fsds_meta_iv`（对照 `fsds_full` 的 ridge-δ）。看的额外效果：
`meta_iv_vs_full`、选中 `n_iv`、`R_PO`、mix_mse 是否更稳。

```bash
python3 scripts/run_fsds_rl_closed_loop.py --compare-meta-iv
```

**3-seed 实测（`--boards all --compare-meta-iv`）**

| board | meta_iv vs full | vs static | n_iv | R_PO | 解读 |
| --- | --- | --- | --- | --- | --- |
| `causal_proxy` | **+0.021** | +0.094 | 6 | 0.004 | 目标翻转：IV 聚焦有额外收益，best module |
| `meta_vision` | +0.011 | +0.172 | 7 | 0.049 | 跨域 hop：持平偏正 |
| `procgen_proxy` | +0.001 | +0.065 | 7 | 0.001 | 布局 OOD + gem 翻转：持平 |
| `mm_hop_proxy` | +0.001 | **+0.313** | 8 | 0.026 | 多模态 V/A/T：大 lift，PO-risk 稳 |
| `recsys_proxy` | −0.011 | +0.166 | 7 | 0.025 | CTR 概念翻转：vs static 大，IV 略保守 |
| `trading_proxy` | −0.000 | +0.051 | 6 | 0.002 | 风险/收益翻转：稳定正 lift |
| `dialogue_proxy` | −0.005 | +0.147 | 7 | 0.028 | 安全 rubric：vs static 正，mix 压紧 |
| bandit_* / maze | −0.010～−0.031 | 仍大幅正 | 6–7 | — | ridge-δ 已够；IV 门控略保守 |

mix_mse 与 `fsds_full` **完全相同**（奖励刷新共享）—— IV 只动 η / descriptor，
不伤合成回报保真。额外产品效果：**在因果/目标翻转板上用 uplift+IV
比纯 ridge 检测更准；在简单 bandit 上不要过度门控；多模态 / RecSys /
对话板把 PO-risk 赋能扩到非机器人域。**

### 6. Diverse boards + PO-risk metrics（本轮）

四块新 proxy，覆盖推荐 / 金融 / 对话 / 多模态：

| proxy | stands in for | concept flip |
| --- | --- | --- |
| `recsys_proxy` | RecSys popularity cov | relevance → CTR |
| `trading_proxy` | finance vol / liquidity | return → risk-first |
| `dialogue_proxy` | dialogue domain vocab | helpfulness → safety rubric |
| `mm_hop_proxy` | MSR-VTT-like V/A/T | retrieval → caption / text |

新增 PO-risk 侧指标（全部进 suite table）：

| metric | 含义 |
| --- | --- |
| `R_PO` / `ATE` | uplift monitor $E[\hat\tau^2]$、$E[\mathrm{po}]$ |
| `ESS` | Kish ESS fraction of self-norm IW ↑ |
| `regret` | batch-mean $R$ − IW value ↓ |
| `IV Jaccard` | hop-to-hop IV 稳定性 ↑ |
| `IV entropy` | IV 质量在 feature groups 上的熵 |

```bash
python3 scripts/run_fsds_rl_closed_loop.py --diverse --compare-meta-iv
python3 scripts/run_fsds_rl_closed_loop.py --boards all --compare-meta-iv
```

## Opportunity matrix

| dataset | FSDS channel | best-fit module | why it is a closed-loop opportunity |
| --- | --- | --- | --- |
| Meta-World ML1/10/45 | task-param cov + sparse reward concept | meta_adapt + region_lr | task descriptor = top-K FSDS features → few-shot inner loop |
| Procgen (train→test levels) | layout/style cov; gem/progress concept | reward_refresh + region_lr | OOD levels are cov; sparse-reward games flip concept mid-curriculum |
| D4RL (umaze→medium/large) | state-support cov | region_lr (shrink) + IW | offline support mismatch; naive η↑ hurts |
| CausalWorld | mass/friction cov; objective concept | all three | causal params → descriptor; safety objective flip → reward mix |
| Meta-Dataset | domain/style cov; class map concept | meta_adapt | cross-domain few-shot; descriptor conditions init |
| MEWA | collab-mode concept | reward_refresh | partner/protocol change = reward-head staleness |
| RecSys / slate bandits | popularity cov; CTR concept | reward_refresh + meta_iv | engagement farming = preference flip |
| Trading / portfolio | vol-regime cov; risk objective | reward_refresh + region_lr | risk↔return flip mid-stream |
| Dialogue / LLM policy | domain vocab cov; safety rubric | reward_refresh + meta_iv | safety-first rubric as concept |
| MSR-VTT multimodal hop | video cov; text concept | meta_iv + region_lr | PO-risk IV on V/A/T groups |
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
| `recsys_proxy` | RecSys popularity + CTR concept flip |
| `trading_proxy` | finance vol/liquidity + risk↔return flip |
| `dialogue_proxy` | dialogue domain + safety rubric flip |
| `mm_hop_proxy` | multimodal V/A/T hop (MSR-VTT-like) |

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

Measured offline IW value + reward-mix MSE + PO-risk monitors on **11 boards**
× 3 seeds (`--boards all --compare-meta-iv`):

| signal | finding |
| --- | --- |
| vs `static` | every FSDS module lifts policy value; largest lifts on `bandit_concept` / `bandit_mixed` / `mm_hop_proxy` |
| reward refresh | cuts `mix_mse` by 10–100× on bandit / maze / meta_vision / recsys / dialogue when holdout accepts the empirical mix |
| `fsds_meta_iv` | best module on `causal_proxy` / `mm_hop_proxy` / `meta_vision` / `procgen_proxy`; PO-risk IV adds under objective / modality flips |
| diverse boards | RecSys / trading / dialogue / multimodal all show **positive** meta_iv vs static (+0.05～+0.31) |
| PO-risk metrics | `R_PO` / ATE / ESS / IV-Jaccard expose when IV is informative vs over-gated |
| `naive_boost` | still competitive on *policy* IW — TSS polarity shows up more cleanly on **mix_mse / w_err** than on policy value in these linear proxies |
| holdout gate | required: unconstrained Ridge refresh overfit on `causal_proxy` and raised mix_mse |

**Practical closed-loop rule:** use FSDS `(c,δ)` to *choose the actuator*
(meta inner-loop vs reward-mix re-fit vs region η), not to invent a new loss.
Policy-η boost alone is an incomplete story under pure covariate boards — pair
it with IW / memory when swapping in real D4RL / Procgen loaders.
