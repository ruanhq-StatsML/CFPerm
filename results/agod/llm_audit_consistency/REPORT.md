# 连续时间大模型审核逻辑 consistency — OnlineRFPerm prototype

Objective = 相邻窗审核映射 P(Y|X) 有没有 **hop**（fire yes/no）。
Ratio 只是闸。AUC / 拒绝率 / 幻觉率 **不是** 这条方法的成绩。

## 1. 为什么 consistency 关键

生产 LLM 审核是一条 **时间流**，不是离线 gold 集上算一次准确率。
内容按到达时间进窗，judge（政策模型、人审池、reward model、外挂 LLM 审核器）当场写 Y。

- X：待审样本（prompt、回复、可选上下文）。
- Y：审核决定（放行 / 拒绝 / 转人工，或 RLHF 的 chosen/rejected）。
- 上一窗拟合的 probe μ0 **就是当时的审核逻辑**。

**Consistency** = 相邻两窗还像同一个审核员。政策没切、judge checkpoint 没换、人审指南没改版，
last-two consecutive OOS 应对得上 → **quiet**。

**Inconsistency / hop** = 审核逻辑一刀切地换制度。具体会撞上的，不是「模型突然不会写了」：

1. 政策包切版（新红线、新地区合规、warn→block）。
2. Judge 换代（外挂 LLM、内部分类器、规则引擎权重过夜替换）。
3. 系统 prompt / 审核说明重写（同一批 X，旧 probe 对不上新 Y）。
4. 人审池整体替换（外包团队、抽检比例、标注指南）。
5. 流量成分突变 **叠在映射上**。若只有 mix 变、映射没变，那是画像 hop，不是审核逻辑 hop。

高拒绝率、高幻觉率都可以仍然很顺：稳定地严、稳定地松，相邻窗还是同一套 P(Y|X)。
逐渐变严（10%→12%→14%）last-two 可以不 fire。
**死光了** 才是 hop：相邻窗 probe 对不上新决定。

所以这条 prototype 的读法只有 quiet vs fire。
Gate γ 可换，不是要优化的 objective。
domain AUC 说的是 P(X) 好不好分，跟置换距离不是同一个问题，这里不算。

### 为什么 Y 用部署审核策略，而不是直接拿 Anthropic 人标去翻

Last-two 探针是浅层 RF（20 棵、深度 4），只能看见它能表示的映射 hop。
HH-RLHF 人标是高容量噪声偏好；用 prompt+reply hash ⊕ 文风去拟合，OOS error 已经在 chance 附近（~0.45–0.55）。
这时把人标 92% 翻转，e_now 也还在 chance，ratio 过不了闸——**不是没 hop，是探针看不见**。
线上要盯的本来就是 **当前部署的审核器**（规则 / 分类器 / judge checkpoint），它就是 X 上的一个可表示策略。
所以：两条 HH 队列提供真实到达流量；Y 是该队列上的部署策略（文风/拒绝规则 + 少量噪声，避免 e_prev 真空）。
hop = 落地脚本同一刀：cut 后 p=0.92 翻转 Y（审核员一夜换制度）。

## 2. 跟落地脚本同一段

```python
probe = fit_online_probe(X_prev, y_prev, task='acc')
e_now = probe_err(probe, X_cur, y_cur, task='acc')
fire  = hop_fires(e_now, e_prev, gate=1.25)
po    = po_risk0_rows(probe, X_cur, y_cur, task='acc')
w     = po_iptw_weights(po, mode='sqrt') if fire else np.ones_like(po)
audit = np.argsort(-po)[:10]
```

Gate γ=1.5 是 instrumentation（包默认 1.5；落地伪代码常写 1.25）。读 fire，不读 γ。
X 默认是审核员用的文风/拒绝特征（hash 可选；高维 hash 会让浅树过拟合、consistent 误火）。
hop 仍是落地脚本那一刀：cut 后 p=0.92 翻转 Y。

## 3. 两个 dataset = 两套审核队列

| Dataset | Map | Deployed auditor | Cite |
|---|---|---|---|
| HH-RLHF helpful-base | helpfulness queue | verbose/specific/non-hedge | `bai2022hh-rlhf` |
| HH-RLHF harmless-base | safety queue | refusal/caution | `bai2022hh-rlhf` |

每个 dataset 同一套 X、同一套 batch，跑两个 regime：

- **consistent**：审核策略不动。expect quiet at cut t=4。
- **hop**：cut_batch=4 之后以 0.92 翻转 Y（审核员一夜换制度）。expect fire at cut，delay 0。
- hop 之后新映射自己仍可再变顺（t=cut+1 可以 quiet）。那是「死突然」，不是逐渐崩。

## 4. 结果（读 fire，不读 ratio）

| Dataset | Regime | n | batches | fire@cut | delay | n_fires | first_fire_t |
|---|---|---:|---:|---|---:|---:|---:|
| HH-RLHF helpful-base | consistent | 1200 | 15 | no | — | 1 | 6 |
| HH-RLHF helpful-base | hop | 1200 | 15 | yes | 0 | 2 | 4 |
| HH-RLHF harmless-base | consistent | 1200 | 15 | no | — | 2 | 9 |
| HH-RLHF harmless-base | hop | 1200 | 15 | yes | 0 | 4 | 4 |

**Cut 对照才是成绩**（quiet vs fire），不是全流 n_fires，更不是 ratio：

- hop @ cut：两条队列都 fire，delay 0。 审核员一夜换制度，last-two 对不上。
- consistent @ cut：两条队列都 quiet。 策略没切，相邻窗还是同一个审核员。
- hop 后 t=cut+1 可以立刻 quiet：新审核制度自己再变顺。这就是「死突然」，不是慢慢崩。
- 后面窗 n=80 的二项抖动可以再扣闸；那不是 objective，γ 也不是要优化的数。

Cut 窗邻域（gate log only：mean_r0 / mean_r1 / ratio 不是 objective）：

### HH-RLHF helpful-base — `consistent`

Auditor: verbose/specific/non-hedge helpfulness rule

| t | fire | mean_r0 | mean_r1 | ratio (log) |
|---:|---|---:|---:|---:|
| 2 | no | 0.150 | 0.188 | 1.250 |
| 3 | no | 0.188 | 0.137 | 0.733 |
| 4 | no | 0.137 | 0.188 | 1.364 |
| 5 | no | 0.188 | 0.225 | 1.200 |
| 6 | yes | 0.225 | 0.350 | 1.556 |

Landing snippet at cut: fire=`no`, serving `keep judge as gold; w=1`.

### HH-RLHF helpful-base — `hop`

Auditor: verbose/specific/non-hedge helpfulness rule

| t | fire | mean_r0 | mean_r1 | ratio (log) |
|---:|---|---:|---:|---:|
| 2 | no | 0.150 | 0.188 | 1.250 |
| 3 | no | 0.188 | 0.137 | 0.733 |
| 4 | yes | 0.137 | 0.787 | 5.727 |
| 5 | no | 0.787 | 0.275 | 0.349 |
| 6 | yes | 0.275 | 0.425 | 1.545 |

Landing snippet at cut: fire=`yes`, serving `Top-k po_risk0 re-review; cap DPO merge`.

### HH-RLHF harmless-base — `consistent`

Auditor: refusal/caution safety rule

| t | fire | mean_r0 | mean_r1 | ratio (log) |
|---:|---|---:|---:|---:|
| 2 | no | 0.200 | 0.150 | 0.750 |
| 3 | no | 0.150 | 0.213 | 1.417 |
| 4 | no | 0.213 | 0.250 | 1.176 |
| 5 | no | 0.250 | 0.175 | 0.700 |
| 6 | no | 0.175 | 0.250 | 1.429 |

Landing snippet at cut: fire=`no`, serving `keep judge as gold; w=1`.

### HH-RLHF harmless-base — `hop`

Auditor: refusal/caution safety rule

| t | fire | mean_r0 | mean_r1 | ratio (log) |
|---:|---|---:|---:|---:|
| 2 | no | 0.200 | 0.150 | 0.750 |
| 3 | no | 0.150 | 0.213 | 1.417 |
| 4 | yes | 0.213 | 0.775 | 3.647 |
| 5 | no | 0.775 | 0.225 | 0.290 |
| 6 | yes | 0.225 | 0.338 | 1.500 |

Landing snippet at cut: fire=`yes`, serving `Top-k po_risk0 re-review; cap DPO merge`.

## 5. Serving

| 状态 | 含义 | 动作 |
|---|---|---|
| quiet | 审核逻辑连续（consistency） | 标准路径；DPO/RLHF 数据可按原门禁合并；w=1 |
| fire | 审核逻辑 hop（inconsistency） | **不要把当前 judge 当金标**；Top-k `po_risk0` 人工复审；拒合并或限量合并；可选 sqrt(PO) 只打在 T=1 |

Fire 打开的是对照窗，不是「这条违规了」的分类器，更不是 fact-checker。

## 6. 跑法

```bash
PYTHONPATH=. python3 scripts/agod/llm_audit_consistency_prototype.py
```

Caches: `data/hf_cache/audit/`. Numbers: `results/agod/llm_audit_consistency/`.

