# OnlineRFPerm → LLM 推理：Justify + Low-hanging-fruit Roadmap

> **一句话**：OnlineRFPerm 不是「再训一个 LLM 监控模型」，而是把已有 **在线制度探针** 挂到助手/智能体 **连续推理日志流** 上——给「答变差了」一个 **可时间戳的 fire**，再开对照窗、路由动作。  
> 这是 bulletin 三条里 **最低挂果** 的一条：包已在、原型已跑、动作接口清晰，差的是生产接线与标签源硬化。

相关：`LLM_STREAM_CHANGEPOINT.md` · `LLM_INFER_ALIGN_PROTOTYPE.md` · `METHOD_LANDING_ROADMAP.md` · `LLM_LANDING_USECASES_BIZ.md` §U1 · `LANDING_BULLETIN_POLISH.md`

---

## 0. 为什么这条是 low-hanging fruit

| 条件 | 现状 | 为什么算「低挂」 |
|------|------|------------------|
| **方法已成包** | `agod.online_rfperm`：`fit_online_probe` / `probe_err` / `hop_fires` / `po_risk0_rows` | 不写新算法，只接线 |
| **LLM 侧 form 已定** | `{t, question, answer, embedding, score}` | HaluEval / HH 已跑通 |
| **硬证据已有** | Halu 流 detection delay = **0 batch**；ratio@fire ≈ 30；infer hop fired | 「能火」不是空话 |
| **动作字典已有** | `model_rollback` / `retrieval_refresh` / `audit_topk` | fire → 工单，不发明新产品 |
| **与 PSI/离线 AUROC 差分清晰** | 盯 \(P(Y\mid X)\) 代理制度 hop，不是边缘分布 | 对外一句话讲得通 |
| **不依赖客服 ¥ 账才能卖** | Phase-0 只交 delay + fire 时间戳 + 动作路由 | CS 账可后接；暂停 CS 也不挡这条 |

**不是 low-hanging 的部分（刻意后置）**：全因果根因、自建事实核查器、无标签纯无监督「语义质量」、整套客服财务闭环。那些不是这条路线的入场券。

---

## 1. Justify（对外怎么讲、对内怎么钉）

### 1.1 业务问题（只答这一句）

> 助手/智能体某一窗起答非所问、幻觉、转人工暴增——**是不是质量制度 hop 了？从哪一窗开始算「出事」？**

没有时间戳 → 无法开对照窗 → 回滚/刷检索/抽审全靠体感。  
OnlineRFPerm 的产品价值就是：**把「出事了」钉成事件**。

### 1.2 方法映射（一张表）

| 业务对象 | 方法对象 | 接口 |
|----------|----------|------|
| 连续 serving / agent 日志 | 时间流 batch \(B_t\) | `Stream` / 日切 / 小时窗 |
| 「这句好不好」的弱标签 | \(Y\)：幻觉标 / 人工纠错 / 裁判分 / 转人工 | `score` 列 |
| 问题·回答表征 | \(X\)：hash ⊕ style / embedding / 路径特征 | `embedding` |
| 「制度变了吗」 | 上一窗拟合探针 μ₀，本窗 OOS 误差比 \(e_{\mathrm{now}}/e_{\mathrm{prev}}\) | `hop_fires(..., gate)` |
| 「先审谁」 | fire 后 `po_risk0` Top-k | `po_risk0_rows` → `audit_topk` |
| 「动哪条臂」 | rag_hit 低 vs 制度火且 rag 够 | 路由轴（§1.4） |

核心公式（落地抄这段）：

```python
probe = fit_online_probe(X_prev, y_prev, task="acc")
e_now = probe_err(probe, X_cur, y_cur, task="acc")
fire  = hop_fires(e_now, e_prev, gate=1.25)   # 打开对照窗
po    = po_risk0_rows(probe, X_cur, y_cur, task="acc")
audit = np.argsort(-po)[:k]                    # 只在 fire 期升配
```

### 1.3 相对常规做法差在哪（justify 核心）

| 常规 | 缺口 | OnlineRFPerm 补什么 |
|------|------|---------------------|
| PSI / KS / embedding 漂移 | 只看 \(P(X)\)，文风漂 ≠ 幻觉律漂 | 盯预测性制度 \(P(Y\mid X)\) 代理 |
| 离线 AUROC / 周报人工抽检 | **无时间戳**，开不了对照窗 | fire 事件 + detection delay |
| 静态置信度 / logprob 阈值 | 制度 hop 后校准失效 | 探针绑在 **上一窗**；hop 后重估 |
| 「再训一个质量分类器」 | 贵、慢、和 serving 解耦 | 浅层 RF 探针，OOS 比，挂现有特征 |

### 1.4 轴拆分（低挂但仍必须写死，否则误动作）

同一「答得差」拆三轴，**三张工单禁混账**：

| 轴 | 信号 | 正确动作 | 禁止 |
|----|------|----------|------|
| **生成制度** | OnlineRFPerm fire，rag 尚可 | `model_rollback` / canary / 限流 | 只刷检索假装修好 |
| **检索缺口** | Top-k `rag_hit` 低 | `retrieval_refresh` | 一刀切回滚整模 |
| **文风 / register** | style domain / FSDS covariate | 改模板 · decoding（**另账**） | 触发客服回滚或偏好合并 |

### 1.5 明确边界（防 overclaim）

- fire ≠「这句话事实错了」——不是事实核查器。  
- fire ≠ 因果根因——只打开窗与门禁。  
- P@10 ≠ 少退款数——少事故来自动作落地后的费率差（有账本时再接）。  
- 文风轴 **永不进** 客服主账。

### 1.6 已有证据（本仓库，可引用）

| 证据 | 落点 | 读法 |
|------|------|------|
| HaluEval 流 delay=0 | `results/agod/llm_changepoint/` | cut 当窗就火 |
| infer hop fired，ratio 大 | `LLM_INFER_ALIGN_PROTOTYPE.md` | 制度跳变可检测 |
| rag 路由 vs 回滚 | `hf_landing` / bulletin landing | 避免误回滚 |
| 包 API 稳定 | `agod.online_rfperm` | 生产可直接 import |

---

## 2. Roadmap（按可交付切片，不估日历）

原则：**先火情时间戳 → 再路由 → 再审计排序 → 最后才接 ¥**。  
每阶段都有 **DoD（完成定义）** 与 **不做清单**。

### Phase 0 — 接线最小闭环（**现在就能交**）✅ 原型已齐

| 项 | 内容 |
|----|------|
| **交付** | serving/合成日志 → `{t,q,a,X,y}` → OnlineRFPerm → `fired` + `ratio` + `delay` |
| **脚本** | `llm_stream_changepoint.py` · `llm_infer_align_prototype.py` · `hf_landing_protos.py` |
| **DoD** | 固定数据可复现 delay；文档可对外引用；quiet 不扰民 |
| **不做** | 不上 ¥ 账；不换真 embedding；不接工单系统 |

**下一刀（仍属 Phase 0 硬化）**

1. 生产日志 adapter：把现有助手 log 列映射到 form（字段字典一页纸）。  
2. Y 源优先级写死：`human_correction` > `judge_score` > `escalation` > 合成幻觉标。  
3. gate / e_floor 按业务窗长扫一版（报 quiet 误火率）。

### Phase 1 — 动作路由上线（**真正的 low-hanging 产品价值**）

| 项 | 内容 |
|----|------|
| **交付** | fire 后自动建议：`retrieval_refresh` vs `model_rollback` vs `audit_topk` |
| **依赖** | 同窗 `rag_hit`（或 overlap 代理）；路径特征可选 |
| **DoD** | 路由表 + 误路由用例（rag 低却回滚 / 制度火却只刷检索）各 ≥3 条对拍 |
| **KPI** | 误大动作次数↓（相对「一律回滚」基线）；人工确认通过率 |
| **不做** | 自动无审批回滚；文风信号进路由主表 |

路由伪码：

```text
if not fire:          quiet → 标准路径
elif rag_hit low:     retrieval_refresh
elif fire & rag ok:   model_rollback / canary
else:                 audit_topk only
```

### Phase 2 — 审计止血（fire 期人力）

| 项 | 内容 |
|----|------|
| **交付** | fire 窗内 `po_risk0` Top-k 进人工 / 二次生成 / 强制检索 |
| **DoD** | quiet 权重退火回均匀；审计件成本可记账（有账时进净贡献） |
| **KPI** | 同人力下捕获坏例率 vs 随机抽；或固定捕获率下审计件数↓ |
| **不做** | 把 P@10 讲成「少了 N 单退款」 |

### Phase 3 — 智能体路径子集（bulletin「连续推理路径」加深）

| 项 | 内容 |
|----|------|
| **交付** | 在 fire 后对 path 维做轻量 FSDS/LOGO：`step_depth` / `tool_call` / `retry` / `route` |
| **DoD** | 输出「先动哪段路径」；与整模回滚互斥纪律 |
| **同构** | 与图谱 L1→L2、多塔「只刷漂的一路」同一句式 |
| **不做** | 新建 agent 产品；ASSOC 边讲因果 |

### Phase 4 — 对照窗经营账（**可选**；CS 暂停则停在此处设计）

| 项 | 内容 |
|----|------|
| **交付** | fire 期内 acted vs ignored → 工单/退款/承接差 → 周归因 / 回本 |
| **前置** | Phase 0–1 稳定；有单价与量 |
| **DoD** | 周归因 gap≈0；AUROC 不进贡献句 |
| **现状** | 客服 SQL 迭代 **已暂停**；本 roadmap **不阻塞** 在 Phase 0–3 |

---

## 3. 优先级看板（接下来先做谁）

| 序 | 切片 | 为何先做 | 依赖 |
|----|------|----------|------|
| **R0** | 生产 log → form adapter + Y 源字典 | 没有这一步后面全是 demo | 日志只读权限 |
| **R1** | Phase 1 路由表 + rag_hit 同窗 | 直接省误回滚，业务体感最强 | R0 |
| **R2** | quiet 误火率 / gate 标定 | 值班敢接；否则会关系统 | R0 |
| **R3** | Phase 2 po_risk0 Top-k 接审计队列 | 人力稀缺时立刻值钱 | R1 |
| **R4** | Phase 3 path 特征进 X / 子集定位 | bulletin「智能体路径」兑现 | R1 |
| **R5** | Phase 4 ¥（若重启 CS） | 财务签字 | R1–R3 稳 |

**刻意不做（本 roadmap 范围外）**

- 用 OnlineRFPerm 替代事实核查 / RAG 评测套件  
- 无 Y 的纯无监督「语义质量分」当主火源  
- 把对齐合并门禁（批式 RFPerm）与推理流混账进同一张客服账  
- 等「完美 embedding」才上线——hash/现有表征够跑 Phase 0–1  

---

## 4. 评估怎么报（别用一个 AUROC 糊弄）

| 层 | 指标 | 及格读法 |
|----|------|----------|
| **检测** | detection delay；quiet 误火率；ratio@fire | delay 小；误火可控 |
| **路由** | 误回滚率；误只刷检索率 | 相对「一律回滚」↓ |
| **审计** | Top-k 捕获率 vs 随机；审计件成本 | 同成本更高捕获 |
| **路径子集** | 注入 path 族 Hit（有演练时） | 与图 L1 Hit 同构 |
| **经营（可选）** | act-vs-ignore 量 × 单价；周归因 gap | gap≈0 才对外讲周贡献 |

禁止：用「探针 AUROC 高」直接等于「贡献了多少 ¥」。

---

## 5. 对外 30 秒话术

> 我们把 OnlineRFPerm 接到 LLM **连续推理日志**：用上一窗浅层探针盯本窗质量制度有没有 hop。  
> **火了** → 打开对照窗，并按检索缺口 vs 生成制度路由到刷检索 / 回滚 / Top-k 审计——而不是 PSI 看长度、也不是一刀切回滚整模。  
> 包和原型已在仓库跑通（Halu 流 delay=0）；下一步是生产日志接线与路由表，不重训监控大模型。

---

## 6. 索引与命令

| 文档 / 脚本 | 用途 |
|-------------|------|
| `agod/online_rfperm.py` | 方法包 |
| `scripts/agod/llm_stream_changepoint.py` | 推理流 delay |
| `scripts/agod/llm_infer_align_prototype.py` | 最小 call 包 |
| `scripts/agod/hf_landing_protos.py` | 幻觉 + rag 路由 |
| `scripts/agod/bulletin_landing_attr.py` | bulletin 三场景板 |
| `METHOD_LANDING_ROADMAP.md` | 创新 A–E + 周归因（经营层） |
| `LLM_LANDING_PROTOTYPE_METHOD_LIST.md` | S1–S10 场景表 |

```bash
PYTHONPATH=. python3 scripts/agod/llm_stream_changepoint.py
PYTHONPATH=. python3 scripts/agod/llm_infer_align_prototype.py
PYTHONPATH=. python3 scripts/agod/hf_landing_protos.py
```

---

## 7. 一句收束

**Justify**：OnlineRFPerm 给 LLM 推理一个 **制度 hop 时间戳**，这是对照窗与正确动作臂的前提；相对 PSI/离线 AUROC，差的就是「能开窗」。  
**Roadmap**：接线 form → 路由表 → Top-k 审计 →（可选）路径子集 →（可选）¥；**前三步就是 low-hanging fruit**，不依赖客服账重启。
