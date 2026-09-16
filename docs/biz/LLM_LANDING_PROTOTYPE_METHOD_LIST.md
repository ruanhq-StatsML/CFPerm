# 大模型落地：Prototype 清单 × 方法论场景（带 Justify）

> **客服助手增量账迭代已暂停**（本文件不扩 CS SQL）。  
> 目标：把「能跑的 prototype」和「方法论场景」收成一张可对外 justify 的清单——  
> **不是**新 LLM/推荐算法，**是**检测 → 子集定位 → 动作/门禁 的生产 use-case。

---

## 0. 怎么读这张清单

| 列 | 含义 |
|----|------|
| **Prototype** | 仓库里可跑的脚本 / 产物（证据链） |
| **方法论场景** | 方法标签只作路由：OnlineRFPerm / RFPerm / PO·BOCPD / FSDS·LOGO / 轴拆分 |
| **Justify** | 为什么值得做：生产痛点 + 若不做会误伤什么 + 钱/KPI 落点 |

**总纪律**

1. 方法回答「火了吗 / 审谁 / 动哪一维」；账本回答「值多少」——别用 AUROC 当贡献。  
2. 文风轴 **另账**，不进客服主账。  
3. 多塔/三塔 = **子集定位形状**，不是做推荐系统。  
4. 风格归因与图谱归因 **先并成特征维度**，再 L1→L2。

**OnlineRFPerm → LLM 推理（low-hanging fruit）**：justify + 分阶段 roadmap 见  
`ONLINERFPERM_LLM_INFER_ROADMAP.md`（接线 form → 路由 → Top-k 审计；不依赖客服账）。  
**真推理流 use-case（transformers / vLLM）**：`ONLINERFPERM_LIVE_INFER.md` / `scripts/agod/online_rfperm_live_infer.py`。

---

## 1. 可跑 Prototype 清单（证据链）

| ID | Prototype | 跑法 | 主产物 | 覆盖的方法论场景 |
|----|-----------|------|--------|------------------|
| **P1** | HF 双面落地（幻觉 + Judge/文风） | `scripts/agod/hf_landing_protos.py` | `results/agod/hf_landing/` | S1 推理制度火；S2 对齐门禁；S3 文风轴；S5 审计排序；S6 RAG 路由 |
| **P2** | 推理/对齐流变点 | `scripts/agod/llm_stream_changepoint.py` | `results/agod/llm_changepoint/` | S1 OnlineRFPerm 在线；S2 对齐对照流；detection delay |
| **P3** | 推理/对齐最小 call 包 | `scripts/agod/llm_infer_align_prototype.py` | `results/agod/llm_infer_align/` + `LLM_INFER_ALIGN_PROTOTYPE.md` | S1+S2 最小接口；PO Top-k；√PO 门控 |
| **P4** | HH 在线 tidy 流 | `scripts/agod/hh_online_rfperm_stream.py` | `results/agod/hh_online_stream/` | S2 偏好流；S3 文风特征准备 |
| **P5** | Bulletin 三类落地板 | `scripts/agod/bulletin_landing_attr.py` | `results/agod/bulletin_landing/` | S1 智能体路径；S2 质检审核；S3 画风；文本/特征归因 |
| **P6** | 特征维度并维归因 | `scripts/agod/feature_dim_attr.py` | `results/agod/feature_dim_attr/` | S7 风格∪图谱→特征维；S8 多维 L1→L2 |
| **P6b** | **画风漂移迭代** | `scripts/agod/style_drift_iter_proto.py` | `results/agod/style_drift_iter/` + `STYLE_DRIFT_ITER_PROTOTYPE.md` | S3 检测→归因→style-only 工单→再检 |
| **P7** | 方法×业务情景对照 | `scripts/agod/method_biz_scenario_prototype.py` | `METHOD_BIZ_SCENARIO_PROTOTYPE.md` | 方法层↔情景层 relevance（含客服账指针，**本轮不迭代 CS**） |
| **P8** | 客服增量账 SQL（已冻结） | `scripts/agod/run_biz_value_sql_demo.py` | `docs/biz/CS_ASSISTANT_CONTRIBUTION.md` | S1 打开窗后的 ¥ 对账（暂停扩层） |

数据：`data/hf_cache/halueval_qa_3000.jsonl` · `hh_rlhf_helpful_base_2500.jsonl`。

---

## 2. 方法论场景清单（S1–S10）+ Justify

### S1 — 智能体 / 助手 **连续推理路径**制度火

| | |
|--|--|
| **方法** | OnlineRFPerm（在线误差比 fire） |
| **Prototype** | P1 / P2 / P3 / P5 |
| **业务问题** | 某一窗起答非所问、幻觉、转人工暴增——是不是 **标签律/质量制度 hop**？ |
| **Justify** | PSI/KS 只看边缘；离线 AUROC 无时间戳。需要把「出事了」钉成可开对照窗的事件，才能谈回滚/刷检索/抽审，否则运营只能凭体感。 |
| **动作** | 开对照窗 → `model_rollback` / `retrieval_refresh` / `audit_topk` |
| **KPI** | 工单 · 退款 · 承接（有账本时）；无账本时至少 delay / fire 时间戳 |
| **不是** | 事实核查器；fire ≠「这句话错了」 |

### S2 — **人工/合成数据**质量 · 对齐合并门禁

| | |
|--|--|
| **方法** | RFPerm / RFPerm-as-Judge（批式；可挂在线连更） |
| **Prototype** | P1 / P3 / P4 / P5 |
| **业务问题** | 这批偏好/合成语料能不能并进生产？裁判映射还跟黄金集同一套吗？ |
| **Justify** | 坏偏好进 DPO/RLHF = 长期防损事故。合并前需要「制度是否 hop」的门禁，而不是只看条数/长度。文风变了 ≠ 偏好变了，必须拆轴。 |
| **动作** | 拒合并 / 限量 / Top-k 人工复核偏好对 |
| **KPI** | 防损（坏批未进生产）；**不进客服工单账** |
| **不是** | 道德裁判；不是唯一因果归因 |

### S3 — **画风 / 文风 / register** 漂移（另账）

| | |
|--|--|
| **方法** | style domain AUC / FSDS 协变量侧 / 文本指标 LOGO |
| **Prototype** | P1 / P5 / P6 / **P6b（迭代闭环）** |
| **业务问题** | 说话味道、素材 register 变了吗？修完再检有没有熄火？ |
| **Justify** | 用户体感强，但钱路径是 CTR/品牌/素材重做。若与幻觉/偏好混账，会误伤偏好头或误触发客服回滚。 |
| **动作** | 改模板 · decoding · 创意批次 |
| **KPI** | CTR · 品牌（**硬规则：永不并进客服主账**） |

### S4 — **检索缺口 vs 生成制度**（动作路由轴）

| | |
|--|--|
| **方法** | rag_hit / Top-k rag 均值 + OnlineRFPerm（生成侧） |
| **Prototype** | P1 / P5（agent path） |
| **业务问题** | 答得差，是知识库撑不住，还是生成制度坏了？ |
| **Justify** | 误回滚整模成本高；误只刷检索会「假装修好」。路由轴直接省一次错误大动作。 |
| **动作** | rag 低 → 只刷检索；制度火且 rag 够 → 回滚/canary 生成侧 |

### S5 — **审计止血**（谁先审）

| | |
|--|--|
| **方法** | PO-risk / BOCPD 实例排序（本仓库 PO 接口） |
| **Prototype** | P1 / P3 |
| **业务问题** | 火情期人力有限，先审哪几条最值？ |
| **Justify** | 随机抽浪费审计件成本；排序把有限人力砸在高风险实例，审计单价进入净贡献扣减才站得住。 |
| **动作** | `audit_topk`；quiet 回均匀 |
| **不是** | 事实核查；P@10 ≠ 少退款数本身 |

### S6 — **Router / 专家 / 工具链** 哪一路掉了

| | |
|--|--|
| **方法** | 任务桶 / 专家路 domain 或 FSDS 路份额 → 子集刷新 |
| **Prototype** | P1（task_shift）；P5（path 指标）；P6（agent_path 维） |
| **业务问题** | 哪类题、哪路专家、哪段工具调用链质量在掉？ |
| **Justify** | 全量刷新专家/工具昂贵；子集定位只动掉点的那一路，与多塔逻辑同构。 |
| **动作** | 只刷新该专家 / 该检索器 / 限步深·审 tool 链 |

### S7 — **特征维度统一归因**（风格 ∪ 图谱先并维）

| | |
|--|--|
| **方法** | 特征维度目录 → L1 Top 维 → L2 FSDS/LOGO |
| **Prototype** | P6；文档 `FEATURE_DIM_UNIFIED_ATTR.md` |
| **业务问题** | 风格指标和图节点是不是两套系统？ |
| **Justify** | 两套归因对不上工单语言。先收成同一层 `D={d1..dk}`（文风=Author/style_register 维），告警与动作统一。 |
| **动作** | 只砸 Top-2/3 维内的刷新/审核 |

### S8 — **商户·作者·商品·订单（+用户）图** 两级定位

| | |
|--|--|
| **方法** | 图 localization（缩维）→ 维内 FSDS/LOGO |
| **Prototype** | P6（`--alias graph`）；`MERCHANT_AUTHOR_PRODUCT_ORDER_GRAPH_ATTR.md`；**评估协议 `GRAPH_USECASE_EVAL.md`** |
| **业务问题** | 近实时先报哪 2–3 个业务维有问题，再维内归因？ |
| **Justify** | 全维细归因噪且贵；L1 缩到 2–3 维再 L2，才是可上线形态。与多塔/特征族同一句式。主考题是 **L1 Hit@K**，不是边故事。 |
| **动作** | 工单落到商户/作者/货盘/交易运营（只砸 S*） |

### S9 — **多塔 / 多路刷新预算**（子集定位，不是推荐）

| | |
|--|--|
| **方法** | FSDS 塔/路 mass → 只刷漂的一路 |
| **Prototype** | 文档 `SUBSET_LOCALIZATION_AUDIT_TOWERS.md`；P6 同构 |
| **业务问题** | 视·音·文（或任意多路）哪一路漂了？ |
| **Justify** | 全塔重算浪费算力；方法不管产品叫推荐还是多模态——形状相同。 |
| **动作** | 只失效/重编码/LoRA 该路 |

### S10 — **安全拒答 / 周更审核 / Embedding 批次**

| | |
|--|--|
| **方法** | RFPerm 批门禁；安全工作点 hop；几何预检 + FSDS |
| **Prototype** | 文档 `WEEKLY_AUDIT_SAFETY_GRAPH_FSDS.md`；对齐侧 P1/P3 |
| **业务问题** | 周更偏好能否放行？拒答策略工作点漂了没？新 embedding 批能否上？ |
| **Justify** | 周更是高频防损；安全要拆漏拒/误拒；坏向量批污染下游全体。 |
| **动作** | quiet 放行 / fire Top-k；策略双跑；拒坏批次 |

---

## 3. Bulletin 三句 ↔ 场景映射（对外主线）

打磨稿见 `LANDING_BULLETIN_POLISH.md`。

| Bulletin 领域 | 场景 | Prototype | Justify 一句 |
|---------------|------|-----------|--------------|
| AI 智能体连续推理路径 | S1 + S4 + S6 | P2 / P5 | 制度火 + 路由轴，避免整模误回滚 |
| 人工/合成数据质检与审核 | S2 + S5 + S10 | P1 / P3 / P5 | 合并门禁 + Top-k，防损优先 |
| 画风/文风漂移 | S3 + S7 | P1 / P5 / P6 / **P6b** | 另账；检测→归因→修→再检 |
| 「归因到文本指标与特征」 | S7 + S8 | P5 / P6 | 先维后特征，动作可落地 |

---

## 4. 方法论创新点 ↔ 场景（justify「方法为什么要存在」）

| 创新点 | 相对常规 | 主要服务场景 | Justify |
|--------|----------|--------------|---------|
| **A OnlineRFPerm** | vs PSI/离线 AUROC | S1, S2 | 给制度 hop **时间戳**，才能开对照窗与门禁 |
| **B 轴拆分** | vs 单一 alignment score | S3, S4, S2 | 一文风/一偏好/一RAG → 三张工单，禁混账 |
| **C PO-risk 排序** | vs 随机抽审 | S5 | 审计人力稀缺时的优先级；进净成本 |
| **D FSDS/LOGO 子集** | vs 全空间 explainer | S7–S9 | 只刷新漂的维/塔/族，省算力少误伤 |
| **E 图→维→特征** | vs 直接全特征 FSDS | S8 | 先缩维再钻，近实时可上线 |

价值句式（不变）：**检测打开窗 → 子集定位 → 动作/门禁**；（有经营面时再）**act-vs-ignore 量 × 单价**。

---

## 5. 优先级建议（下一轮做 LLM，不扩 CS）

| 优先级 | 做什么 | 为什么 |
|--------|--------|--------|
| **P0** | 保持 P1+P5+P6+**P6b** 可复现，对外用 §3 bulletin 表 | 三类落地 + 并维 + **画风迭代**已齐 |
| **P0** | **OnlineRFPerm→LLM 推理 roadmap R0–R1**：生产 log adapter + 路由表 | low-hanging；详见 `ONLINERFPERM_LLM_INFER_ROADMAP.md` |
| **P0** | **真推理流 use-case** `online_rfperm_live_infer.py`（transformers/vLLM） | HaluEval 真生成 → fire delay=0；见 `ONLINERFPERM_LIVE_INFER.md` |
| **P1** | 加强 S4/S6：真实 embedding 替换 hash；path/tool 指标接生产日志 | 路由轴与智能体路径是误动作成本最高处 |
| **P1** | S8 用一窗真实商户·作者·商品·订单切片跑 P6 `--alias graph` | 把「能实时」从示意变成业务维证据 |
| **P2** | S10 周更/安全双跑最小门禁板（批式 RFPerm + 漏拒/误拒分列） | 防损高频，但不依赖客服账 |
| **暂停** | 客服助手 SQL 再挖层 | 按你的要求先停；P8 保持现状数字即可 |

---

## 6. 索引

| 文档 | 用途 |
|------|------|
| `STYLE_DRIFT_ITER_PROTOTYPE.md` | 画风漂移迭代 Before→After |
| `LLM_LANDING_USECASES_BIZ.md` | U1–U10 业务总图 |
| `LANDING_BULLETIN_POLISH.md` | 对外 bulletin 打磨稿 |
| `FEATURE_DIM_UNIFIED_ATTR.md` | 风格∪图谱→特征维 |
| `SUBSET_LOCALIZATION_AUDIT_TOWERS.md` | 子集定位 / 多塔 / 审核 |
| `WEEKLY_AUDIT_SAFETY_GRAPH_FSDS.md` | 周更·安全·图→FSDS |
| `MERCHANT_AUTHOR_PRODUCT_ORDER_GRAPH_ATTR.md` | 四维图两级归因 |
| `GRAPH_USECASE_EVAL.md` | **图谱 use-case：刻画 + 评估（L1 Hit@K / L3 子集动作）** |
| `METHOD_LANDING_ROADMAP.md` | 创新 A–E + 落地接口 |
| `ONLINERFPERM_LLM_INFER_ROADMAP.md` | **OnlineRFPerm→LLM 推理 justify + low-hanging roadmap** |
| `METHOD_BIZ_SCENARIO_PROTOTYPE.md` | 方法×情景 relevance |
| `LLM_INFER_ALIGN_PROTOTYPE.md` | 最小 call 包结果 |

```bash
PYTHONPATH=. python3 scripts/agod/hf_landing_protos.py
PYTHONPATH=. python3 scripts/agod/bulletin_landing_attr.py --synth
PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py --synth
PYTHONPATH=. python3 scripts/agod/feature_dim_attr.py --alias graph --synth
PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py          # 画风迭代
PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py --synth
```
