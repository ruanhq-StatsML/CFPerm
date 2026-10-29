# 落地场景 × 可安排的数据

九条场景不要揉成一个「漂了就归因」。先定车道，再定 **一张预测表**：一列 Y，一堆 X。

| 车道 | 问的问题 | 方法 | 数据形态 |
|---|---|---|---|
| **A 映射监控** | P(Y\|X) 还是不是同一套（consistency / hop） | OnlineRFPerm last-two fire | 流：`(X, y, batch)` |
| **B 画像** | P(X) 变了没（文风 / 画风 / mix） | 分块 last-two / 域对照 | 只有 X，或 X 的塔 |
| **对照** | 哪座塔 / 哪族特征的画像先不像了 | 同一时钟，切块并排 | X 分块；**不是归因** |

A 不是归因。塔/族是画像对照。B 不能替代 A 的 fire。

---

## 总表：每条场景用什么数据

| # | 场景 | 车道 | Y（一个） | X（一堆） | 现在就能安排 | 时钟 |
|---|---|---|---|---|---|---|
| 1 | 机器学习模型上线监控 | A | 线上标签：点击、违约、转化、纠错 | 当时 serving 特征 | Amazon Reviews、Covertype、Criteo 子集、你们自己的 Tencent-GR | 按日 / 按 batch |
| 2 | 大模型推理监控 | A | 幻觉 / 忠实 / 人工纠错 / rag_hit 二值 | 问题+回答(+知识) 特征 | **HaluEval**、TruthfulQA、HotpotQA、SQuAD、RAGTruth | 到达窗 |
| 3 | 大模型自动数据质量检测 | A（批式门禁） | 质检通过 / 该不该进训练 | 文本质量特征、来源、长度、重复 | FineWeb-Edu 分、C4 vs Wiki、HH 候选批 vs 金标批 | 合并批 T=0/T=1 |
| 4 | 大模型审核 | A | 审核决定：过 / 不过 | 待审样本特征 | **HH helpful / harmless**、BeaverTails、WildGuard、ToxicChat | 审核到达窗 |
| 5 | 画风检测 | **B** | 通常 **没有** 业务 Y；要预测时 Y=风格类 | 图像/文案风格特征 | WikiArt、DiffusionDB、HH 文风向量、BAM | 素材批次 |
| 6 | 多模态连续学习 | A | 该任务标签：检索命中、匹配、caption 对错 | 各塔表征拼起来 | **MSR-VTT**、Fashion-IQ、COCO caption、WebVid | 任务/时间切窗 |
| 7 | 多模态塔对照 | **对照** | 无（或另有任务 Y 只给 A） | X 分塔并排 | 与 6 同一批，切开视频/文本 | 同窗 |
| 8 | 特征族 / 流量切片对照 | **对照** | 无 | 按族切开的 X | 审核 `x_*` 族、Tencent-GR 族 | 日 / 窗 |
| 9 | 模型迭代更新 | A | 新模型输出是否还服从旧映射；或「可合并」 | 同一批 X，新旧预测或偏好 | 同 prompt 两版本输出、HH 合并门禁、AlpacaEval 对 | 版本切面 |

下面按条展开：Y 是什么、X 从哪来、公开数据怎么切成流、不要把哪条误接到哪条车道。

---

## 1. 机器学习模型上线监控 — 车道 A

业务：模型已经在线上打分。要盯的是 **P(Y|X) 还跟发布时一样吗**，不是 PSI 看边缘，也不是事后归因「哪个特征害了 KPI」。

**一张表**

- Y：当时就能拿到、或短滞后拿到的标签。点击、下单、违约、欺诈、人工纠错。一个数。
- X：请求打到模型时的特征，必须是 serving 当时那一版，不能用未来才算好的衍生。
- batch：日切或流量窗。

**可安排的数据**

| 数据 | 怎么变成 (Y, X, batch) | 备注 |
|---|---|---|
| Amazon Reviews（按 `review_time`） | Y=是否 ≥4 星（或 helpful）；X=文本 hash/长度/类目 | 有自然时间，适合流 |
| Covertype / Airlines delay | Y=类别；X=表特征；人为切窗或注入 cut 后标签律翻转 | 经典 ML 流，验证 OnlineRFPerm 本身 |
| Criteo / Avazu 子集 | Y=click；X=稀疏特征 | 近广告/推荐 serving |
| 仓库已有 Tencent-GR 路 | Y=转化/点击；X=序列+统计特征 | 真推荐形态，优先于玩具表 |

**不要**：把「CTR 掉了所以归因到特征 17」写进这条。那是第 8 条。这条只问 fire 不 fire：映射 hop 了没有。quiet → 继续当金标；fire → 开对照窗、限流、回滚。

---

## 2. 大模型推理监控 — 车道 A

业务：助手每天在答。突然胡编 / 答非所问，是 **生成制度 hop**，不是单条 badcase，也不是幻觉率高低本身。

**一张表**

- Y：这条回答「坏没坏」。幻觉是/否、人工纠错、忠实度、或 rag_hit 过阈。一个二值（或分箱后仍当分类）。
- X：问题、回答、可选 knowledge 的特征。hash ⊕ 文风 ⊕ rag overlap。
- batch：serving 到达顺序。

**可安排的数据**

| 数据 | 状态 | Y | 切流 |
|---|---|---|---|
| `pminervini/HaluEval` `qa_samples` | **已有** landing + infer-bench | hallucination yes/no | 前半忠实、cut 后幻觉质量（制度 hop） |
| TruthfulQA / HotpotQA / SQuAD n=300 | **已有** `data/hf_cache/infer_bench_export/` | 对/错、可答/不可答 | 弱流；适合表征，切窗要人为 pack |
| RAGTruth | 可拉 | 幻觉跨任务标签 | 同时带 retrieval 轴，和生成 hop 拆开 |
| LMSYS-Chat-1M / WildChat | 可拉 | 需弱标：拒答、毒性、长度异常当 proxy Y | 最像真实 mix；Y 要自己定 |

HaluEval 已经证明：cut 上 fire = 幻觉 **标签律** hop。幻觉率从 8% 到 94% 只是旁证，objective 仍是 fire。

推理侧 X 可以和审核侧不同，Y 必须是「这次推理好不好」，不是审核过不过。

---

## 3. 大模型自动数据质量检测 — 车道 A，批式

业务：这批语料 / 偏好对 **能不能进训练**。对照的是黄金参考批，不是线上 last-two（也可以 last-two 扫连续采集日）。

**一张表**

- Y：质检员或规则给的「可进 / 不可进」，或相对金标的 chosen。一个数。
- X：文档/对话特征（重复率、长度、教育分、来源、语言）。
- 两个 batch：T=0 金标，T=1 候选合并。

**可安排的数据**

| 数据 | Y | 用法 |
|---|---|---|
| FineWeb-Edu（或 DCLM 过滤分） | 高/低教育分过阈 | 新爬取批 vs 已过检批 |
| C4 / The Pile 切片 vs Wikipedia | 来源当弱 Y，或人工抽检 Y | 脏批 vs 干净批门禁 |
| Cosmopedia 合成 vs 真人 Wiki | 合成=1 | 合成污染：映射是否把制度带偏 |
| HH 候选 vs helpful-base 金标 | 人标 chosen | **已有** 对齐门禁形态 |

和「推理监控」的差别：这里 fire → **拒合并 / 限量**，不是回滚线上模型。钱路径是防损，不是客服工单。

---

## 4. 大模型审核 — 车道 A（当前 prototype）

业务：审核员（人 / 政策模型 / 外挂 LLM）当场写决定。盯 P(Y|X) 还是不是同一套审核逻辑。

**一张表（已经落地）**

- Y：**一个**审核决定，过=1 / 不过=0。不是 13 维归因。
- X：待审回复上的一堆特征 `x_n_toks … x_refuse`。
- 文件：`results/agod/llm_audit_consistency/xy_*.csv`
- 流量来源：`data/hf_cache/audit/helpful_base_600.jsonl`、`harmless_base_600.jsonl`

HH 的 chosen/rejected **不当 Y**（浅探针预测不了）。chosen/rejected 只说明流量从哪条队列来。

**还可以安排、Y 更像真审核的数据**

| 数据 | Y | 为什么值得补 |
|---|---|---|
| BeaverTails / PKU-SafeRLHF | 安全类别或 is_safe | 真安全审核标签，不是模拟政策 |
| WildGuardMix / Aegis | 违规 / 安全 | 更近线上审核器 |
| ToxicChat / OpenAI Moderation eval | 毒性/违规 | 审核队列，X=对话特征 |
| CivilComments | toxicity | 长流，可按时间切 |

补真 Y 之后，仍然是预测：probe 学 P(Y|X)。hop = 政策包/judge 换代。不要把 13 个 x 当成「审核逻辑的归因维」。

---

## 5. 画风检测 — 车道 B（画像），偶尔才有 Y

业务：素材味道变了没有。CTR/品牌账，**不要**并进客服主账，也 **不要** 当成对齐失败。

默认没有审核那种 Y。对象是 **P(X)**：颜色、笔触、构图、文案正式度。

| 数据 | 当 X 还是当 Y | 用法 |
|---|---|---|
| WikiArt | X=图表征；可选 Y=artist/style 只为了可分性旁证 | 画风域移：早期流派 vs 后期 |
| DiffusionDB / Midjourney dumps | X=生成图+prompt 特征 | 生成画风 hop；Y 可缺省 |
| BAM / Danbooru 标签 | 标签可当弱 Y | 风格分类监控，仍是 B 不是审核 |
| HH 文风 13 维 | 已有，**这是 X 的风格子集** | 文风 AUC 旁证；不能替代审核 fire |

若业务硬要预测「这张图是否品牌合规」，那才进入车道 A：Y=合规，X=视觉特征。那是审核的图像版，不是「画风归因」。

---

## 6. 多模态连续学习 — 车道 A

业务：视频/图/文塔在连续学。某一窗任务或匹配律 hop 了，该刷新哪条学习路径，但 **先** 问映射还在不在。

**一张表**

- Y：这一条多模态样本的任务标签。检索命中、配对对错、caption 是否匹配。一个数。
- X：各塔向量拼在一起（或投影后）。预测 Y，不在这条里拆塔贡献。
- batch：时间或任务阶段（pretrain → 检索 → 生成）。

**可安排的数据**

| 数据 | Y | 切连续学习 |
|---|---|---|
| MSR-VTT | 检索/匹配；caption 对视频 | 按视频 id 序或 split 当时间；仓库其它 branch 已做过 |
| Fashion-IQ | 相对描述下的目标图 | 属性/品类切阶段 |
| COCO captions | 图文匹配 | 年份/split 当窗 |
| WebVid | 视频-文本 | 真大规模流，子集即可 |

连续学习的 hop：cut 后改标签律（例如从「检索命中」改成「风格匹配」），或改标注规范。quiet 表示当前头还能当教师。

---

## 7. 多模态塔对照 — 并排，不是归因

和 6 同一批数据，X 切开：`X_video | X_text | X_audio`。  
每块只用自己的 X，问上一窗和这一窗还像不像。两块并排：谁动了、谁没动。

不把「视频动了」写成「视频导致了 Y」。没有分塔表征就空着这一列。细节：`docs/reports/Tower_Style_Mix_Contrast.md`。

---

## 8. 特征族 / mix 对照 — 并排，不是归因

审核表按长度 / 语气 / 拒绝切开；推荐按序列 / 交叉 / 上下文切开。  
mix 是来源和任务占比，不是风格向量。

三列可以各动各的：mix 换了、文风还在、映射仍 quiet。只有映射 fire 才动审核/回滚。细节同上对照笔记。

---

## 9. 模型迭代更新 — 车道 A（门禁）

业务：新 checkpoint / 新 DPO 合并 / 新 judge **能不能替换旧的**。

**一张表**

- 同一批 X（同一组 prompt / 同一组请求）。
- Y_old = 旧模型（或旧 judge）的决定；Y_new = 新的。
- 把 Y_old 当 T=0 的标签学 probe，拿 Y_new 当 T=1 的标签打 OOS：fire ⇒ 新制度，拒全量替换，灰度 + Top-k。

**可安排的数据**

| 数据 | 怎么当「两版本」 |
|---|---|
| HH 金标批 vs 翻转/合成偏好批 | **已有** hop 构造 |
| 同一 HaluEval 题，忠实答 vs 幻觉答 | 两个「模型版本」的输出 |
| AlpacaEval / Arena 同 prompt 两模型 | Y=哪边赢，或各模型自己的裁判分 |
| 自有 A/B：v1/v2 serving 日志 | 最真；X 对齐 request_id |

迭代之后若要看「视频塔动了还是文案动了」，用第 7/8 条对照，不要写成版本更新的根因。

---

## 已经在本仓库排上的，和还该拉的

**已经在盘上、能直接喂 A**

| 包 | 路径 | 覆盖场景 |
|---|---|---|
| HH 审核流 | `data/hf_cache/audit/*_600.jsonl` + `xy_*.csv` | 4，兼 9 的门禁 |
| HaluEval / SQuAD / HotpotQA / TruthfulQA | `data/hf_cache/infer_bench_export/` | 2，兼 3 的弱质检 |
| HaluEval landing 大流 | 其它 branch `halueval_qa_3000.jsonl` | 2 的主数字 |

**值得下一轮按场景拉（仍是一个 Y + 一堆 X）**

1. **BeaverTails 或 WildGuard** → 审核 Y 换成真人/专业安全标签（场景 4）。
2. **RAGTruth** → 推理监控带上检索轴，避免把 RAG 缺口当成生成 hop（场景 2）。
3. **WikiArt 或 DiffusionDB 子集** → 画风只走车道 B（场景 5）。
4. **MSR-VTT 子集** → 场景 6 的 (Y, X) 流；同一份切开视频/文本做塔对照（场景 7）。
5. **Amazon Reviews 带时间戳** → 场景 1 的经典 ML 上线流。

---

## 接线（不要串台）

```text
样本流 (Y, X, batch)
        │
        ├─ A OnlineRFPerm ── fire? ── 监控 / 审核 / 推理 / 质检 / 迭代门禁
        │                         quiet → Y 仍可当金标，w=1
        │                         fire  → 对照窗；Top-k po_risk0；cap 合并
        │
        ├─ B P(X) 文风/画风/mix ── 动/不动并排（另账）
        │
        └─ 塔 / 特征族切开 ── 同一时钟谁的画像不像了（对照，不是归因）
```

场景 1–4、6、9：先把 **一个 Y、一堆 X** 的 csv/parquet 排出来，和现在的 `xy_hh_*.csv` 同一形态。  
场景 5、7、8：文风/画风/mix 和塔/族，都是画像对照，详见 `docs/reports/Tower_Style_Mix_Contrast.md`。
