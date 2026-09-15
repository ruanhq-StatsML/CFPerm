# HH 数据整理 · 画风/文风检测 · 后续 procedures 详解

> 聚焦 **data manipulation**（不是连续 batch 检测本身）。  
> 代码主入口：`scripts/agod/hh_online_rfperm_stream.py`  
> 画风对照面更完整版：`scripts/agod/hf_landing_protos.py`（`style_vector` + `domain_auc`）

---

## 0. 总流程（你要的那条）

```text
HH raw dialog (chosen/rejected)
        │
        ├─ human_prompt(chosen)  → question
        ├─ assistant_reply(*)    → answer
        ├─ label chosen=1/rej=0  → y_raw
        └─ style_feature(answer) → 文风/画风 10 维
                │
                ▼
        tidy: t_idx | batch | question | answer | y | feature
                │
                ├─ (可选) style_domain_auc：低/高正式度 tercile 可分性
                └─ (可选) OnlineRFPerm：用 feature+y 做制度火 / po_risk0
```

---

## 1. `human_prompt` / `assistant_reply` — 为什么要拆

HH-RLHF 一行是整段对话：

```text
\n\nHuman: <问题>\n\nAssistant: <回答>
```

`chosen` / `rejected` 各自带同一 Human、不同 Assistant。整理时必须拆成 **question / answer**，否则：

- feature 会混进 Human 文本，文风轴被问题长度污染；
- 同一 prompt 两行 answer 对不上，后续偏好标签无锚点。

### 1.1 `human_prompt`

```94:96:scripts/agod/hh_online_rfperm_stream.py
def human_prompt(dialog: str) -> str:
    m = re.search(r"Human:\s*(.*?)(?:\n\nAssistant:|$)", dialog, flags=re.S)
    return (m.group(1).strip() if m else dialog.strip())[:800]
```

| 点 | 说明 |
|----|------|
| 匹配 | `Human:` 到下一个 `\n\nAssistant:`（或文末） |
| 截断 | 800 字符：防超长 prompt 撑爆 hash/表 |
| 取自 | **只用 `chosen` 的 Human**（与 rejected 同题，避免双写不一致） |
| 失败 | 无匹配则整段 strip（兜底） |

### 1.2 `assistant_reply`

```99:101:scripts/agod/hh_online_rfperm_stream.py
def assistant_reply(dialog: str) -> str:
    parts = re.split(r"\n\nAssistant:", dialog)
    return parts[-1].strip() if len(parts) >= 2 else dialog.strip()
```

| 点 | 说明 |
|----|------|
| 取法 | 按 `\n\nAssistant:` 切，**最后一段** = 最终助手回复 |
| 为何 last | 多轮时只要最终 answer；HH helpful-base 多为单轮 |
| 用途 | ① 写入 `answer` 列 ② **画风特征只打在 answer 上** |

### 1.3 一对 → 两行（manipulation 核心）

```164:171:scripts/agod/hh_online_rfperm_stream.py
    for r in rows:
        q = human_prompt(r["chosen"])
        for key, lab in (("chosen", 1), ("rejected", 0)):
            a = assistant_reply(r[key])
            questions.append(q)
            answers.append(a)
            labels.append(lab)
            styles.append(style_feature(a))
```

结果：同一 `question`，两行 `answer`，`y∈{0,1}`。这就是 tidy 表的原子行。

---

## 2. 画风 / 文风检测逻辑（你说不错的那块）

### 2.1 对象：\(P(X)\) 画像，不是偏好对错

- **文风/画风轴**：回答「说话语气 / register 变了没」→ `style_feature` / `style_domain_auc`
- **偏好轴**：回答「chosen/rejected 映射变了没」→ OnlineRFPerm on `y`
- **禁止**：用文风 AUC 解释对齐失败；文风工单 ≠ 客服工单

### 2.2 `style_feature` 十维（打在 **answer** 上）

| 维 | 含义 | Intuition |
|----|------|-----------|
| 0 | 词数 / 200 | 长短：啰嗦 vs 短答 |
| 1 | 字符 / 800 | 体积 |
| 2 | 均词长 / 10 | 用词复杂度 |
| 3 | `?` 密度 | 反问/不确定 |
| 4 | `!` 密度 | 情绪/促销感 |
| 5 | hedge 词（maybe/perhaps/…） | 含糊语气 |
| 6 | formal 词（therefore/however/…） | **正式度**（画风主轴之一） |
| 7 | ` i ` 密度 | 第一人称口语 |
| 8 | 换行密度 | 结构/列表感 |
| 9 | 大写占比 | 喊麦/标题党 |

**为什么好用：** 便宜、可审计、streaming；对齐漂移常先表现为「口气变了」。  
**不够什么：** 不解语义；同义正式改写可能漏检 → 生产叠 embedding。

### 2.3 还有：`style_domain_auc`（hf_landing 对照面）

在 `hf_landing_protos.py` 的 Judge demo 里，不只存 10 维，还做：

1. `formality ≈ style[:,6] + 0.5*style[:,1]`（正式度 + 长度）
2. 取低/高 tercile 两堆 style 向量
3. `domain_auc(lo, hi)`：RF 分类「低正式 vs 高正式」的 AUC  
   - **高** ⇒ 文风画像可分（\(P(X)\) 有结构）  
   - 与 `judge_err_ratio` / OnlineRFPerm fire（\(P(Y\mid X)\)）**分列**

这是「画风检测」的 **可分性指标**；`style_feature` 是 **逐条画像**。

### 2.4 和 hash 怎么配

| featurizer | 内容 | 何时 |
|------------|------|------|
| `style_only` | 仅 10 维画风 | 只盯语气轴 |
| `hash_only` | 仅 HashingVectorizer | 不推荐生产 |
| `style_hash`（默认） | style ⊕ hash | demo 闭环；生产换 embedding 替换 hash 段 |

---

## 3. 后面的 procedures（拆步）

### Procedure A — tidy 表（data manipulation 本体）

1. HH jsonl → 逐对展开 chosen/rejected  
2. `human_prompt` / `assistant_reply`  
3. `style_feature(answer)`  
4. **permute**：打散 chosen/rejected 交错，避免假时间结构  
5. 切整 batch：`t_idx`，`batch = t_idx // n_per`  
6. 拼 `feature`（style / hash / 二者）  
7. （仅 demo）cut 后翻转部分 `y` 模拟偏好 hop；线上用真标签  
8. 落盘：`stream_table.parquet` + `feature_matrix.npy`

产出列：`t_idx | batch | question | answer | y | feature`（另有 `y_raw_pref` / `flipped` 审计用）。

### Procedure B — 画风评估（可选，不对齐主账）

1. 对 `feature` 里 style 段做 formality  
2. tercile 切分 → `style_domain_auc`  
3. 高 AUC → 开 **素材/decoding** 工单，**不进**客服¥账

### Procedure C — 制度检测（你已清楚，一笔带过）

`feature + y + batch` → OnlineRFPerm：`fit_online_probe` → fire？→ `po_risk0` / `√PO`。  
**输入已是 Procedure A 的表**；换 embedding 只换 `feature` 列。

### Procedure D — 对齐动作

| 信号 | 动作 |
|------|------|
| style_auc 高、concept 未火 | 只动文风/素材 |
| OnlineRFPerm fire | 拒合并 / 限量 / Top-k 复核 |
| fire + 高 po | √PO 重加权或人工审偏好对 |

### Procedure E — 接到业务¥（另一条线）

HaluEval hop knobs → CS seed → 少工单/退款/**多承接**。  
文风轴默认 **不**进这条。

---

## 4. 还有其他的吗？

| 模块 | 文件 | 作用 |
|------|------|------|
| HH Judge + style_auc | `hf_landing_protos.py` | 对齐门禁 + 画风可分性 |
| HH tidy + style10 | `hh_online_rfperm_stream.py` | **你要的 form** |
| HaluEval 幻觉/RAG | 同 `hf_landing_protos.py` | 推理制度；rag_hit 路由 |
| 客服¥账 | `sql/biz_value/04–14` | act-vs-ignore × 单价 |
| 包 API | `agod.online_rfperm` / `po_iptw` | 检测与 √PO |

**和 data-manipulation 最相关、可继续加深的：**

1. 多轮对话：`human_prompt` 取首轮 vs 末轮 Human 的策略  
2. 画风维扩展：emoji / 中英混排 / 模板句比例  
3. `feature = [style ‖ embedding]` 固定接口，检测不动  
4. style_auc 与 preference fire 的 **联合门禁票**（已在 Roadmap 写过）

---

## 5. 怎么自己跑 manipulation

```bash
PYTHONPATH=. python3 scripts/agod/hh_online_rfperm_stream.py --featurizer style_only
# → results/agod/hh_online_stream/stream_table.parquet
# → stream_table_preview.csv（人眼看 question/answer）
```

只关心表：读 parquet 的 `t_idx, batch, question, answer, y`；`feature` 可换成你们向量后再接 Procedure C。
