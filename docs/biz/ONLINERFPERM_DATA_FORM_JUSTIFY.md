# Justify: LLM infer → tidy serving form → OnlineRFPerm

> **一句**：这条 use-case 的硬 justify 是 **data formulation**——把 HaluEval / SQuAD / HotpotQA / TruthfulQA（以及生产助手日志）压成与 HH 偏好流同构的 tidy 表；探测器、gate、路由都是复用。

主文：`docs/biz/ONLINERFPERM_INFER_DRAFT.tex` §Formulation of data

---

## 原子行

```
t_idx | batch | question | answer | y
```

| 列 | 含义 | OnlineRFPerm 角色 |
|----|------|-------------------|
| `t_idx` | serving 序 | 在线钟 |
| `batch` | `t_idx // n_per` | 窗 $B_t$（fit / OOS） |
| `question` | 查询 | 进 $X$ |
| `answer` | 模型输出 | 进 $X$；生产换 serving log |
| `y` | 坏例 $\{0,1\}$ | $Y$ |

审计/路由列（`hopped`, `rag_hit`, `faith`, `knowledge`, `gold`, `system`）**不进 fire 规则**，只解释 $y$ 怎么来、火了动哪条臂。

两层同一 schema：

| 层 | 文件 | 内容 |
|----|------|------|
| source | `*_source.parquet` | question / knowledge / gold |
| stream | `*_stream_table.parquet` | answer + $y$ + batch clock |

路径：`results/agod/infer_dataframes/`

---

## 与 HH 的同构

HH 已是 `t_idx | batch | question | answer | y | feature`。  
推理 QA 同一骨架；变的是 $y$ 语义（偏好 → 忠实度 / 纠错 / 转人工），不是表形。  
换 embedding 只换 $X$，不换表。

**对外一句**：LLM serving 进 OnlineRFPerm，靠的是变成 HH 形 tidy 流。

---

## 为什么四套语料能压在一起

知识包装不同（snippet / passage / multi-para / answer bank），但每一 turn 都有 query、生成答案、可打分质量位。  
压成一表 → 一套 clock（$n_{\mathrm{per}}$, cut, $n_{\mathrm{ref}}$）+ 一套 `first_k`，而不是四个 ad-hoc 脚本。

Quiet/hop 是**列纪律**，不是第二种文件格式：

- quiet（batch < cut）：带 knowledge，`hopped=False`
- hop（batch ≥ cut）：去 knowledge + invent，`hopped=True`
- $n_{\mathrm{ref}} = n_{\mathrm{per}} \times \mathrm{cut}$（20×5=100）

---

## 生产换列（不换形）

保留 parquet 列；`answer` ← serving log；`y` 优先级：

`human_correction` > judge > escalation > synthetic faithfulness

---

## 这张表故意不装什么

单句事实对错、因果根因、¥。那些挂在 fire **之后**（对照窗 / 路由 / 可选账本）。  
formulation 停在「让 monitor 能接上」。

其余（探测器对比、路由、path 子集）都好说——先钉住这张表。
