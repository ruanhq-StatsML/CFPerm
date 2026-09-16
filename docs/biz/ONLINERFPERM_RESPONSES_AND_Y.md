# Responses · answer-precision · $Y\in\{0,1\}$ · streaming note

> Prototype 钉死四件事：  
> 1) 各 dataset 的 **response 长什么样**  
> 2) faith 用 **answer-precision**（端答案 / 常温稀释），**不用 Jaccard**  
> 3) **$Y$ 是二值 $\{0,1\}$**；faith 本身才是 $[0,1]$ 连续分  
> 4) 默认 **$n_{\mathrm{per}}=20$ batch 钟**；`n_per=1` streaming 会变，见 §5

代码：`scripts/agod/online_rfperm_live_infer.py`（`answer_precision`）  
表：`results/agod/infer_dataframes/*_stream_table.parquet`

---

## 1. 各 dataset 的 response 是什么

同一 protocol：quiet = knowledge 进 prompt + faithful system；hop = 去 knowledge + invent system。  
Mock quiet 回声知识前 ~18 token（端答案形）；hop 固定 invent 句（Atlantis / Dr. Fabricatus）。  
真模型（transformers/vLLM）同 system，只是生成器不同。

| dataset | knowledge 包装 | quiet response（典型） | hop response（典型） | gold |
|---------|----------------|------------------------|----------------------|------|
| **HaluEval** | 短 knowledge snippet | 复述 / 紧贴 snippet 的一两句 | 自信胡编，与 knowledge 脱节 | 短答案句 |
| **SQuAD** | 段落 context | 从段落抽出的短答或段落前缀 | 同上 invent | 短 span（如队名） |
| **HotpotQA** | 多段 distractor（很长） | quiet 仍短；端答案嵌在长 K 里 | invent | 很短（yes / 地名） |
| **TruthfulQA** | correct-answer bank（`\|` 拼接） | 贴近正确答法银行 | invent | best_answer |

样例（mock）：

**HaluEval quiet**  
- Q: Which magazine was started first Arthur's Magazine or First for Women?  
- A: *Arthur's Magazine (1844–1846) was an American literary periodical…*（知识前缀）  
- gold: First for Women was started first.  
- faith (answer-precision) ≈ 1.0 → **$Y=0$**

**HaluEval hop**  
- A: *I am certain it was founded in Atlantis-42 in 1899 by Dr. Fabricatus regarding: …*  
- faith ≈ 0.2 → **$Y=1$**

**SQuAD quiet**  
- Q: Which NFL team represented the AFC at Super Bowl 50?  
- A: 段落前缀；gold = `Denver Broncos`；faith≈1 → $Y=0$

**HotpotQA quiet（关键：常温稀释）**  
- A 仍短，K 却是 ~10 段。  
- Jaccard ≈ 0.07 → 旧路径 quiet $Y\approx1$（假坏）。  
- answer-precision：端答案 token 几乎都在 K 里 → faith≈1 → **$Y=0$**，hop 后才 $Y=1$。

**TruthfulQA quiet**  
- A 贴近 `correct_answers` 拼接串；hop 仍 invent。

---

## 2. 为什么不用 Jaccard：端答案 + 常温稀释

| 度量 | 公式 | 长 K 时短 A |
|------|------|-------------|
| Jaccard（旧） | $\|A\cap K\|/\|A\cup K\|$ | $\|A\cup K\|\approx\|K\|$ → 分数被稀释到 ≈0 |
| **answer-precision（现）** | $\|A\cap K\|/\|A\|$ | 短答若嵌在 K 里 → ≈1 |

「端答案」= 短、贴近 gold / 抽取 span 的 serving 输出。  
「常温稀释」= 把短 A 丢进长 ambient K 用并集归一，把 grounded 答也打成「不忠实」。

Hotpot 对照（同生成）：legacy Jaccard-full quiet→hop ≈ 1→1；prec quiet→hop ≈ 0→1（见 `ONLINERFPERM_HOTPOT_PROTO.md`）。

默认阈值：`faith_thr=0.45`（precision 尺度；旧 Jaccard 的 0.18 不再适用）。

```text
faith = |tokens(answer) ∩ tokens(knowledge)| / |tokens(answer)|   ∈ [0,1]
Y     = 1{ faith < thr }  (hop 窗可叠加 invent 启发式)
```

---

## 3. $Y$ 是不是 0–1？——拆两层

| 量 | 取值 | 用途 |
|----|------|------|
| **`faith` / `rag_hit`** | **连续** $[0,1]$ | answer-precision 分数；审计、画图、调 thr |
| **`y` / `y_bad`** | **二值** $\{0,1\}$ | OnlineRFPerm / OOB 探针的分类标签 |

- $Y=0$：本 turn「好」（faith 够）  
- $Y=1$：本 turn「坏」（faith 低，或 hop+invent）  
- **不是**把 continuous faith 直接当回归 $Y$。探针 task = accuracy / 分类。  
- parquet：`y` int ∈ {0,1}；`faith` float ∈ [0,1]；`faith_metric=answer_precision`

生产换 $Y$ 源仍保持二值（或先阈值化）：  
`human_correction` > judge-thresholded > escalation > synthetic precision。

---

## 4. Prototype 验收（默认 batch，$n_{\mathrm{per}}=20$）

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_multi_datasets.py \
  --backend mock --n-per 20 --cut-batch 5 --n-batches 10 --faith-thr 0.45
```

| dataset | fire | delay | $y$ quiet→hop | faith quiet→hop（mock） |
|---------|-----:|------:|---|---|
| HaluEval | 5 | 0 | 0→1 | ≈1.00→0.24 |
| SQuAD | 5 | 0 | 0→1 | ≈1.00→0.20 |
| HotpotQA | 5 | 0 | 0→1 | ≈1.00→0.38 |
| TruthfulQA | 5 | 0 | 0→1 | ≈1.00→0.21 |

表形不变：`t_idx|batch|question|answer|y` + `faith`/`faith_metric`。

---

## 5. `n_per=1` = streaming detection —— 会变，默认不切

对：`n_per=1` 就是把每个 turn 当成一个 batch，detection 变成 observation 级 streaming。

**会变什么**

| 项 | $n_{\mathrm{per}}=20$（默认） | $n_{\mathrm{per}}=1$（streaming） |
|----|-------------------------------|-------------------------------------|
| 钟 | batch 窗 | 单点 turn |
| `n_ref≥100` | cut=5 | cut≥100 |
| ORF 探针 | 每窗 fit 上一 batch | 若仍「上一窗=1 点」则 RF 极噪 |
| `first_k` | trail **batch** 下标 | trail **obs** 下标 |
| 经验（mock Hotpot） | delay **0** | 同套 ratio-gate **可不火**（单点 fit 不稳） |

对照：

```bash
# streaming-shaped clock (same tidy form)
PYTHONPATH=. python3 scripts/agod/online_rfperm_multi_datasets.py \
  --backend mock --n-per 1 --cut-batch 100 --n-batches 200 \
  --datasets hotpotqa --out results/agod/online_rfperm_stream_nper1
# → y 仍 0→1，但 ORF fire 可能 ---（单点 hop_fires 不稳）
```

若要真正的 streaming 与 manuscript `onlinePermOOB` 对齐：  
**冻结** reference 上的探针 → 对 trail 打 OOB/error 标量流 → BOCPD/PH/ADWIN/`first_k`（`batch_size=1`），而不是每步用 1 个点重 fit ORF。

**就这样**：batch 主表仍钉 $n_{\mathrm{per}}=20$；连续 streaming 另开正式对照（见下）。

---

## 6. Streaming testing（连续弄法）— 已 prototype

脚本：`scripts/agod/online_rfperm_streaming_test.py`  
报告：`docs/biz/ONLINERFPERM_STREAMING_TEST.md`

| 模式 | 连续怎么走 |
|------|------------|
| **orf_freeze** | quiet 上 fit 浅层 RF；trail 逐点打 0/1 error；滚动均值 / quiet ≥ γ → fire |
| **orf_slide** | 窗长 `win`、步长 1；OnlineRFPerm 连续 OOS ratio-gate |
| **smooth_control** | 把 hop 行换成 quiet 重采样；应 **不火** |

当前四套（win=20, γ=1.25）：freeze delay=0，slide delay=1；smooth 全静。

### 怎么理解（你说的那句）

> 用 Random Forest 做 component model，有质量制度 hop 就该检出；  
> **哪个不火，说明在这个探针下推理链路很丝滑**（没有可用的 $P(Y\mid X)$ 跳变）。

对，可以这样读——附加两句边界：

1. 前提是 $Y$/特征诚实（answer-precision 等）。标签烂也会「假丝滑」。  
2. smooth_control 不火 = gate 没胡乱报警；hop 火 + control 静 = RF 真看到了制度差。

---

## 7. Fracture（推理断裂）— 简单但实用

> 不火 ≈ 推理正常；中间突然 fracture / concept drift → 给断裂打时间戳。

脚本：`scripts/agod/online_rfperm_fracture_perturb.py`  
报告：`docs/biz/ONLINERFPERM_FRACTURE_PERTURB.md`

扰动（$t=n_{\mathrm{ref}}$ 起突然变）：

| kind | 断裂类型 |
|------|----------|
| `invent_fracture` | 生成制度断：invent + 去 grounding |
| `label_flip` | concept drift：$Y$ 翻转 |
| `answer_corrupt` | 忠实度断：答案打乱 |

六套卡（HaluEval / SQuAD / Hotpot / TruthfulQA / BoolQ / NQ-open）× 三种扰动：freeze delay=0，slide delay=1；smooth 全静。
