# 已做好的特征在哪 + 样例 Prototype

## 路径

| 文件 | 内容 |
|------|------|
| `results/agod/hh_online_stream/stream_table.parquet` | tidy 表：`t_idx, batch, question, answer, y, feature` |
| `results/agod/hh_online_stream/feature_matrix.npy` | 矩阵 **(2000, 138)** = style10 ⊕ hash128 |
| `results/agod/hh_online_stream/stream_table_preview.csv` | 无向量的人眼预览 |
| `results/agod/hh_online_stream/meta.json` | featurizer / cut / n_per |
| `data/hf_cache/hh_rlhf_helpful_base_2500.jsonl` | 原始 HH-RLHF 缓存 |

meta：`{"cut_batch": 4, "n_per": 80, "featurizer": "style10+hash128", "feature_dim": 138}`

## 特征怎么拼的

```text
HH dialog
  → human_prompt(chosen)  → question
  → assistant_reply(chosen/rejected) → answer
  → style_feature(answer) → 10 维画风/文风
  → HashingVectorizer(question+"\n"+answer) → 128 维弱语义
  → feature = concat(style10, hash128)     # 共 138
```

**画风检测 intuition：** `formal` / `hedge` / 长度 / 标点 描述的是 \(P(X)\) 说话画像，不是偏好对错。  
高 `style_domain_auc`（见 `results/agod/hf_landing/hh_judge_style.json`）说明文风可分；偏好 hop 另走 OnlineRFPerm + `y`。

## 样例

### 样例 t_idx=0（batch=0, y=1）

- **question:** How can you learn to be polite
- **answer:** I have an idea about a system I’d like to try.  I’d first model politeness, then teach you to copy it, and then reward you for good behavior.…
- featurizer=`style10+hash128` · dim=138 · flipped=0

| style 维 | 含义 | 值 |
|----------|------|-----|
| `n_tok` | 词数/200 | 0.1500 |
| `n_char` | 字符数/800 | 0.1762 |
| `avg_word` | 平均词长/10 | 0.3567 |
| `qmark` | ? 密度 | 0.0000 |
| `bang` | ! 密度 | 0.0000 |
| `hedge` | maybe/perhaps/... 犹豫词 | 0.0000 |
| `formal` | therefore/however/... 正式词 ← 画风主轴之一 | 0.0000 |
| `i_count` | 第一人称 i | 0.0000 |
| `newlines` | 换行密度 | 0.0000 |
| `upper_ratio` | 大写比例 | 0.0213 |

hash128 ‖·‖₂ ≈ 1.000；head8 = `[0.09759000729485333, 0.09759000729485333, 0.09759000729485333, 0.09759000729485333, 0.0, 0.0, 0.0, 0.09759000729485333]`


## 怎么读（三行 Python）

```python
import pandas as pd, numpy as np
df = pd.read_parquet("results/agod/hh_online_stream/stream_table.parquet")
X  = np.load("results/agod/hh_online_stream/feature_matrix.npy")  # (2000, 138)
# X[i, :10] = style；X[i, 10:] = hash
```

生产：只用 parquet 的 `question/answer/y`，把 `feature` 换成你们 embedding 即可。

## 还有其他相关产物吗？

| 产物 | 路径 | 用途 |
|------|------|------|
| HH Judge + style_auc | `results/agod/hf_landing/hh_judge_style.json` | 文风 AUC / judge_err_ratio |
| Halu 推理 hop | `results/agod/hf_landing/halu_regime_rag.json` | 幻觉制度 knobs |
| 连续检测 hop 表 | `results/agod/hh_online_stream/hop_history.json` | fire / po Top-k |
| 生成脚本 | `scripts/agod/hh_online_rfperm_stream.py` | data-manipulation |

```bash
PYTHONPATH=. python3 scripts/agod/show_prepared_features.py
```
