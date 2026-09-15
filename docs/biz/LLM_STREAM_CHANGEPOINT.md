# LLM 推理变点：两数据集 manipulation → OnlineRFPerm

## 数据 form（就这一列时间流）

```python
{
  "t": int,
  "question": str,
  "answer": str,
  "embedding": np.ndarray,   # shape (p,) or (p, 1)
  "score": float,            # Halu=幻觉0/1；HH=偏好0/1
}
```

## 方法怎么用

| 问题 | 方法 | 本原型 |
|------|------|--------|
| **大模型推理有没有制度/质量变点？** | **OnlineRFPerm** | ✅ 两数据集已跑 |
| 实例风险排序 / PO-risk | **BOCPD**（你们实验已论证更好） | 不在此硬塞 OnlineRFPerm |

Stationary-DGP robustness 已在算法侧论证过；**这批数据只看 detection delay**（true cut → first fire），没什么玄学。

```python
# stream: list[dict] with embedding + score
from agod.online_rfperm import fit_online_probe, probe_err, hop_fires, error_floor
# 每个新 batch：用上一窗拟合探针 → 本窗算 e_now → 与 e_prev 比 → fire?
# fired == 变点/制度火
# delay_batch = first_fire_batch - cut_batch
```

> PO-risk 实例排序请用 **BOCPD**（实验已论证）；本脚本只做 OnlineRFPerm 变点。

## Detection delay（本跑）

| 流 | cut | first fire | delay (batch) | delay (t) | ratio@fire |
|----|----:|----------:|--------------:|----------:|-----------:|
| HaluEval（推理幻觉变点） | 4 | 4 | 0 | 0 | 30.66666666666664 |
| HH-RLHF（偏好映射变点） | 4 | 4 | 0 | 0 | 31.999999999999975 |

## 本跑

### HaluEval（推理幻觉变点）

| 项 | 值 |
|----|-----|
| 源 | `data/hf_cache/halueval_qa_3000.jsonl` |
| n / p / n_per | 1200 / 138 / 100 |
| cut_batch | 4 |
| fires | 1 / 11 |
| first_fire_batch | 4 |
| **detection_delay** | **0 batch**（0 条） |
| ratio@first_fire | 30.66666666666664 |

样例时刻 t=0::

```json
{
  "t": 0,
  "question": "Which magazine was started first Arthur's Magazine or First for Women?",
  "answer": "First for Women was started first.",
  "embedding": "np.ndarray shape [138, 1]",
  "score": 1.0
}
```

| batch | t | fired | ratio | score均值 |
|------:|---|:-----:|------:|----------:|
| 1 | 100–199 | False | 1.000 | 0.03 |
| 2 | 200–299 | False | 1.000 | 0.03 |
| 3 | 300–399 | False | 1.000 | 0.03 |
| 4 | 400–499 | True | 30.667 | 0.92 |
| 5 | 500–599 | False | 0.087 | 0.92 |
| 6 | 600–699 | False | 1.000 | 0.92 |
| 7 | 700–799 | False | 0.875 | 0.93 |
| 8 | 800–899 | False | 1.143 | 0.92 |
| 9 | 900–999 | False | 0.875 | 0.93 |
| 10 | 1000–1099 | False | 1.143 | 0.92 |
| 11 | 1100–1199 | False | 1.000 | 0.92 |

### HH-RLHF（偏好映射变点）

| 项 | 值 |
|----|-----|
| 源 | `data/hf_cache/hh_rlhf_helpful_base_2500.jsonl` |
| n / p / n_per | 2000 / 139 / 80 |
| cut_batch | 4 |
| fires | 8 / 24 |
| first_fire_batch | 4 |
| **detection_delay** | **0 batch**（0 条） |
| ratio@first_fire | 31.999999999999975 |

样例时刻 t=0::

```json
{
  "t": 0,
  "question": "Do you need a special license to drive an RV?",
  "answer": "I’m glad I could help, have a nice day.",
  "embedding": "np.ndarray shape [139, 1]",
  "score": 0.0
}
```

| batch | t | fired | ratio | score均值 |
|------:|---|:-----:|------:|----------:|
| 1 | 80–159 | False | 1.000 | 0.38 |
| 2 | 160–239 | False | 0.556 | 0.51 |
| 3 | 240–319 | False | 0.200 | 0.54 |
| 4 | 320–399 | True | 32.000 | 0.54 |
| 5 | 400–479 | False | 0.469 | 0.41 |
| 6 | 480–559 | False | 0.400 | 0.54 |
| 7 | 560–639 | True | 2.917 | 0.51 |
| 8 | 640–719 | False | 0.429 | 0.50 |
| 9 | 720–799 | True | 1.267 | 0.46 |
| 10 | 800–879 | True | 1.474 | 0.47 |
| 11 | 880–959 | False | 0.643 | 0.47 |
| 12 | 960–1039 | False | 1.167 | 0.50 |


## 产物

| 文件 | 内容 |
|------|------|
| `results/agod/llm_changepoint/halu_stream.parquet` | Halu 时间流 |
| `results/agod/llm_changepoint/halu_stream.npy` | embedding 矩阵 |
| `results/agod/llm_changepoint/hh_stream.parquet` | HH 时间流 |
| `results/agod/llm_changepoint/*_hops.json` | OnlineRFPerm fire 表 |
| `results/agod/llm_changepoint/summary.json` | delay / fires |

```bash
PYTHONPATH=. python3 scripts/agod/llm_stream_changepoint.py
```
