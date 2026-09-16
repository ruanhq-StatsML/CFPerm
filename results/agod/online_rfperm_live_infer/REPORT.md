# OnlineRFPerm × Live LLM Infer（HaluEval）

> Use-case：真生成流上挂 OnlineRFPerm → fire → 动作路由。  
> **不是**新 LLM 算法；**是**推理框架输出 → 制度火 → 工单。

## Setup

| 项 | 值 |
|----|-----|
| backend | `transformers` |
| model | `/workspace/data/models/SmolLM2-135M-Instruct` |
| data | HaluEval QA (`halueval_qa_3000.jsonl`) |
| n / n_per / batches / cut | 96 / 16 / 6 / 3 |
| X dim | 73 (hash64 ⊕ style8 ⊕ rag1) |
| gate / audit_k | 1.25 / 5 |

## Regime

- **quiet** (batch < cut)：system=faithful + knowledge in prompt  
- **hop** (batch ≥ cut)：system=invent + **no knowledge**（生成制度 + 检索缺口同时跳）

## Quality shift

| 窗 | mean y_bad | mean rag_hit |
|----|------------|--------------|
| quiet | 0.167 | 0.307 |
| hop | 0.333 | 0.218 |

## OnlineRFPerm

| 项 | 值 |
|----|-----|
| first_fire_batch | 3 |
| **detection_delay** | **0 batch** |
| cut_batch | 3 |

| batch | regime | fired | ratio | y_bad | rag_hit | action |
|------:|:------:|:-----:|------:|------:|--------:|--------|
| 0 | quiet | False |  | 0.12 | 0.314 | `quiet_pass` |
| 1 | quiet | False | 4.000 | 0.25 | 0.268 | `quiet_pass` |
| 2 | quiet | False | 0.500 | 0.12 | 0.340 | `quiet_pass` |
| 3 | hop | True | 2.500 | 0.31 | 0.227 | `model_rollback_or_audit_topk` |
| 4 | hop | False | 0.800 | 0.25 | 0.227 | `quiet_pass` |
| 5 | hop | True | 1.250 | 0.44 | 0.200 | `model_rollback_or_audit_topk` |

## 对外一句

真推理流（`transformers`）上 OnlineRFPerm：cut=3 → first_fire=3，
**delay=0**；quiet→hop 时 y_bad 0.17→0.33，
动作按 rag 路由到 `model_rollback_or_audit_topk`。

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_live_infer.py
# vLLM:
PYTHONPATH=. python3 scripts/agod/online_rfperm_live_infer.py \
  --backend vllm --base-url http://127.0.0.1:8000/v1 --model <served>
```
