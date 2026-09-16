# LLM 推理 / 对齐 Prototype（call ``agod``）

## 你 call 哪个包

```python
from agod.online_rfperm import (
    fit_online_probe, hop_fires, po_risk0_rows,
    probe_err, run_rfperm_stream,
)
from agod.po_iptw import po_iptw_weights
from agod.po_refit import Stream
```

## 数据在哪

| 用途 | 路径 | 源 |
|------|------|----|
| 对齐 / Judge | `data/hf_cache/hh_rlhf_helpful_base_2500.jsonl` | Anthropic/hh-rlhf |
| 推理 / 幻觉 | `data/hf_cache/halueval_qa_3000.jsonl` | pminervini/HaluEval |

## PO 核心逻辑（贴这段）

```python
probe = fit_online_probe(X_prev, y_prev, task='acc')
e_now = probe_err(probe, X_cur, y_cur, task='acc')
fire  = hop_fires(e_now, e_prev, gate=1.25)
po    = po_risk0_rows(probe, X_cur, y_cur, task='acc')
w     = po_iptw_weights(po, mode='sqrt') if fire else np.ones_like(po)
audit = np.argsort(-po)[:10]   # Top-k 人工 / 升配
```

## 本跑结果

### 推理（幻觉制度）

| 项 | 值 |
|----|-----|
| 数据 | `data/hf_cache/halueval_qa_3000.jsonl` |
| package | `agod.online_rfperm / agod.po_iptw` |
| n | 1200 |
| stream fires | 2 |
| hop fired / ratio | `True` / `7.80240229115313` |
| po_risk0 P@10 / AUROC | `1.0` / `1.0` |
| 单步 gate fired / mean_w | `True` / `0.9999989075843194` |
| 动作 / 门禁 | `model_rollback_or_audit_topk` |

- serving：`{'on_fire': '升配检索 / 人工审计 Top-k(po_risk0) / 可选 √PO 重加权下一跳训练', 'on_quiet': '标准路径，w=1'}`

### 对齐（RFPerm-as-Judge）

| 项 | 值 |
|----|-----|
| 数据 | `data/hf_cache/hh_rlhf_helpful_base_2500.jsonl` |
| package | `agod.online_rfperm / agod.po_iptw` |
| n | 1600 |
| stream fires | 0 |
| hop fired / ratio | `False` / `0.9807708676608308` |
| po_risk0 P@10 / AUROC | `0.4` / `0.4896291640477688` |
| 单步 gate fired / mean_w | `False` / `1.0` |
| 动作 / 门禁 | `ACCEPT` |

- DPO/RLHF：`{'on_fire': '拒合并或限量合并；po_risk0 Top-k 人工复核偏好对；可选 √PO 重加权', 'on_quiet': '放行合并，w=1'}`

## 怎么跑

```bash
PYTHONPATH=. python3 scripts/agod/llm_infer_align_prototype.py
# → results/agod/llm_infer_align/
# → docs/biz/LLM_INFER_ALIGN_PROTOTYPE.md
```

## 还有别的吗？

| 还要什么 | 包 / 文件 |
|----------|-----------|
| 流式一整段 hop 历史 | `run_rfperm_stream(Stream(...))` |
| 门控 √PO 权重 | `po_iptw_weights(..., mode='sqrt')` |
| 更重的 refit 流 | `agod.po_refit.run_refit_stream` |
| 接到客服¥账 | `scripts/agod/run_biz_value_sql_demo.py` |
| 方法×情景异同 | `docs/biz/METHOD_BIZ_SCENARIO_PROTOTYPE.md` |

边界：不是事实核查器；不是因果归因；fire 只打开对照窗 / 门禁。
