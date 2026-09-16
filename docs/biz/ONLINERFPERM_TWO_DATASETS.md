# OnlineRFPerm × two LLM infer datasets (CPU)

LaTeX: `docs/biz/ONLINERFPERM_TWO_DATASETS.tex` · `results/agod/online_rfperm_two_datasets/two_datasets.tex`

| dataset | delay | fire | y_bad | action |
|---|---|---|---|---|
| halueval | 0 | 3 | 0.19→0.31 | `model_rollback_or_audit_topk` |
| squad | 2 | 5 | 0.75→0.92 | `model_rollback_or_audit_topk` |

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_two_datasets.py
```
