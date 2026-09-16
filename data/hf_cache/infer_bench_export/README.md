# Infer datasets (n_ref >= 100)

| dataset | file | knowledge |
|---|---|---|
| halueval | `halueval.jsonl` | HaluEval knowledge |
| squad | `squad.jsonl` | SQuAD context |
| hotpotqa | `hotpotqa.jsonl` | Hotpot paragraphs |
| truthfulqa | `truthfulqa.jsonl` | correct answers |

Recommended clock: `n_per=20`, `cut_batch=5` → **n_ref=100**, `n_batches=10` → trail=100.

```bash
PYTHONPATH=. python3 scripts/agod/prepare_infer_datasets.py
PYTHONPATH=. python3 scripts/agod/online_rfperm_multi_datasets.py --n-per 20 --cut-batch 5 --n-batches 10
```
