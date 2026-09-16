# Fracture perturbations — detect inference breaks

> **一句（简单但实用）**：不火 ≈ 推理链路在探针下正常；中间突然 fracture / concept drift → OnlineRFPerm 给断裂打时间戳。

## Setup

| item | value |
|------|-------|
| n / n_ref / win / gate | 200 / 100 / 20 / 1.25 |
| faith | answer-precision → binary $Y$ |
| smooth control | all-quiet stream (no mid break) |
| fracture | sudden change at $t=n_{\mathrm{ref}}$ |

## Perturbations

| kind | what breaks |
|------|-------------|
| `invent_fracture` | serving break: invent + drop grounding |
| `label_flip` | concept drift on $Y$ |
| `answer_corrupt` | faithfulness break: scrambled answers |

## Results

Caught **18/18** fracture cases.

| dataset | perturbation | y quiet→frac | freeze delay | slide delay | smooth first1 | read |
|---------|--------------|--------------|-------------:|------------:|---------------|------|
| halueval | `invent_fracture` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| halueval | `label_flip` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| halueval | `answer_corrupt` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| squad | `invent_fracture` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| squad | `label_flip` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| squad | `answer_corrupt` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| hotpotqa | `invent_fracture` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| hotpotqa | `label_flip` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| hotpotqa | `answer_corrupt` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| truthfulqa | `invent_fracture` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| truthfulqa | `label_flip` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| truthfulqa | `answer_corrupt` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| boolq | `invent_fracture` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| boolq | `label_flip` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| boolq | `answer_corrupt` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| nq_open | `invent_fracture` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| nq_open | `label_flip` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |
| nq_open | `answer_corrupt` | 0.00→1.00 | 0 | 1 | None/None | fracture_caught_smooth_quiet — inference break timestamped |

## Read

- Smooth quiet + fracture fire → **can timestamp an inference break**.
- Your experience matches: under concept-drift-style flips the RF component
  usually fires; no fire on the quiet control means the chain looked fine.

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_fracture_perturb.py \
  --n 200 --n-ref 100 --win 20 --gate 1.25
```
