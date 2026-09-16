# HotpotQA faithfulness prototype (OnlineRFPerm)

## Problem

Default multi-dataset path labels with **Jaccard(answer, full distractor)**.
Hotpot `knowledge` is ~10 paragraphs; quiet answers are short → Jaccard ≈ 0.07
→ quiet `y≈1` → no hop contrast → OnlineRFPerm does not fire.

## Handling logic

1. Load Hotpot `supporting_facts` → `knowledge_support` (cited sentences only).
2. Keep `knowledge_full` for the legacy baseline.
3. Quiet prompt uses **support**, not the full distractor dump.
4. Compare labelers on the same generations:

| labeler | rule |
|---|---|
| `y_jaccard_full` | Jaccard(A, K_full) < thr (legacy) |
| `y_jaccard_support` | Jaccard(A, K_support) < thr |
| `y_prec_support` | **recommended**: answer-precision \|A∩S\|/\|A\| < thr |
| `y_prec_and_gold` | audit: precision + gold-hit (can add quiet noise) |

## Results

backend=`mock`, n_ref=`100`

| labeler | y quiet→hop | fire | delay |
|---|---|---|---|
| `y_jaccard_full` | 1.00→1.00 | None | None |
| `y_jaccard_support` | 0.01→0.99 | 5 | 0 |
| `y_prec_support` | 0.00→1.00 | 5 | 0 |
| `y_prec_and_gold` | 0.00→1.00 | 5 | 0 |

### Score means (quiet / hop)

- jaccard_full: 0.079 → 0.041
- jaccard_support: 0.447 → 0.106
- prec_support: 0.918 → 0.281
- gold_hit: 1.00 → 0.03

## Read

- Legacy `jaccard_full` stays saturated (quiet≈hop≈1) — this is the documented bug.
- `y_prec_support` opens a quiet→hop gap and OnlineRFPerm fires at cut (delay 0).
- `y_prec_and_gold` is optional audit; gold-miss alone should not drive the ORF clock.
- Formulation unchanged: still `t_idx|batch|question|answer|y`; only how `y` / knowledge are built changes.

```bash
PYTHONPATH=. python3 scripts/agod/online_rfperm_hotpot_proto.py --backend mock
```

