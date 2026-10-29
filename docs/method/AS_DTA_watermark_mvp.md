# Watermark evidence MVP — six real-data kinds

## Same vs Diff

### Same

- All kinds emit sparse evidence into the same strength enum
- Gray / low margin → vague → human (never hard strong-Y)
- strong_neg / missing → never write Y=human; hand off to PO+FN

### Diff

- api_log: strongest provenance; robust≈score (no wash)
- hard/soft_wm: p-value + robust probe; wash → gray
- fingerprint: never strong_pos alone; high FP under OOD
- user_declare: policy bit, not cryptographic recovery
- missing: structural absence, not a negative test result

## Counts by kind

| kind | n | strong_pos | gray | strong_neg | missing | vague_rate |
| --- | --- | --- | --- | --- | --- | --- |
| api_log | 80 | 21 | 29 | 27 | 3 | 0.45 |
| hard_wm | 80 | 26 | 16 | 35 | 3 | 0.30 |
| soft_wm | 80 | 14 | 33 | 27 | 6 | 0.51 |
| fingerprint | 80 | 0 | 42 | 32 | 6 | 0.57 |
| user_declare | 80 | 20 | 5 | 0 | 55 | 0.15 |
| missing | 80 | 0 | 0 | 0 | 80 | 0.14 |

![](watermark_mvp_kinds.png)

```bash
python3 scripts/run_as_dta_watermark_mvp.py
```
