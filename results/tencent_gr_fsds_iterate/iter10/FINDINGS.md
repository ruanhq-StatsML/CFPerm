# iter10 — DS handoff wrap (CLI + overnight summary)

Overnight loop **closed**. 讲武德: `W`=period; PO feeds FSDS — **not** ATE.

## Shipped

CLI:

```bash
PYTHONPATH=. python3 scripts/tencent_gr/po_help_fsds.py \
  --select-k 15 --seed 0 \
  --out-dir results/tencent_gr_fsds_iterate/iter10/cli_smoke_s0
```

Runs `po_help_select` → official FSDS; writes VIMP / blend / selected / `summary.json`.

## CLI smoke (cached grids)

| seed | W2 HGB AUC | AP | PO-risk |
|---:|---:|---:|---:|
| 0 | 0.759 | 0.0061 | 1.19e-6 |
| 1 | 0.529 | 0.0012 | 1.25e-6 |
| 2 | 0.871 | 0.0136 | 1.15e-6 |
| **mean** | **0.720** | | |

k=18 seed0: W2 **0.765** (slight bump vs k=15 on that seed).

## Locked recommendations

See `ITERATION_LOG.md` → **OVERNIGHT_SUMMARY**. Short form:

1. Default: **PO-help → FSDS @k=15** (CLI above).  
2. Optional k=18 for plain PO-VIMP; rare-π for lower σ; LOO-pos maj to freeze one list.  
3. Always seed mean±std. Never claim ATE.

20-min iterate timer **unsubscribed** after this wrap.
