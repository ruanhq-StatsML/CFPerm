# Multi-dataset OnlineRFPerm ($n_{\mathrm{ref}}=100$)

Clock: $n_{\mathrm{per}}=20$, cut$=5$, $n_{\mathrm{ref}}=100$. Backend: mock scaffold.  
Write-up: `docs/biz/ONLINERFPERM_INFER_DRAFT.tex` · all-method table: `docs/biz/ALL_METHODS_ONE_TABLE.tex`

| dataset | n | n_ref | fire | delay | $y_{\mathrm{bad}}$ quiet→hop |
|---|---:|---:|---:|---:|---|
| HaluEval | 200 | 100 | 5 | 0 | 0.00→1.00 |
| SQuAD | 200 | 100 | 5 | 0 | 0.00→1.00 |
| HotpotQA | 200 | 100 | --- | --- | 1.00→1.00 |
| TruthfulQA | 200 | 100 | 5 | 0 | 0.00→0.98 |

HotpotQA: quiet already saturated under overlap labels — negative control, not method ranking.

Export: `data/hf_cache/infer_bench_export/` · streams: `results/agod/infer_dataframes/`
