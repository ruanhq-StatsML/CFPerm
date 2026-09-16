# Multi-dataset OnlineRFPerm ($n_{\mathrm{ref}}=100$, answer-precision)

Clock: $n_{\mathrm{per}}=20$, cut$=5$, $n_{\mathrm{ref}}=100$. Backend: mock.  
Faith: **answer-precision** $|A\cap K|/|A|$（不是 Jaccard）。$Y\in\{0,1\}$。

| dataset | n | n_ref | fire | delay | $y$ quiet→hop |
|---|---:|---:|---:|---:|---|
| HaluEval | 200 | 100 | 5 | 0 | 0.00→1.00 |
| SQuAD | 200 | 100 | 5 | 0 | 0.00→1.00 |
| HotpotQA | 200 | 100 | 5 | 0 | 0.00→1.00 |
| TruthfulQA | 200 | 100 | 5 | 0 | 0.00→1.00 |

Responses / $Y$ / streaming note: `docs/biz/ONLINERFPERM_RESPONSES_AND_Y.md`  
Export: `data/hf_cache/infer_bench_export/` · streams: `results/agod/infer_dataframes/`
