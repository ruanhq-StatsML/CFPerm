# Infer DataFrames

Tidy form (same spirit as HH `stream_table.parquet`):

`t_idx | batch | dataset | question | answer | y | hopped | rag_hit | faith | ...`

| file | role |
|------|------|
| `*_source.parquet` | question / knowledge / gold |
| `*_stream_table.parquet` | generated answer + $y$ + batch clock |

Datasets: HaluEval, SQuAD, HotpotQA, TruthfulQA.

Clock used in the note: $n_{\mathrm{per}}=20$, cut$=5$, $n_{\mathrm{ref}}=100$.

Write-up: `docs/biz/ONLINERFPERM_INFER_DRAFT.tex`  
All-method table: `docs/biz/ALL_METHODS_ONE_TABLE.tex`
