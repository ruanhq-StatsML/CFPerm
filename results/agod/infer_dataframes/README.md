# Infer DataFrames — tidy serving form

Atomic row (same skeleton as HH preference stream):

```
t_idx | batch | question | answer | y
```

Audit / route (not required for fire): `dataset`, `hopped`, `rag_hit`, `faith`, `knowledge`, `gold`, `system`.

| file | role |
|------|------|
| `*_source.parquet` | `question / knowledge / gold` |
| `*_stream_table.parquet` | generated `answer` + $y$ + batch clock |

Datasets: HaluEval, SQuAD, HotpotQA, TruthfulQA.  
Clock in the note: $n_{\mathrm{per}}=20$, cut$=5$, $n_{\mathrm{ref}}=100$.

**Justify (formulation):** `docs/biz/ONLINERFPERM_DATA_FORM_JUSTIFY.md`  
**Full note:** `docs/biz/ONLINERFPERM_INFER_DRAFT.tex`  
**All-method table:** `docs/biz/ALL_METHODS_ONE_TABLE.tex`
