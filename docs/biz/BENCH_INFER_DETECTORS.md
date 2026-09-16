# Detector bench — manuscript first_k (trail batch)

Logic: `docs/biz/FIRST_K_LOGIC.py`  
Faith: answer-precision → binary $Y$  
Clock: $n_{\mathrm{per}}=20$, cut$=5$, $n_{\mathrm{ref}}=100$  
Narrative: `docs/biz/ONLINERFPERM_RESPONSES_AND_Y.md` · `ONLINERFPERM_INFER_DRAFT.tex`

Index `0` = first trail batch after cut; `---` / `None` = never.

## Table 1 — first1 (onset)

| method | HaluEval | SQuAD | HotpotQA | TruthfulQA |
|---|---:|---:|---:|---:|
| OnlineRFPerm | 0 | 0 | 0 | 0 |
| BOCPD | 0 | 0 | 0 | 0 |
| Page--Hinkley | 0 | 0 | 0 | 0 |
| ADWIN | 0 | 0 | 0 | 0 |
| DDM | 0 | 0 | 0 | 0 |
| STEPD | 0 | 0 | 0 | 0 |
| HDDMA | 0 | 0 | 0 | 0 |
| ECDDWT | 0 | 0 | 0 | 0 |

## Table 2 — first2 / first3 (sustained)

| method | HaluEval k=2 | HaluEval k=3 | SQuAD k=2 | SQuAD k=3 | HotpotQA k=2 | HotpotQA k=3 | TruthfulQA k=2 | TruthfulQA k=3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OnlineRFPerm | --- | --- | --- | --- | --- | --- | --- | --- |
| BOCPD | --- | --- | --- | --- | --- | --- | --- | --- |
| Page--Hinkley | --- | --- | --- | --- | --- | --- | --- | --- |
| ADWIN | --- | --- | --- | --- | --- | --- | --- | --- |
| DDM | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| STEPD | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HDDMA | --- | --- | --- | --- | --- | --- | --- | --- |
| ECDDWT | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Under answer-precision, Hotpot onset matches the other three. Default clock remains $n_{\mathrm{per}}=20$ (not streaming $n_{\mathrm{per}}=1$).

LaTeX: `docs/biz/BENCH_INFER_DETECTORS.tex` · one-page: `docs/biz/ALL_METHODS_ONE_TABLE.tex`
