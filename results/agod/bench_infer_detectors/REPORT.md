# Detector bench — manuscript `first_k` (trail batch)

Logic: `docs/biz/FIRST_K_LOGIC.py`  
Datasets: `results/agod/online_rfperm_multi_datasets` (export: `data/hf_cache/infer_bench_export`)  
Clock: $n_{\mathrm{per}}=20$, cut$=5$, $n_{\mathrm{ref}}=100$  
Narrative: `docs/biz/ONLINERFPERM_INFER_DRAFT.tex` · one-page table: `docs/biz/ALL_METHODS_ONE_TABLE.tex`

Index `0` = first trail batch after cut; `---` / `None` = never.

## Table 1 — first1 (onset)

| method | HaluEval | SQuAD | HotpotQA | TruthfulQA |
|---|---:|---:|---:|---:|
| OnlineRFPerm | 0 | 0 | --- | 0 |
| BOCPD | 0 | 0 | 4 | 0 |
| Page--Hinkley | 0 | 0 | --- | 0 |
| ADWIN | 0 | 0 | --- | 0 |
| DDM | 0 | 0 | --- | 0 |
| STEPD | 0 | 0 | 0 | 0 |
| HDDMA | 0 | 0 | --- | 0 |
| ECDDWT | 0 | 0 | --- | 0 |

## Table 2 — first2 / first3 (sustained)

| method | HaluEval k=2 | HaluEval k=3 | SQuAD k=2 | SQuAD k=3 | HotpotQA k=2 | HotpotQA k=3 | TruthfulQA k=2 | TruthfulQA k=3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OnlineRFPerm | --- | --- | --- | --- | --- | --- | --- | --- |
| BOCPD | --- | --- | --- | --- | --- | --- | --- | --- |
| Page--Hinkley | --- | --- | --- | --- | --- | --- | 0 | --- |
| ADWIN | --- | --- | --- | --- | --- | --- | --- | --- |
| DDM | 0 | 0 | 0 | 0 | --- | --- | 0 | 0 |
| STEPD | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HDDMA | --- | --- | --- | --- | --- | --- | --- | --- |
| ECDDWT | 0 | 0 | 0 | 0 | --- | --- | 0 | 0 |

## Read

- Onset ($k=1$): OnlineRFPerm aligns with selective baselines at index 0 on HaluEval / SQuAD / TruthfulQA.
- Sustained ($k=2,3$): OnlineRFPerm fires once at the hop then usually returns quiet; DDM / STEPD / ECDDWT stay hot.
- HotpotQA: quiet $y$ saturated under overlap labels — treat as negative control, not method ranking.

LaTeX: `docs/biz/BENCH_INFER_DETECTORS.tex`
