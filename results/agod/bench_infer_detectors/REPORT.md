# Detector bench — manuscript first_k (trail batch)

Logic: `docs/biz/FIRST_K_LOGIC.py`

Datasets: `/workspace/data/hf_cache/infer_bench_export`

## Table 1 — first1 (trail batch)

| method | halueval | squad | hotpotqa | truthfulqa |
|---|---:|---:|---:|---:|
| OnlineRFPerm | 0 | 0 | None | 0 |
| BOCPD | 0 | 0 | 4 | 0 |
| PageHinkley | 0 | 0 | None | 0 |
| ADWIN | 0 | 0 | None | 0 |
| DDM | 0 | 0 | None | 0 |
| STEPD | 0 | 0 | 0 | 0 |
| HDDMA | 0 | 0 | None | 0 |
| ECDDWT | 0 | 0 | None | 0 |

## Table 2 — first2 / first3

| method | halueval k=2 | halueval k=3 | squad k=2 | squad k=3 | hotpotqa k=2 | hotpotqa k=3 | truthfulqa k=2 | truthfulqa k=3 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OnlineRFPerm | None | None | None | None | None | None | None | None |
| BOCPD | None | None | None | None | None | None | None | None |
| PageHinkley | None | None | None | None | None | None | 0 | None |
| ADWIN | None | None | None | None | None | None | None | None |
| DDM | 0 | 0 | 0 | 0 | None | None | 0 | 0 |
| STEPD | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HDDMA | None | None | None | None | None | None | None | None |
| ECDDWT | 0 | 0 | 0 | 0 | None | None | 0 | 0 |

LaTeX: `docs/biz/BENCH_INFER_DETECTORS.tex`

