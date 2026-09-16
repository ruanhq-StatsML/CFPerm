# Detector bench — manuscript first_k (trail batch)

Logic: `docs/biz/FIRST_K_LOGIC.py`

Datasets: `/workspace/data/hf_cache/infer_bench_export`

## Table 1 — first1 (trail batch)

| method | halueval | squad |
|---|---:|---:|
| OnlineRFPerm | 0 | 2 |
| BOCPD | 0 | 1 |
| PageHinkley | 0 | 1 |
| ADWIN | None | None |
| DDM | 0 | 1 |
| STEPD | 2 | None |
| HDDMA | None | None |
| ECDDWT | 0 | 1 |

## Table 2 — first2 / first3

| method | halueval k=2 | halueval k=3 | squad k=2 | squad k=3 |
|---|---:|---:|---:|---:|
| OnlineRFPerm | None | None | None | None |
| BOCPD | None | None | 1 | None |
| PageHinkley | 0 | None | 1 | None |
| ADWIN | None | None | None | None |
| DDM | 0 | 0 | 1 | None |
| STEPD | None | None | None | None |
| HDDMA | None | None | None | None |
| ECDDWT | 0 | None | 1 | None |

LaTeX: `docs/biz/BENCH_INFER_DETECTORS.tex`

