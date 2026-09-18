# OnlineRFPerm + LLM routing --- FAR

n_reps=25, n_batches=20.
Routing: high-VIMP RF + low-VIMP TabPFN-style kNN.

| method | stationary FAR | random-noise FAR |
|---|---|---|
| addis | 0.0% | 0.0% |
| saffron | 0.0% | 0.0% |
| fix_alpha | 60.0% | 80.0% |
| hop | 96.0% | 92.0% |
| page_hinkley | 64.0% | 40.0% |
| ewma | 8.0% | 0.0% |
| cusum | 48.0% | 32.0% |
| ddm | 4.0% | 4.0% |
| adwin | 4.0% | 4.0% |
| martingale | 20.0% | 36.0% |
