# Prediction table — one Y, many X

Not attribution. Not VIMP. A supervised table:

- `y`: binary audit decision (the thing to predict)
- `x_*`: features of the reply (the predictors)
- `batch`: arrival window

OnlineRFPerm: `fit_online_probe(X, y)` predicts y from x on the last window,
then checks whether that map still predicts the next window.

Files:

- `results/agod/llm_audit_consistency/xy_hh_helpful_consistent.csv`  n=1200  regime=consistent
- `results/agod/llm_audit_consistency/xy_hh_helpful_hop.csv`  n=1200  regime=hop
- `results/agod/llm_audit_consistency/xy_hh_harmless_consistent.csv`  n=1200  regime=consistent
- `results/agod/llm_audit_consistency/xy_hh_harmless_hop.csv`  n=1200  regime=hop

Schema: `y,batch,x_n_toks,x_n_chars,...,x_thank`

