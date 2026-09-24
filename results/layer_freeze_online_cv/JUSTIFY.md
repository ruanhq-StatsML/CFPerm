# Justify: OnlineRFPerm onset × PO / MSE / MMD vs ref

WHEN = frozen `RandomForestRegressor().predict(X_new)`, T = MSE − E_ref, last-two hop.
WHAT = PO-risk (RF nuisances) × serving MSE × MMD²(X_new, X_ref).
Concept DGP: P(X) fixed, β rotates after labeled onset → MMD stays quiet, MSE/PO move.
Covariate DGP: P(Y|X) fixed, μ(X) walks → MMD fires, PO stays quieter; MSE-only is X shift.

| dataset | n_new | onset_true | onset_hat | onset_rank | board | n_watch | n_x_shift | n_tricky | n_freeze | n_hop |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| dgp_concept | 2000 | 2 | None | None | watch | 5 | 0 | 0 | 0 | 0 |
| dgp_covariate | 2000 | 2 | None | None | keep_training | 0 | 0 | 0 | 0 | 0 |
| bankmarketing | 5000 | None | 4 | None | x_shift | 0 | 2 | 0 | 0 | 1 |
| eeg | 1500 | None | 2 | None | x_shift | 0 | 3 | 0 | 0 | 1 |

Freeze only when PO and MSE both break. No online-bootstrap.

