# One simplex, several batch procedures

\(\pi\) is a frozen two-sample localization over human blocks \(B_m\). It is **not** an online-learning iterate. Nothing here updates \(\pi\) or the reader from a query stream, a loss sequence, or a regret bound.

The same simplex is a **bundle** of statistical maps. They differ in which object \(\pi\) multiplies.

## The bundle

| map | what \(\pi\) multiplies | finite-sample object |
|---|---|---|
| same-expert packing | the observation: \(E=\{\arg\max\pi\}\) | σ-algebra of one reader |
| concat packing | observation: top-\(k(\pi)\) raw blocks | one reader, several sections |
| query update of \(E\) | \(r=\lambda\pi+(1-\lambda)s(x)\) | still \(E\), not \(\pi\) |
| opinion pool | specialist *scores* \(\sum\pi_m\hat e_m\) | not packing |
| π-RF / BAWF | feature *sampling* \(p_j\propto\pi_m/|B_m|\) | trees (scale-invariant) |
| adaptive group ridge | penalties \(\lambda_m\propto 1/\pi_m\) | linear / GD equivalent of column scale |
| **next-batch group LR** | steps \(\eta_m=\eta_0 M\pi_m\) | finite-step GD on the *next* batch |
| **block boosting round** | next weak learner lives in \(B_{\hat m}\) | gradient boosting, still batch |

All of them are mean-preserving reallocations of a budget (tokens, splits, ridge, steps). Average effort stays \(\eta_0\); only *where* it is spent changes.

## Next-batch learning rates (the GD twin)

Trees ignore monotone column scaling, so “feature weights” do nothing to RF splits unless you change sampling (π-RF). Differentiable / linear models *do* feel a step size.

Batch A: estimate \(\pi\) from \(W\). Freeze it.
Batch B: the next batch, not a stream.
\[
\theta_m \leftarrow \theta_m - \eta_m\nabla_{\theta_m}L_B,\qquad
\eta_m=\eta_0 M\,\tilde\pi_m.
\]
High \(\pi_m\): that block drifted, so \(\theta_m\) is more stale → larger step on B.
Low \(\pi_m\): leave it slow. Damp \(\eta_m\propto 1/\pi_m\) is the ablation (treat the drifted block as nuisance).

Why finite steps: at convergence of unregularized ERM the allocation of \(\eta\) washes out. The statistic is the **path**. Early stopping / ridge makes the preconditioner \(\mathrm{diag}(\eta_m)\) identifiable.

Same math as adaptive ridge: boosting \(\eta_m\) is like shrinking \(\lambda_m\). Same simplex as packing: packing *drops* low-\(\pi\) raw blocks; group LR *slows* them but still sees them. Choose packing when there is a token budget; choose LR when the next batch still trains on all coordinates.

## What else is sitting in the same drawer

- **Block gradient boosting:** next round’s base learner is fit only on \(B_{\hat m}\) (or sampled with \(p\propto\pi\)). Friedman GBM, not AdaBoost-online.
- **Two-timescale batch EM:** slow \(\theta_{\mathrm{low}\,\pi}\), fast \(\theta_{\mathrm{high}\,\pi}\), one pass on B.
- **Block-diagonal Fisher / Newton:** \(H_m^{-1}\) scaled by \(\pi_m\) as a batch preconditioner.
- **Proximal pull to the previous fit:** \(L_B+\sum_m(1-\pi_m)\|\theta_m-\theta_m^{(A)}\|^2\) — stable blocks stay put.
- **Evaluation** is always A/B/C slices: localization on A, train on B, holdout on C. Switch-rate / path AUC / mass on GT. Shuffle \(m^\star\) as negative control. VIMP-boost as the width-biased competitor.

None of that is a bandit over experts.

Finite-step GD on the next batch (3 seeds): synthetic valence, holdout AUC uniform **0.927** vs π-boost **0.943** vs damp **0.751**; coefficient mass on GT 0.46 → **0.67**. Inject: VIMP-boost 0.783 vs π-boost 0.805; damp 0.738.

```bash
python3 -m pytest -q tests/test_pi_next_batch.py
python3 scripts/run_pi_next_batch.py
```
