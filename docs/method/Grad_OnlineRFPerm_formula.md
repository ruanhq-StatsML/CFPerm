# Grad-OnlineRFPerm — formula writeup (not AGOD)

Source: [`docs/method/Grad_OnlineRFPerm_formula.tex`](Grad_OnlineRFPerm_formula.tex)  
PR: https://github.com/ruanhq-StatsML/CFPerm/pull/64

## Scalar gate (one stream)

Unfrozen set $\theta_U=\{\theta_i:\texttt{requires\_grad}\}$. One scalar:

$$
g_t=\bigl\|\nabla_{\theta_U} L_t\bigr\|_2
=\Biggl(\sum_{\theta_i\in\theta_U}\bigl\|\nabla_{\theta_i} L_t\bigr\|_2^{2}\Biggr)^{1/2}.
$$

**Not** a mean of per-layer norms; **not** per-layer OnlineRFPerm + `any(reject)`.

## OnlineRFPerm

$$
T_t=g_t-e_{\mathrm{ref}},
\quad
e_{\mathrm{ref}}=\tfrac1b\sum_{s<b}g_s.
$$

EWMA $p$-value (large $T$ → small $p$), then alpha-investing online FDR → one reject $R_t^{\mathrm{grad}}$.

## Diagnostic shares (no FDR)

$$
\pi_{t,\ell}=\frac{\|g_{t,\ell}\|_2}{g_t},\qquad
\sum_\ell\pi_{t,\ell}=1
$$

only after global reject (freeze-depth ranking).

## Lead

$$
\mathrm{Lead}=t_{\mathrm{grad}}-t_{\mathrm{MSE}}
\quad(\text{negative = Grad earlier}).
$$

MVP (5×5): mean Lead $=-3.24$, $P(\mathrm{Lead}<0)=72\%$.

## Supplementary experiments

See [`Grad_OnlineRFPerm_extras.md`](Grad_OnlineRFPerm_extras.md):

- Null / grace FPR
- α sensitivity
- Freeze closed-loop (MSE vs FLOPs)
- Extra stream packs (metro / beijing / stocks / waymo)
