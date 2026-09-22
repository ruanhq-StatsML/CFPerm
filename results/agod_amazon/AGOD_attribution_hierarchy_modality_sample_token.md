# Attribution grains: modality / sample / patch-token

This note separates **what each grain measures** from **what the Amazon
prototype actually controls**. Only the modality grain is closed-loop.

## 1. Three grains = three questions

| Grain | Question | Statistic | Matching action |
|---|---|---|---|
| **Modality** | Which channel drifts? | MSG \(g_m\) on block \(X_m\) | **LR\(_m\) / gate \(\partial L/\partial\theta_m\)** ← implemented |
| **Sample** | Which examples carry the shift? | domain score / residual on row \(i\) | sample weight / replay — *not* LR\(_m\) |
| **Patch / token** | Where inside a modality? | LOO / occlusion \(\Delta_k\) | patch mask / token distill weight — *not* LR\(_m\) |

## 2. Upstream vs downstream (different jobs)

- **Upstream** (shift diagnostics → control): AUC\(_m\), VIMP, PO, \(g_m\), \(\alpha_m\), LR×, adapt FLOPs.
- **Downstream** (task utility): holdout Acc/AUC Δ after the update; reference Acc Δ (forgetting).

They need not move together. High \(g_{\text{image}}\) does not guarantee larger holdout ΔAcc.

## 3. Modality layer (implemented control law)

```
g_m = Normalize(AUC_m · VIMP_m + γ · PO_m)
α   = Softmax(g / τ)          # EMA
LR_m = lr0 · (β + (1-β) · α_m · |M|)
if α_m < θ: zero grads of modality-m projection
```

## 4. Sample layer (logic only)

On a fixed modality representation (or concat):
- domain RF: reference vs current;
- score \(s_i = P(W{=}1\mid x_i)-1/2\) (or residual risk).
High \(s_i\) = row looks most out-of-reference.

**Do not** fold \(s_i\) into Softmax→LR\(_m\). That mixes grains.
If ever wired, action = \(w_i=w(s_i)\) inside \(\sum_i w_i \ell_i\) or replay buffer.

## 5. Patch / token layer (logic only)

Inside the already-chosen modality:
- leave-one-patch-out / token occlusion → \(\Delta_k\) of domain score or PO;
- rank / visualize top patches or tokens.

**Do not** re-decide LR\(_m\). If ever wired, action = mask / reweight that
modality's local distill terms.

## 6. Why this prototype stops at modality

The actuator is already \(\alpha_m\to\mathrm{LR}_m\).
Sample/patch scores answer finer localization questions; without a matching
actuator they are diagnostics, not “secondary AGOD”.
