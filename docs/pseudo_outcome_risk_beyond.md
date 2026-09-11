# Adapt for Distribution Shift, Then Predict the Next Batch

Nothing new. This is not online learning. It is the ordinary shift-adaptation loop: from batch $t$ localize where $D_t$ differs from $D_0$, adapt **only there**, use that to predict batch $t{+}1$.

The LLM (or any foundation encoder) is a frozen generator $p_\theta$. Pseudo-outcome (PO) risk is the discrepancy; VIMP / LOCO of $\hat\tau$ names the subset $S^\star$ that is worth the adaptation budget. Then you predict on the next batch. Unique CS / CD explanation is still impossible; subset localization is still all you need.

No prototype in this pass.

---

## 1. The loop (the whole idea)

\[
D_0 \;\xrightarrow{\ \varphi,\;\hat\tau,\;\mathrm{VIMP}\ }\; S^\star
\;\xrightarrow{\ \text{adapt on }S^\star\text{ only}\ }\;
\widehat{Y}_{t+1}.
\]

1. Score the gap between a reference batch $D_0$ and the current batch $D_t$ with the DR / R-learner residual $\varphi$ (distance, not CATE).
2. Localize: rank modalities / rows / spans by how much they carry $\hat\tau$. That ranking is $S^\star$.
3. Adapt for the shift **on $S^\star$** — reweight, residual on those coordinates, extra labels, retrieve those spans, re-encode that modality. Do not update $\theta$.
4. Predict the next batch with frozen $p_\theta$ composed with that localized adaptation.

That is domain adaptation / batch-shift correction with an identified *where*, not a new learning paradigm. Online learning would update $\theta_t$ as a martingale of losses. Here $\theta$ is fixed; what changes is a small correction aimed at the next batch.

---

## 2. Why localization is the only extra sentence

CFPerm already says you cannot uniquely decompose a performance drop into CS vs CD vs noise. So “adapt for distribution shift” cannot mean “retrain the generator until the story fits.” It means: **find a subset that accounts for the discrepancy, correct there, predict ahead.**

\[
S^\star
\;\in\;
\arg\max_{S:\;\mathrm{cost}(S)\le B}
\;
\mathrm{VIMP}_{\hat\tau}(S).
\]

$B$ is the adaptation budget for the next batch (compute, labels, a few residual coordinates). $S^\star$ is which part of $x$ that budget should hit. Layer 1 = modality, Layer 2 = instances / windows, Layer 3 = tokens / boxes. Same ranking, three resolutions, same use: better $\widehat{Y}_{t+1}$.

$\varphi=(Y-\hat\mu)(W-\hat e)$ (R-learner; DR is the efficient sibling) is the usual orthogonal residual. $\hat\mu$ and $\hat e$ partial out pooled outcome and propensity; $\hat\tau(x)$ puts the leftover gap back on $x$. Coordinates that still move $\hat\tau$ are the ones that will still hurt the next batch if you ignore them. That is the statistical reason the ranking is a map for *next-batch* adaptation rather than a heatmap of $p_\theta$.

---

## 3. Frozen generator is the default in this loop

$p_\theta$ was fit on some old mixture. The next batch is a different mixture. Classical response: importance weight, residual corrector, re-estimate a small head — not re-estimate the whole conditional.

If you unfreeze $\theta$ you are doing a different problem (continual / online learning). You also lose the CFPerm object: $\varphi$ is no longer a discrepancy *beside* a frozen predictor, orthogonality mixes with the generator’s optimizer, and a permutation audit no longer tests the same thing. None of that is required for “adapt on $S^\star$, predict $t{+}1$.”

Recomputing $\varphi$ when batch $t{+}1$ becomes batch $t{+}2$ is the same diagnostic on the next window. Still not online learning.

---

## 4. What we are willing to claim

- PO-risk localization is how this repo points an adaptation budget at a shift, so the next batch is cheaper to predict.
- The generator stays frozen. “Adapt” means a correction supported on $S^\star$, not $\theta\leftarrow\theta-\eta\nabla\theta$.
- There is no new theorem here and no unique CS/CD story. Subset localization, then next-batch prediction.

What we are not claiming: online learning; a new LLM training method; causal effects under confounding; that $p_\theta$’s own saliency is a substitute for VIMP of $\hat\tau$.
