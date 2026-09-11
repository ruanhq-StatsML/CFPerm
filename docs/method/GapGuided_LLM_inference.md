# Gap-guided reading budget

This line is a **two-sample reading problem** on \(X\). It is not TMLE, not PO-risk, and not causal reasoning. There is no treatment effect, no \(Y\)-intervening, no clever covariate \(H_m\) in the loop.

\(W\in\{0,1\}\) is a batch/domain tag. \(\hat e(x)\approx P(W=1\mid X=x)\). Human blocks \(B_m\) partition coordinates. Scores \(\pi\) and \(s_m(x)\) localize *which block carries the two-sample difference*. They do not say that \(W\) causes \(Y\), or that a channel is a mediator.

## π as a prior (shrinkage, not a causal prior)

Channel identity \(C\in\{1,\ldots,M\}\) is the unobserved index “which block localizes \(P_0\neq P_1\) on \(X\)”.

- \(\pi\) is estimated from the **training two-sample only**. It is a Dirichlet-mean style estimate of \(P(C=m\mid\text{domain pair})\). After Stage-1 it is frozen. It does not depend on the query \(x\).
- \(s_m(x)=|\hat e(x)-\hat e(x_{-m})|/\sum|\cdot|\) is a noisy, per-row measurement of the same \(C\) (LOMO on the propensity, still two-sample).
- Blend \(r=\lambda\pi+(1-\lambda)s\) is **linear shrinkage**. \(\lambda=1\) ignores this \(x\) and uses the batch localization. \(\lambda=0\) is a hard instance gate.

That is the only sense in which \(\pi\) is a prior: a population distribution over *which block to read*, updated by a noisy row-level localization. It is not \(P(Y\mid do(X))\), not a propensity of a treatment, not TMLE.

**Same expert** is the \(\lambda=1\) rule: \(\hat m=\arg\max_m\pi_m\), then every later query is read by the same specialist. **Instance expert** is \(\hat m(x)=\arg\max_m s_m(x)\). Those are different evaluands. Do not mix them.

## Abstain: you may not drop a block

Loss of a reading policy, stripped of LLM poetry:

\[
L(E,C)=\mathbf{1}\{C\notin E\}\,L_{\mathrm{miss}} + |E|\,L_{\mathrm{token}}
\]

If the posterior over \(C\) is flat, the action that avoids \(L_{\mathrm{miss}}\) is \(E=\{1,\ldots,M\}\): **do not subset**. That is abstain-from-dropping.

Two operationalizations, same logic:

| if you must produce a score | if you may refuse |
|---|---|
| **pack all** — concatenate every \(B_m\) (or enable every tool) | **HITL** — do not pack a subset; ask for the top channel |

Abstain is **not** “the LLM is unsure about the answer”. It is “the two-sample gap is not localized, so a subset window is an unjustified restriction of the σ-algebra”.

Gate (on \(\pi\) for the same template; on \(r\) if you allow instance updates): entropy \(\ge 0.92\log M\), or \(\max\pi<\tau_{\mathrm{lo}}\). When \(\max\pi\ge\tau_{\mathrm{hi}}\), \(|E|=1\) (same expert). In between, \(|E|=k\) concat.

Forced top-1 on a **diffuse** shift (equal signal in every block) is the failure mode abstain exists for: you discard three-quarters of the two-sample signal.

## Context packing and concatenating enabled blocks

Opinion pool \(\sum_m\pi_m\hat e_m(X_{B_m})\) **looks at every specialist’s score**. That is not packing.

Context packing is one reader on a **subset of raw channels**:

\[
E=E(\pi)\qquad\text{(same template for every }x\text{)}
\]
\[
\mathrm{pack}(x)=\mathrm{concat}_{m\in E}^{\text{read-order}}\;\mathrm{serialize}(x_{B_m})
\]

Rules:

1. Order is read-order (descending \(\pi\) or \(r\)). Concatenation is the prompt string; the statistical analogue is column-bind of those slices, **one** RF/LLM on the packed design.
2. Disabled blocks are **absent**. Do not emit empty `### text` headers (that still spends tokens and leaks that the channel exists).
3. Sidecar \(\pi,s,r,Z\) may sit in the system prompt. \(W\) does not. \(H_m\) does not (it is not part of this line).
4. \(|E|=1\) packing **is** the same expert, written as a window: the user message contains only that section.
5. \(|E|=k>1\) is still **one** forward pass, several sections. It is not \(k\) experts voting.

`pack_user_context` does (4)–(5) as strings. `pack_blocks` + `holdout_packed_auc` do the same map on columns.

## “同一个专家” — characterize, then evaluate

Let \(\mathcal{G}_m\) be functions of \(X_{B_m}\) only (measurable w.r.t. that block). Training data may be used; **this** \(x\)’s other blocks may not.

**Same expert.** Choose \(\hat m=\hat m(\pi)\) from the training two-sample. Deploy one \(f\in\mathcal{G}_{\hat m}\) for every subsequent query. Equivalent packing: \(E=\{\hat m\}\) for the whole deployment window.

Not the same object:

| object | \(m\) depends on | what the reader sees |
|---|---|---|
| same expert | \(\pi\) only | \(X_{B_{\hat m}}\) for all \(x\) |
| instance expert | \(\pi,s(x)\) | \(X_{B_{\hat m(x)}}\) — different specialist per row |
| concat-\(k\) | \(\pi\) only | \(\mathrm{concat}_{m\in E} X_{B_m}\), \(|E|=k\) |
| pack all / abstain | — | full \(X\) |
| opinion pool | \(\pi\) | all \(\hat e_m\), not raw concat |

Evaluate four numbers, in order:

1. **Localization.** \(1\{\hat m=m^\star\}\) for same-expert π vs VIMP vs a *single* random specialist (not a per-row coin flip — that would not be 同一个). Instance hit \(P(\hat m(X)=m^\star)\) is a separate line.
2. **Specialist utility.** Held-out two-sample AUC of the frozen \(\hat e_{\hat m}\). Ceiling: oracle \(m^\star\). This is “hire one specialist”.
3. **Packed-reader utility.** Held-out AUC of **one** RF trained on `pack_blocks(X, E)`, \(E\in\{\{\hat m_\pi\},\{\hat m_{\mathrm{VIMP}}\},\{m^\star\},\mathrm{top}\text{-}2(\pi),\mathrm{all}\}\). This is “one reader, concatenated sections”. It is the eval for packing, not (2).
4. **Abstain check.** Diffuse / shuffle: pack-all \(\ge\) forced same-expert. Concentrated inject: same-expert π \(\approx\) oracle and **beats** VIMP when a wide nuisance block exists. Negative control: shuffle \(m^\star\) → hit on that name collapses.

Headline from `scripts/run_gap_guided_inference.py` (3 seeds):

- Inject valence, wide text: same-expert π packed AUC **0.753** (hit 1); VIMP packed AUC **0.600** (hit 0).
- Diffuse equal shift: pack-all **0.876** vs one specialist **0.710**. Abstain-from-drop is pack-all here.
- Shuffle valence: packed AUC ~0.5, hit collapses.

## Code

```bash
python3 -m pytest -q tests/test_gap_guided_inference.py
python3 scripts/run_gap_guided_inference.py
```

`population_pack` → \(E(\pi)\). `pack_user_context` → string concat. `holdout_packed_auc` → one reader on packed columns.
