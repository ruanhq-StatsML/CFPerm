# Frozen Generator, Localization, and Where the Adaptation Budget Goes

**This is not online learning.** The LLM stays frozen. Pseudo-outcome (PO) risk plus variable importance tell you **which part of $x$ is most worth adapting** — which modality, which instances, which spans — under a scarce intervention budget. That is the statistical intuition. The rest of this note is only that sentence, unpacked.

No prototype in this pass.

---

## 0. Yes: “which part is most worth adapting”

That reading is the right one.

CFPerm’s doctrine is already: unique CS / CD / dependence decomposition is impossible; **subset localization is all you need.** PO-risk is the discrepancy. Localization of $\hat\tau$ is a ranking of subsets by how much they carry that discrepancy. If you are allowed to change *one* thing — re-encode one modality, retrieve on one span, collect labels on one slice, attach a tiny residual on one block of coordinates — you should change the subset that ranks first.

That is not:

- online learning of the generator,
- continual pretraining,
- LoRA on the 70B every $N$ users,
- “the model adapts itself because it saw new tokens.”

Those are *parameter updates of a huge conditional*. The object here is *an identified ranking of subsets*. The generator can stay frozen; the ranking tells the *policy around the generator* where to spend the budget.

---

## 1. Two objects, not one

Write the deployed system as

\[
\underbrace{p_\theta(y\mid x)}_{\text{frozen generator}}
\qquad\text{and}\qquad
\underbrace{\mathcal{L}(\varphi, \hat\tau; D_0, D_t)}_{\text{localization of the discrepancy}}.
\]

- $p_\theta$ is the LLM (or ViViT / CLAP / GPT-2 encoders). $\theta$ does **not** move. It is a fixed conditional, trained on some mixture that is already gone.
- $D_0$ is a reference batch (calibration users, trusted corpus, last week, no-tool). $D_t$ is the current batch (this user, this session, this prompt, this video window). $W$ is the batch label.
- $\varphi$ is the DR / R-learner residual — a Neyman-orthogonal score for “$D_t$ is distinguishable from $D_0$ after partialling out $\mu$ and $e$.” It is a **distance**, not a CATE (ignorability fails under the alternative; README).
- $\hat\tau(x)$ is $\varphi$ regressed on $x$. Variable importance / LOCO / permuCATE of $\hat\tau$ is a **ranking of coordinates, groups, instances, or spans**.

The statistical job is the second object. The first object is a sampler you are allowed to call.

If you update $\theta$ online, you have thrown away the split: the generator and the discrepancy become the same high-capacity fit, the ranking is no longer a readout of a frozen conditional against a moving world, and you are back to ordinary continual learning with no orthogonality and no localization theorem. That is someone else’s bill.

---

## 2. Why localization, not explanation, and not a full adapt

Three facts from the CFPerm impossibility, reused as inference logic:

1. A drop in performance (or a rise in PO-risk) has **infinitely many** CS / CD / noise realizations. “The concept drifted” is not identified. Adapting the whole generator is a response to an unidentified story: $\theta$ can absorb any of those realizations. You will not know *what* you adapted.
2. A change in VIMP likewise has infinitely many coordinate-plus-dependence realizations. So even the ranking is not a unique “root cause.” It *is* a **sufficient statistic for a budgeted intervention**: which subset, if you were allowed to touch only that subset, would move the discrepancy most.
3. Therefore the honest action is not “fine-tune everything until the loss comes back.” It is: **localize, then spend the adaptation budget on the localized set, leave $p_\theta$ alone.**

“Most worth adapting” is an optimization sentence with a constraint:

\[
S^\star
\;\in\;
\arg\max_{S:\; \mathrm{cost}(S)\le B}
\;
\text{VIMP}_{\hat\tau}(S)
\quad\text{or}\quad
\text{LOCO-risk of }S.
\]

$S$ can be a modality block (video vs audio vs text), a set of windows, a set of tokens, a bounding box, a retrieved passage vs a continuation. $B$ is compute, labels, latency, or “how many coordinates of a residual head you are willing to fit.” The generator is not in the argmax. It is the thing you generate *from* after the intervention on $S$.

That is the same logic as post-hoc subgroup localization in the batch test. The LLM does not change the math. It only makes the cost of adapting the *wrong* $S$ (or adapting $\theta$ instead of $S$) obvious.

---

## 3. What the residual actually says (statistical intuition)

Nuisances: $\hat\mu(x)\approx\mathbb{E}[Y\mid X=x]$, $\hat e(x)\approx\mathbb{P}(W=1\mid X=x)$. Cross-fit, clip $e$.

R-learner product (what the code stores as `residual_y * residual_t`):

\[
\varphi
=
\bigl(Y-\hat\mu(X)\bigr)\,
\bigl(W-\hat e(X)\bigr).
\]

DR / AIPW is the efficient sibling; same intuition. Then $\hat\tau(x)$ fits $\varphi$ on $X$, and PO-risk is $\mathbb{E}[\hat\tau(X)^2]$ or $\mathbb{E}[\varphi^2]$.

Read this without causal language:

- $\hat\mu$ partials out “what $Y$ looks like as a function of $x$ in the pooled data.”
- $\hat e$ partials out “how batch membership is already predictable from $x$.”
- Whatever is left in $\varphi$ is **discrepancy that is not already explained by the pooled outcome surface or by easy propensity**.
- $\hat\tau(x)$ puts that leftover *back onto $x$*. Coordinates where $\hat\tau$ is sensitive are coordinates where the two batches still disagree.

So localization is not saliency of the LLM, and it is not “this token caused the answer.” It is: **after you have removed the parts of $Y$ and $W$ that a pooled nuisance already knows, which parts of $x$ still carry the batch gap.** Those are the parts that are worth an intervention, because an intervention anywhere else is, to first order, already in $\mu$ or $e$.

Orthogonality is why this is a statistical method rather than a heuristic heatmap. First-order errors in one nuisance do not wreck $\mathbb{E}[\varphi]$ if the other is consistent. That is the only reason you can keep $p_\theta$ frozen and still trust a ranking computed from small, refittable nuisances sitting *beside* the generator.

---

## 4. Frozen generator + three layers = three scales of “where to adapt”

The generator proposes. Localization ranks. A tiny policy spends $B$ on $S^\star$. Nothing in $\theta$ moves.

```
frozen p_θ(y | x)
        │
        │  call / decode / retrieve if asked
        ▼
   current x  vs  reference D0
        │
        │  φ, τ̂, VIMP / LOCO
        ▼
   S* = most expensive-to-ignore subset
        │
        ├── Layer 1  modality / block     (video ⊕ audio ⊕ text)
        ├── Layer 2  instance / window    (which clip, which turn)
        └── Layer 3  token / box          (which span to fetch or mask)
        │
        ▼
   spend budget B on S* only
   (re-encode, retrieve, label, residual-on-S*, extra verify)
   then call p_θ again — still frozen
```

**Layer 1 — which block of $x$.**
On a concat embedding, VIMP shares say whether the gap sits in video, audio, or text. “Audio is the share that moves PO-risk” means: if you can adapt only one encoder, or run only one expensive encoder, or collect only one extra annotation stream, adapt **audio**. It does *not* mean fine-tune the LLM, and it does *not* mean “the concept of the video drifted.”

**Layer 2 — which rows.**
High $|\varphi|$ or high $|\hat\tau(x)|$ rows are the instances that carry the gap. If you can label, retrieve, or inspect only $k$ examples, take those. That is orthogonal-score targeting, not uncertainty sampling of $p_\theta$. The generator’s entropy can be high on a region that is already well represented in $D_0$; $\varphi$ is high where $D_t$ is still distinguishable from $D_0$. Those are different sets. The second set is the one worth adapting to.

**Layer 3 — which span.**
Map the same $\hat\tau$ back onto tokens or boxes. The spans that rank are the spans worth a local intervention: retrieve a grounding passage for *that* span, mask *that* span, ask a clarifying question about *that* mention, spend verify compute on *that* continuation. Integrated gradients of $p_\theta$ will generally disagree. They answer “what did the frozen generator use.” Localization answers “what still distinguishes the batches.” You adapt to the second, because the first is a property of $\theta$, which you have chosen not to touch.

Across layers the sentence is the same: **$S^\star$ is the cheapest sufficient intervention.** The LLM is the thing you call after you have intervened on $S^\star$.

---

## 5. What “adapt” is allowed to mean

Because $p_\theta$ is frozen, “adapt” cannot mean $\theta\leftarrow\theta-\eta\nabla_\theta$. It means an intervention whose support is $S^\star$:

| Budget $B$ | Intervention on $S^\star$ | What you do *not* do |
|---|---|---|
| Compute | Run the expensive encoder / verifier only on the localized modality or span | Re-encode everything; extra CoT on every token |
| Retrieval | Fetch docs / tools targeted at localized spans | Always-on RAG; LLM-as-router with no score |
| Labels | Annotate high-$\|\varphi\|$ rows, or the localized modality | Uniform relabel; full RLHF on the new user |
| Parameters | A residual / gate on the coordinates of $S^\star$ only | LoRA / full FT of the generator |
| Attention of the user | Highlight the localized span (“this is where the session moved”) | Ask the model to narrate a unique root cause |

The residual-on-$S^\star$ row is the only place a parameter moves, and it is not the LLM. It is a small regression on an orthogonal target, restricted to the localized coordinates. That is ordinary statistics. Calling it “online learning of the LLM” is a category error: the learning, if any, is of $\hat\tau\mid_{S^\star}$, which is a localization readout, not a new generator.

---

## 6. Why a frozen generator is the point, not a limitation

If you unfreeze $\theta$, three bad things happen at once.

1. **Identification.** The generator can fit any member of the impossibility class. You lose the right to say “we adapted the part that carried the discrepancy.” You adapted *a* interpolating conditional.
2. **Orthogonality.** Nuisances $\mu,e$ and the generator become entangled. $\varphi$ stops being a score *beside* the model and becomes another training loss. You no longer have a control variate; you have a second optimizer.
3. **Audit.** Permutation / CF-split type-I control is a statement about a discrepancy given a frozen (or at least cross-fit, not jointly updated) procedure. If $\theta$ is chasing $D_t$ online, the permutation no longer tests the same object.

Keeping $p_\theta$ frozen is what makes localization a *statistical* answer to “where to adapt.” The generator is a sampler with a fixed law. The world ($W$, the user mixture, the multimodal $x$) moves. $\varphi$ measures the gap. VIMP names the subset. A small policy spends $B$. That is the whole loop. It looks modest next to continual learning, and it is the one loop in which the CFPerm object is actually used as designed.

---

## 7. Recomputing the score is not online learning

$\varphi$ can be evaluated on a new batch, a new session, a new window. Nuisances can be refit slowly; $\hat\tau$ can be refit on a buffer; a permutation audit can be run rarely. That is **the same batch procedure, applied again**. It is not an online-learning algorithm.

- Online learning: $\theta_t$ is a martingale of losses; the predictor *is* the thing being updated.
- Here: $\theta$ is constant. What is recomputed is a **diagnostic ranking** of subsets. The predictor you ship is still $p_\theta$, possibly composed with a gate that reads $(S^\star, \hat\tau)$.

If it helps: think of PO-risk localization as influence diagnostics / LOCO for a frozen model against a moving batch label — not as SGD. The “update” is the update of the *map of where the gap is*, so that the next intervention is aimed. The generator does not get a vote.

---

## 8. What we are willing to claim

- Localization of PO-risk is a ranking of subsets by contribution to an orthogonal discrepancy. Under a budget, that ranking is exactly “which part is most worth adapting.”
- The LLM (or any foundation encoder) stays a frozen generator. Adapting $p_\theta$ is a different, unidentified, non-orthogonal problem.
- Unique CS/CD explanation remains impossible. The identified product surface is: **name $S^\star$, intervene only there, call $p_\theta$ again.**

What we are not claiming: online learning; causal user-level effects under confounding; that the model “understands” the shift; that a heatmap of $p_\theta$ is a substitute for VIMP of $\hat\tau$.

The generator is frozen. The statistic points at the subset. That is the idea.
