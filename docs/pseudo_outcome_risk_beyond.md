# Pseudo-Outcome Risk Beyond Batch Testing

**What else the score can do — especially online update, and statistically grounded LLM inference.**

This is a writeup, not a prototype. The implementation stays out of this pass on purpose: the object is already in CFPerm; what is missing is a map of *uses* in which a statistical method is the inner loop, not a dashboard after the model has already decided.

The claim in one sentence: **once you have a Neyman-orthogonal, per-example discrepancy, you can update it online and you can hand it to an LLM as a control variate.** Almost everything people currently ask prompting, RAG, RLHF, or “the model’s own confidence” to do is a *policy*. Pseudo-outcome (PO) risk is a *score*. Policies should be functions of scores. That inversion is where statistics plays a core role rather than a supporting one.

---

## 0. Doctrine that does not change

CFPerm already says this, and none of the uses below walk it back:

- Under the alternative, ignorability fails. $\hat\tau$ is **not** a causal CATE. The DR / R-learner residual is a **distance / discrepancy** between two batches.
- Unique decomposition of that discrepancy into covariate shift vs. concept drift is **impossible** (README: three levels — performance drop, $P(X)$ vs $P(Y\mid X)$, and VIMP change each have infinitely many realizations). The identifiable target is **localization**, not disentanglement.
- The permutation / CF-split procedure is the type-I device. Online monitors and LLM gates do **not** inherit that guarantee automatically.

What changes is only the *use* of the residual: from a one-shot test statistic to a streaming control signal.

---

## 1. The object, in one place

Batch indicator $W$ is “old vs new”, “user A vs user B”, “session $t-1$ vs $t$”, “trusted corpus vs current prompt”, “no-tool vs with-tool” — **not** a treatment unless you actually randomized it.

Cross-fit nuisances $\hat\mu(x)\approx\mathbb{E}[Y\mid X=x]$ and $\hat e(x)\approx\mathbb{P}(W=1\mid X=x)$. Two siblings already in the repo:

- **DR / AIPW-style pseudo-outcome** (efficient residual for a contrast):
  \[
  \varphi^{\mathrm{DR}}
  =
  (\hat\mu_1-\hat\mu_0)
  +
  \frac{(Y-\hat\mu_W)\,(W-\hat e)}{\hat e(1-\hat e)}.
  \]
- **R-learner / product residual** (what `DRPerm` / `cfperm_vimp` store as `residual_y * residual_t`):
  \[
  \varphi^{\mathrm{R}}
  =
  (Y-\hat\mu(X))\,(W-\hat e(X)).
  \]

Then $\hat\tau(x)$ is a regression of $\varphi$ on $X$, and

\[
\text{PO-risk}
\;=\;
\mathbb{E}[\hat\tau(X)^2]
\quad\text{or}\quad
\mathbb{E}[\varphi^2].
\]

Four properties that make this the right primitive for everything below:

1. **Orthogonality.** First-order errors in $\hat\mu$ or $\hat e$ drop out of $\mathbb{E}[\varphi]$ when the other nuisance is consistent. Streams and LLMs are exactly the setting where one nuisance will be stale and the other still usable.
2. **One number per example.** After nuisances exist, a new point $(x,w,y)$ produces a scalar $\varphi$ without refitting the permutation test.
3. **Localization is free.** LOCO, permuCATE, RF VIMP on concatenated modalities — anything we already do to $\hat\tau$ — is a readout of the same residual. Layer-1 modality shares, Layer-2 instance selection, Layer-3 token/box localization are *views* of $\varphi$, not a second scientific claim.
4. **It is not an explanation.** $\varphi$ does not say why the user changed their mind. It says *where the two batches remain distinguishable after outcome and propensity have been partialled out*. That is the statistically honest sentence, and it is already enough to act.

The R-loss form makes the same point operationally: $\tau$ is fit by residual-on-residual regression. Online, that is a *fast* regression. The nuisances are a *slow* regression. Two timescales, one orthogonal score.

---

## 2. Online update: the residual as a monitor, not a batch $p$-value

The current CFPerm test is batch-vs-batch: permute $W$, refit, threshold the VIMP quantile. That is the right **validity** device. It is the wrong **runtime** device. You cannot permute 150 times per token, per session, or per video window. You also should not wait until a full new batch has been labeled and then run a test that tells you, after the fact, that the operating point already moved.

Online-update is the statement that **the same $\varphi$ can be computed on arrival**, accumulated, and only occasionally audited by the expensive test.

### 2.1 Split the timescales (this is the architecture)

| Object | Role | Update frequency |
|---|---|---|
| $\hat\mu,\hat e$ | nuisances | slow — every $B$ points, or when a CUSUM on their own residuals fires |
| $\varphi_i$ | per-point score | every arrival (or when $Y$ completes) |
| PO-risk / $\hat\tau$ | discrepancy + localization | fast — SGD, recursive least squares, streaming forest, EWMA of $\varphi^2$ |
| Permutation / CF-split | confirmatory type-I control | rare — end of window, or when the monitor exceeds a bound |

Two-timescale stochastic approximation is the statistical name for this, not “fine-tune more often.” Nuisances move on a slow clock because they are high-capacity and because orthogonality *forgives* them being slightly stale. $\hat\tau$ moves on a fast clock because it is a regression on an already-orthogonalized target.

Double robustness is what makes the split respectable:

- If $\hat e$ is still good and $\hat\mu$ is stale, the product residual still has signal through $(Y-\hat\mu)(W-\hat e)$.
- If $\hat\mu$ is still good and the mixture $W$ is drifting, signal comes through $W-\hat e$.
- If both are badly wrong, the monitor is allowed to go quiet — that is the honest failure mode — and the rare permutation audit is what catches it.

### 2.2 What you actually accumulate

Do not store a new copy of the world. Store a sketch:

\[
\mathcal{S}_t
=
\bigl\{\,
\varphi_i,\;
\hat\tau(x_i),\;
\text{top-$k$ VIMP coordinates},\;
W_i,\;
\text{window id}
\,\bigr\}_{i\le t}.
\]

Operational statistics on $\mathcal{S}_t$:

- **Level.** EWMA of $\varphi_i^2$ or of $\hat\tau(x_i)^2$. This is the current PO-risk.
- **Change.** CUSUM / Shiryaev–Roberts on the same series. This is “the operating point moved,” not “there is shift in the abstract.”
- **Where.** Rolling VIMP / LOCO shares on a reservoir of high-$|\varphi|$ rows. This is Layer-1 localization on a sliding window.
- **Which rows.** Keep the top residual instances (Layer 2). Raw video/audio is optional and expensive; the residual is the memory.

A sufficient online memory for a multimodal product is this sketch. That is how a statistical method stays in the loop without becoming another foundation-model training run.

### 2.3 Delayed labels, missing $Y$, partial information

LLM and preference logs are not i.i.d. batches with $Y$ on arrival.

- $W$ and $X$ are often immediate: this session vs last week; this user cluster vs the calibration pool; this prompt vs the trusted corpus.
- $Y$ may be delayed (thumbs, retention, task success, human rating) or missing.
- Then $\hat e$ can be updated on $(X,W)$ continuously; $\varphi$ is *completed* when $Y$ shows up. Until then, $\hat\tau(X)$ from the last fitted residual regression is still a **prediction of discrepancy**, which is enough to gate.

This is ordinary two-phase sampling / missing-at-random completion of an AIPW score. It is not a new causal story, and it does not require pretending that delayed thumbs are a CATE.

A useful split in products:

1. **Immediate gate** — $\hat\tau(x_{\text{now}})$ and current EWMA, no $Y$ yet.
2. **Completed update** — when $Y$ arrives, write $\varphi$, take a SGD step on $\hat\tau$, maybe a slower step on $\hat\mu$.
3. **Audit** — when the CUSUM fires, or every $K$ windows, run permute-then-refit on a buffer.

The immediate gate is what inference needs. The completed update is what learning needs. The audit is what validity needs. Mixing those three clocks is how people either (a) never ship or (b) ship a $p$-value with no type-I meaning.

### 2.4 Rolling localization on concatenated modalities

Given an MSR-VTT-style matrix $(N\cdot T,\, d)=\mathrm{video}\oplus\mathrm{audio}\oplus\mathrm{text}$, the online object is a sliding window of rows. At each window:

- PO-risk: is this window distinguishable from the reference batch?
- VIMP shares: is the discrepancy sitting in video, audio, or text?
- Keep only the salient rows (Layer 2) and, if needed, tokens/boxes (Layer 3).

Because unique CS/CD decomposition is impossible, the online system should publish **shares and localizers**, not “concept drift = 0.6.” The monitor is allowed to say: *this session’s shift is audio-heavy relative to the calibration pool*. That sentence is actionable (re-encode audio, fetch a better transcript, spend CLAP compute). The sentence “the concept has drifted” is not.

### 2.5 Sequential testing vs monitoring (do not confuse them)

Online, two different questions get collapsed:

- **Monitoring.** Has the operating point moved enough that the *policy* should change? EWMA / CUSUM on $\varphi^2$. No type-I claim. This is the inner loop.
- **Testing.** Can we reject exchangeability of $W$ at level $\alpha$ on this buffer? Permutation / CF-split. Type-I claim. This is the audit.

You can sequentialize the test (alpha-spending, always-valid $p$-values, e-values / betting scores on $\varphi$). That is a real statistical project. It is **not** required to start using the monitor. The mistake is to treat the EWMA as a $p$-value, or to refuse to gate until a batch test has rejected.

Always-valid sequential tests on the orthogonal score are a natural next paper. The product does not have to wait for that paper to use the score as a gate.

### 2.6 What “online-update” is not

It is not “fine-tune the LLM every 100 users.”
It is not sequential Bayesian causal identification.
It is not a replacement for the permutation test.

It is **orthogonal score + slow nuisances + fast residual regression**, with the permutation test as an occasional audit. Validity lives in the audit; power and latency live in the monitor. That split is the whole online story.

---

## 3. Empowering LLM inference: the residual as a control variate

This is larger than a shift detector. An LLM at inference time is a sequence of decisions under a user/context mixture that is not the pretraining mixture. The industry currently stuffs that problem into prompting, RAG, or RLHF. Those are **function classes**. They are not scores.

PO-risk is a score. The LLM can stay frozen. The statistician owns the gate.

That is the inversion: **do not ask the 70B model to know that it is OOD; compute $\varphi$ and change the policy.** The LLM remains a generator. The residual is the controller. This is where a statistical method occupies the inner loop of inference, which is precisely where it is usually absent.

### 3.1 The interface

```
reference batch D0          current context Dt
     |                           |
     +------ nuisances μ, e -----+
                     |
              pseudo-outcome φ
                     |
         τ(x), PO-risk, VIMP localization
                     |
      ┌--------------┼--------------┐
      ▼              ▼              ▼
   route/abstain   reweight       localize
   (which model,   (decode,       (which
    RAG, tools)     retrieve)      span / modality)
```

Three control channels, all statistically named:

1. **Scalar gate.** If PO-risk of $(x_{\text{prompt}}, W=\text{now})$ vs $D_0$ exceeds a threshold (or a conformal quantile of $\varphi^2$ on a calibration stream), change the *policy*: bigger model, retrieve, ask a clarifying question, refuse, or spend more test-time compute. The LLM does not have to “know” it is OOD; the residual says so.
2. **Coordinate gate.** VIMP / LOCO on $\hat\tau$ says *which block of $x$ carries the discrepancy* — image tokens vs transcript vs user metadata. That is Layer-1 modality attribution used as an encoder / attention routing prior: spend ViViT compute on the video block if video share dominates; don’t.
3. **Span gate.** Layer-3 localization (token / box) on the same residual is a **soft mask** over the prompt or over candidate continuations: downweight spans that look like the shift drivers, or conversely *highlight* them for the user. This is not integrated gradients of the LLM. It is a meta-learner residual mapped back onto the sequence. The two heatmaps will disagree; that disagreement is the point.

The LLM never has to “explain the drift.” The impossibility result already says that explanation is not identifiable. Localization of $\varphi$ is.

### 3.2 Concrete inference jobs where a statistic is the core, not a helper

These are not metaphors. Each one is: define $W$, compute $\varphi$, let a small policy read $(\hat\tau(x), \text{shares})$.

**Mixture / preference shift (the CFPerm $W$).**
Users are batches. A session is a new $W=1$. Online PO-risk vs a calibration cohort is a *user-mixture detector*. The LLM’s next action (style, verbosity, tool use, safety template, system prompt) should be a function of $(\hat\tau(x), \text{modality shares})$, not of a vibes-based persona vector.

DPO / RLHF already encode preferences as implicit rewards. They answer “what does this user like *on average in the training mixture*.” PO-risk answers a different, statistically sharper question: **has this user’s mixture moved relative to the policy’s training mixture, and along which coordinates?** If it has, you do not need to re-RLHF. You gate, retrieve, or update a residual head.

**Retrieval vs parametric knowledge.**
Treat retrieved passages and the model’s ungrounded continuation as two batches. High PO-risk localized on the continuation is a hallucination localizer. High PO-risk localized on retrieval is a corpus-mismatch localizer (wrong index, stale docs, query drift). The *fix* is different — generate less vs retrieve better — but the *score* is the same object. Today those two failures are both called “hallucination.” The residual distinguishes them by localization.

**Draft / verify and test-time compute.**
Speculative decoding and “think longer” are currently heuristics (always verify; verify when the draft is long; verify when entropy is high). A statistically honest rule is: spend extra verify steps when $\varphi^2$ of the draft continuation against the user’s current mixture is large. Cheap draft, expensive verify *conditional on the residual*. Entropy of the LLM is not a discrepancy against $D_0$. $\varphi$ is.

**Tool routing and agents.**
An agent that always searches is a waste; one that never searches hallucinates. Route to a tool when PO-risk says the current $x$ is distinguishable from the no-tool calibration set. The tool policy can be a small $\hat\tau$-driven classifier — contextual bandits with an orthogonal score as the context — not another LLM-as-router. The LLM-as-router has no orthogonality, no delayed-$Y$ completion, and no audit.

**Calibration of “confidence.”**
LLM softmax is not a probability of being exchangeable with a trusted batch. $\varphi$ (or a cross-fit DR score of correctness) *is* a residual with an asymptotic expansion. You can conformalize $\varphi^2$ on a stream and get a coverage statement about “this generation is exchangeable with the trusted batch,” which is a weaker and truer claim than “the model is 0.92 sure.” Abstention, refusal, and “ask a clarifying question” should hang off that coverage statement, not off a verbalized confidence token.

**Continual personalization without touching weights.**
Keep $\hat\mu,\hat e$ global. Let $\hat\tau$ be per-user or per-cohort, updated online. The LLM is a frozen generator; personalization is a residual head. That is the statistical analogue of LoRA, except the object being updated has an influence function, can be sketched, can be audited by permutation, and does not require storing 70B gradients. Multi-turn conversation is a rolling $W$: last $k$ turns vs the user’s calibration history. The residual head is what should move turn-by-turn, not the transformer.

**Evaluation and A/B without lying about causality.**
When you A/B two prompts or two models, $W$ is the arm. If assignment is randomized, $\varphi$ recovers a causal contrast and you may say so. If assignment is confounded (power users get the new model), CFPerm’s own warning applies: $\varphi$ is still a discrepancy, not an effect. Publishing both sentences is the method. The industry currently publishes only the causal sentence. The statistical contribution is to keep the other one in the room.

**Mixture-of-experts / encoder routing.**
On concatenated multimodal $x$, Layer-1 shares of $\hat\tau$ are a routing prior over experts (video encoder vs audio encoder vs text). This is cheaper than running every encoder always, and it is not “the LLM decided the video mattered.” It is VIMP of an orthogonal residual. For streaming video, the same shares can turn encoders on and off window-by-window.

**Safety and policy templates.**
A safety classifier that always fires is unusable; one that never fires is a scandal. Treat “trusted safe batch” vs “current prompt” as $W$. Gate the strict template when PO-risk is high and localized on the spans the safety policy actually cares about. Again: the LLM does not classify the shift; $\varphi$ does, and the template is a policy.

### 3.3 Why the LLM cannot replace the score

People imagine the model will introspect the shift (“I notice this user is different; I should retrieve”). It will not, identifiably.

The impossibility result already says: a given performance drop, a given distributional distance, a given VIMP change, each has infinitely many CS / CD / dependence realizations. Asking GPT to “explain the drift” is asking it to pick a point in that equivalence class. It will pick a fluent one. Fluency is not identification.

The move is the opposite:

- **Do not explain, localize.**
- **Do not backprop the 70B, update the residual.**
- **Do not permute every token, monitor $\varphi$ and audit occasionally.**
- **Do not ask the LLM whether it is OOD; compute a discrepancy against $D_0$.**

That is a statistical method occupying the inner loop of inference. Prompting, RAG, and RLHF remain available as *actions the gate may choose*. They stop being the detector.

### 3.4 What the residual is allowed to say to the LLM

A small, honest API — the whole interface:

```
gate(x_now) -> {
  po_risk:          scalar,          # EWMA / τ(x)²
  alarm:            {none, watch, act},
  shares:           {video, audio, text, ...},  # Layer 1
  local_rows:       ids,             # Layer 2
  local_spans:      token / box ids, # Layer 3
  exchangeable:     conformal interval for φ² vs D0,
  causal:           false unless W was randomized
}
```

The LLM (or a 10k-parameter policy in front of it) reads this structure. It does not get a paragraph of chain-of-thought about “root cause.” Shares, localizers, a conformal bit, and a flag that this is not a CATE. That is enough to route, retrieve, abstain, or spend compute. It is also all that is identified.

---

## 4. Other things the same score can still do

Once $\varphi$ is in the inner loop, a list that is easy to underestimate:

- **Federated / on-device $\hat\tau$.** Nuisances stay on the server; each device updates a local residual head. Communication is sketches of $\varphi$, not raw text.
- **Replay and debugging.** A production incident is a buffer of high-$|\varphi|$ rows with VIMP shares. That is a statistically named incident report. “The model got worse” is not.
- **Data valuation / collection.** Collect more labels where $|\varphi|$ is large and localization is unstable. That is orthogonal-score active learning, not uncertainty sampling of the LLM.
- **Teacher–student / distillation.** Distill only on rows where PO-risk against the teacher batch is small; treat high-PO-risk rows as OOD and do not pretend the student should match them.
- **Pretrain vs posttrain vs deploy mixtures.** Three batches, three pairwise $\varphi$’s. Localization tells you whether a deploy failure sits in the SFT mixture, the preference mixture, or the live user mixture. Today those are collapsed into “alignment.”
- **Multilingual / regional routing.** $W$ = locale. Same machinery; shares will sit on script, ASR, or retrieval index rather than on “the model is biased,” which is not a statistical sentence.
- **Always-valid sequential tests** (e-values on $\varphi$) if a paper-facing type-I guarantee is required in streaming. Optional; the monitor does not depend on it.

The pattern is the same every time: **name $W$, form $\varphi$, localize $\hat\tau$, act with a small policy, audit rarely.** The imagination is not in inventing new neural modules. It is in noticing how many product decisions are unnamed discrepancies.

---

## 5. Named experiments (for later; not this pass)

These are names of experiments, not an implementation plan. The prototype is intentionally left to a from-scratch pass.

1. **Streaming PO-risk on a concatenated multimodal matrix.** EWMA of $\varphi^2$; alarm vs permutation audit every $K$ windows. Report modality shares over time. Do not report “% concept drift.”
2. **Delayed-$Y$ completion.** Update $\hat e$ every prompt; complete $\varphi$ when feedback arrives; compare to an oracle that had $Y$ immediately. The gap is the cost of two-phase sampling, which is a statistical quantity.
3. **Gate a frozen LLM.** Two policies: always-RAG vs RAG-when-PO-risk-high. Metrics: tokens, latency, and a proper scoring rule on a held-out mixture shift. The residual should win on the scoring rule *and* on tokens if the gate is any good.
4. **Span mask from Layer-3.** Use token-level localization of $\hat\tau$ as a decode constraint; compare to attribution heatmaps of the LLM itself. They will disagree. Report both; do not average them.
5. **User-mixture routing.** $W$ = cohort. Online $\hat\tau$ selects system prompt / tool set. Keep the LLM frozen. Ablate against an LLM-as-router with the same action set.
6. **Conformal $\varphi^2$.** Coverage of “exchangeable with $D_0$” on a stream with mixture shifts. Compare to verbalized confidence and to softmax entropy. The statistical method should be calibrated; the LLM internals should not be.

Validity checks that should travel with every prototype: cross-fitting of nuisances, clipped propensity, a permutation or CF-split audit on a slower clock, and an explicit sentence that $\varphi$ is not a CATE unless $W$ is randomized.

---

## 6. What we are willing to claim

- PO-risk is an orthogonal, per-example discrepancy that can be updated online and read out as localization.
- That residual is a legitimate control variate for frozen-LLM inference: routing, retrieval, test-time compute, span/modality gating, and residual-head personalization.
- Unique source decomposition remains impossible. The product surface is **localization + a monitor + an occasional test**.
- This is a place where a statistical method is the core: the LLM is a generator; $\varphi$ is the controller.

What we are not claiming: that meta-learners give causal user-level effects under confounding; that an LLM with a residual head “understands” the shift; that online $\varphi$ replaces a properly permutated type-I guarantee; that the model can identifiably explain drift.

The statistical method is the inner loop. The LLM is the generator. That split is the whole idea.
