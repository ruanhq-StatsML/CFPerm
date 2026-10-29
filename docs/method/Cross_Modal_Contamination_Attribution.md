# Multi-Modal Attribution use-case: cross-modal contamination, layer by layer

**Opportunity.** Multimodal systems fail when one modality’s shift **leaks into**
another’s scores or into the joint predictor — not only when a single block’s
marginal $P(X_m)$ moves. The FSDS stack already measures modality shares
(RF-Domain / MMD / PO-risk) and hop actuators (modality-π, PO masks). The
missing product story is **contamination**: attribute *which layer* the leak
lives at, then pick the cheapest actuator. That is the use-case.

Grounded in the MSR-VTT window extract already in-repo
(`docs/method/MSRVTT_Multimodal_Attribution_tables_only.tex`,
`results/po_hop_attribution/`) and the locked hop recipe
(Amazon `hop_ridge`; multimodal modality-π; PO-risk = **mask only**).

---

## 1. What “cross-modal contamination” means here

Write the fused representation / score as a composition of modality blocks
$X=(X_v,X_a,X_t)$ and a learner $f$ (retrieval, captioning head, cross-modal
Ridge, uplift-style meta probe — same geometry).

**Contamination** is any of:

| Pattern | Intuition | Typical false read if you only look at one block |
| --- | --- | --- |
| **Support leak** | $P(X_m)$ shifts and the domain classifier / MMD still “sees” it through correlated $X_{m'}$ | Blame the wrong modality |
| **Concept leak** | $Y\mid X$ or the cross-modal map $X_{m'}\leftarrow g(X_m)$ drifts while marginals look quiet | “Everything stable” until ΔMSE / AUUC blows |
| **Fusion soak** | Late fusion / attention lets a noisy block dominate the joint score | High joint error, flat per-block VIMP |
| **Null structural** | A block is window-invariant by construction (e.g. clip-level text on this extract) | Language “doesn’t matter” — wrong conclusion |

On the MSR-VTT sliding-window board, batch $W$ = early vs late windows
*inside* each video. Empirically the **video** block dominates covariate
shift shares; **text** is near a structural null for $W$-shift (caption
constant across windows). Contamination questions become: when audio moves,
how much of the *text prediction* error (video+audio→text) is audio’s fault
vs video’s; when video dominates RF shares, is that pure $P(X_v)$ or fusion
soak into the joint probe.

---

## 2. Layer stack (attribution hierarchy)

Attribute **outside-in**. Each layer answers one question and gates the next.
Do not jump to patch masks before modality shares; do not reweight features
by PO VIMP for training (already lost to `hop_ridge` on Amazon).

| Layer | Question | FSDS / repo probe | Output |
| --- | --- | --- | --- |
| **L0 · Global hop** | Did anything move? Cov vs concept? | MMD$^2$, domain RF AUC on full $X$; PO-risk / DRPerm on stacked batches | typed drift $(c,\delta)$; go/no-go |
| **L1 · Modality mass** | Which block carries the mass of the shift? | RF-Domain VIMP → $\pi_m$; coord-MMD VIMP; PO-risk VIMP / LOGO | simplex shares $(\pi_v,\pi_a,\pi_t)$ |
| **L2 · Within-block localization** | Where inside the dominant block? | Region / patch bins; per-video Friedman–Wilcoxon on shares; Holm cells on $\|d\|$ | region list / patch set |
| **L3 · Cross-modal coupling** | Does block $m$’s shift **contaminate** prediction that uses $m'$? | Cross-modal Ridge / hop (e.g. $X_v{+}X_a\to X_t$); leave-one-modality-out ΔMSE; PO **mask** on donor patches | contamination score $C_{m\to m'}$ |
| **L4 · Actuator** | Cheapest fix given the layer | modality-π hop weights; NarrowSupport-style overlap; PO mask; EWMA λ from `(ĉ,δ̂,Δπ)`; full relearn last | action ticket |

**Contamination score (operational).** On a hop $t{-}1\to t$, for donor
modality $m$ and target head $m'$ (or joint $Y$):

\[
  C_{m\to m'}
  =
  \underbrace{\Delta\mathrm{MSE}^{(-m)}}_{\text{LOGO / drop donor}}
  \;-\;
  \underbrace{\Delta\mathrm{MSE}^{(\mathrm{mask}\,m)}}_{\text{PO-mask donor patches}}
\]

or, in share form, the change in target error when donor is
π-downweighted vs zero-masked. Sign convention: large positive LOGO with
small mask Δ ⇒ **fusion soak** (need fusion/π, not pixel surgery); large
mask Δ with stable π ⇒ **localized donor concept** (mask / refresh that
block); both large ⇒ **compound contamination** (π + mask + consider
relearn).

---

## 2b. How contamination arises · how to characterize it · intermediate CF

This is the methodological core: **contamination is a counterfactual claim
about an intermediate representation**, not a gradient highlight on pixels.

### How contamination arises (mechanisms)

Write a hop $t{-}1\to t$ with modality encoders $h_m=\phi_m(X_m)$ and a
fusion / head $f(h_v,h_a,h_t)$ (or cross-modal Ridge $h_t \approx A[h_v;h_a]$
on the MSR-VTT board). Contamination appears when a change that *lives in*
$h_m$ changes the risk of a head that is *supposed* to read $m'$ (or the
joint), through one of four paths:

| Mechanism | Path | Intermediate object that moved |
| --- | --- | --- |
| **Encoder drift** | $P(X_m)$ or $\phi_m$ shifts ⇒ $h_m$ moves | block embedding / patch tokens of $m$ |
| **Correlation proxy** | $h_m$ correlated with $h_{m'}$; domain RF attributes mass to $m'$ | wrong-block VIMP; true driver is still $h_m$ |
| **Fusion soak** | attention / late concat lets $\|h_m\|$ dominate $f$ | fusion hidden $z=F([h_v;h_a;h_t])$ |
| **Head map drift** | $f(\cdot\mid h_{m'})$ changes while $h_{m'}$ looks stable | head weights / PO surface on $z$ |

Input-only LOGO (drop raw $X_m$) **confounds** these: dropping $X_v$ also
kills every downstream $h_v$ *and* anything that was proxying through $h_v$.
That is why the use-case insists on **intermediate counterfactuals**.

### How to characterize contamination

Three complementary characterizations (use together):

1. **Mass (associational).** L1 shares $\pi_m$ from RF / MMD / PO VIMP —
   “where does the hop *look* active?” Necessary, not sufficient.
2. **Risk delta (interventional).** Change one intermediate, hold the rest,
   remeasure target risk — “what *would* the head do without that path?”
3. **Edge signature.** Compare LOGO-at-input vs CF-at-intermediate vs
   PO-mask-at-coords; the pattern diagnoses mechanism (table below).

Operational contamination on an edge $m\to m'$:

\[
  C_{m\to m'}
  \;=\;
  R\!\bigl(f(h^{\mathrm{obs}})\bigr)
  \;-\;
  R\!\bigl(f(h^{\mathrm{cf}}_{m\leftarrow \mathrm{ref}})\bigr)
\]

with $h^{\mathrm{cf}}_{m\leftarrow\mathrm{ref}}$ the counterfactual intermediate
(defined next). Positive $C$ ⇒ donor $m$’s observed intermediate *hurts*
the target relative to the reference world.

### Intermediate counterfactual attribution (the right do-operator)

**Claim.** Attribute contamination by intervening on $h_m$ (or on fusion
$z$), not only on raw coordinates of $X_m$.

Fix a hop and a target risk $R$ (MSE of cross-modal Ridge, retrieval loss,
AUUC of a monitor head — same geometry). Let $h^{\mathrm{obs}}$ be the
observed intermediates on the evaluation slice. Choose a **reference
world** for donor $m$ (pick one; report which):

| CF world | do-operator on intermediate | Intuition | Repo analogue |
| --- | --- | --- | --- |
| **Past block mean** | $h_m \leftarrow \mu_{t-1}^m$ (keep other $h_{m'}$) | “as if modality $m$ had not hopped” | `block_means` + modality-π cosines |
| **Zero / noise / missing** | $h_m \leftarrow 0$ or noise on top-PO coords *inside* $h_m$ | “kill localized donor concept” | `apply_mask` / `run_msrvtt_patch_mask` |
| **LOGO block** | drop $h_m$ entirely from $f$ | “no path from $m$” | LOGO / leave-one-group-out |
| **Fusion splice** | $z\leftarrow F([h_v^{\mathrm{obs}}; h_a^{\mathrm{cf}}; h_t^{\mathrm{obs}}])$ | intervene *after* encoders | fusion-adapter EWMA target |

Then the **intermediate CF attribution** for donor $m$ on target $m'$ is

\[
  \mathrm{Attr}_m
  =
  R(f(h^{\mathrm{obs}})) - R\!\bigl(f(\mathrm{do}(h_m=h_m^{\mathrm{ref}}))\bigr).
\]

Summing $\mathrm{Attr}_m$ over donors does **not** have to recover total
risk (interactions / non-additive fusion). Report the vector
$(\mathrm{Attr}_v,\mathrm{Attr}_a,\mathrm{Attr}_t)$ and the pairwise
$C_{m\to m'}$, not a forced simplex unless you L1-normalize for display.

**Why this is better than input Grad-CAM / raw LOGO alone**

- Grad-CAM is local linearization inside one forward; hop contamination is
  a **batch-counterfactual** across $W$.
- Raw-input LOGO mixes encoder drift with fusion soak; intervening on
  $h_m$ *after* $\phi_m$ isolates fusion/head paths.
- Past-mean CF (`$h_m\leftarrow\mu_{t-1}^m$`) matches the locked hop
  geometry already used for IW: the causal reference is always the last
  observed batch, never the future.

### Mechanism signature (read the CF table)

| Input LOGO $\Delta R^{(-X_m)}$ | Intermediate CF $\mathrm{Attr}_m$ (past-mean) | PO-mask $\Delta R$ inside $h_m$ | Diagnosis |
| --- | --- | --- | --- |
| large | large | small | **fusion soak** — whole block pollutes $f$; use π / fuse-gate |
| large | large | large | **localized concept in donor** — PO-mask / refresh $\phi_m$ |
| large | small | small | **correlation proxy** — mass on $m$ is spurious; check $m'$ |
| small | large | * | **head/fusion path** without raw-coord mass — audit $F$ |
| ~0 on text | ~0 | ~0 | **structural null** — report, do not conclude “unused” |

MSR-VTT anchors: video dominates L1 mass; hop+π_m moves MSE
$0.0898\to 0.0896$; video patch mask $\Delta$MSE $+0.0009$ vs audio
$-0.0004$. That signature is closer to **mild fusion/π** than to a raging
localized audio concept — exactly the kind of intermediate-CF read the
board is for. Harder hops (injected audio concept, varying text) are
where $\mathrm{Attr}_a$ and mask Δ separate.

### Minimal recipe to compute on a hop (ops)

1. Encode $h_m^{\mathrm{obs}}=\phi_m(X_m)$ on eval batch $t$ (or use the
   stored block features as $h_m$ when $\phi_m$ is frozen).
2. Build $h_m^{\mathrm{ref}}=\mu_{t-1}^m$ from `block_means`.
3. Score $R^{\mathrm{obs}}=R(f(h^{\mathrm{obs}}))$ and
   $R_m^{\mathrm{cf}}=R(f(\mathrm{do}(h_m=\mu_{t-1}^m)))$.
4. $\mathrm{Attr}_m = R^{\mathrm{obs}}-R_m^{\mathrm{cf}}$; also run PO-mask
   inside $h_m$ for the mechanism table.
5. Emit L1 $\pi_m$ × intermediate $\mathrm{Attr}_m$ cross-check → actuator
   (π-hop / PoMask / EwmaFusion / RelearnHead).

**Product line.** Intermediate CF attribution is how you say
“audio contaminated the text head on this hop” in a way that is
*falsifiable* (swap the do-world, rerun) and *actionable* (mask vs π vs
relearn), without inventing a new multimodal loss.

---

## 3. Opportunity (where the product value is)

### O1 · Stop blaming the wrong modality

Single-share dashboards (only RF VIMP) over-credit the tallest block.
MSR-VTT: video≈0.77 / audio≈0.19 / text≈0.04 under RF-Domain, but text’s
near-zero is **structural** on this extract — not “language is useless.”
Layer L1+L3 separates “primary for $W$-shift” from “primary for the
cross-modal head.” Opportunity: **one contamination ticket per hop** that
names donor → target, not a single pie chart.

### O2 · Cheapest actuator before retrain

Locked recipe already says: Amazon → heatmap hop-ridge; multimodal →
modality-π; PO-risk → **mask only** (do not amplify/downweight features by
PO VIMP — both lost to `hop_ridge`). Contamination layers map to that
recipe:

| Diagnosed layer | First actuator (cheap) | Avoid |
| --- | --- | --- |
| L1 mass on $m$, L3 quiet | modality-π hop $w_s=\exp(\gamma\sum_m\pi_m(\cos\mu_s^m,\mu_{t-1}^m-1))$ | full backbone retrain |
| L3 mask ΔMSE large on $m$ | PO zero/noise/missing mask on top-PO patches | feature reweight by VIMP |
| L0 concept, L1 flat | PO-risk / DRPerm → refresh head / EWMA λ↑ | silent hold |
| L0 cov, support poor | refresh REF / narrow overlap support | LOCO cuts under mix shift |

Live MSR-VTT cross-modal hop: $0.0898\to 0.0896$ with π_m (video≈0.71,
audio≈0.26) — small but **signed** improvement without a new loss.
Amazon top-16 PO mask: +0.13 ΔMSE — clear donor localization.
Opportunity: **contamination routing = pick π vs mask vs λ**, not a new
objective.

### O3 · Compound vs pure shift (same as uplift monitor geometry)

Uplift monitoring taught: AUUC localizes ranking; IV tiers localize
identification; business gates when they disagree. Multimodal analogue:

| L1 (modality mass) | L3 (coupling) | Read | Act |
| --- | --- | --- | --- |
| $m$ dominates | $C_{m\to m'}$ small | Pure donor shift | π on $m$; hold $m'$ |
| $m$ dominates | $C_{m\to m'}$ large | Contaminating donor | π **and** mask / fuse-gate |
| flat shares | $C$ large for one edge | Fusion soak | attention/fusion audit; NarrowSupport |
| text structural null | any | Do not conclude “language unused” | report structural zero |

Opportunity: reuse the **cross-check matrix** pattern from the uplift
actionables — same product language, new physics (modalities instead of
LIVE quintiles).

### O4 · Online EWMA from attribution, not from vibes

`(ĉ, δ̂, Δπ) → attribution_ewma_rate → λ_new`:
covariate hop → small λ (keep past ensemble); concept / π jump → large λ
(refresh mix). Contamination adds $\Delta C_{m\to m'}$ as an optional
fourth signal: a spike in coupling without π change ⇒ fusion soak ⇒ raise
λ on the **fusion adapter**, not on every backbone.

### O5 · Eval / paper surface that is still scarce

Most multimodal “attribution” stops at Grad-CAM / attention rollout inside
one forward pass. FSDS contamination is **batch-hop causal geometry**:
domain RF + MMD + PO-risk + LOGO/mask across $W$. Opportunity for the
writeup: a tables-only use-case (like MSR-VTT attribution) that answers
“which layer leaked?” with shares, $C_{m\to m'}$, and an actuator column —
auditable for hourly windows.

---

## 4. Empirical anchors (already measured)

### Modality mass on $P(X\mid W)$ (MSR-VTT windows)

| Modality | RF-Domain | MMD | PO-risk | Role in $W$-shift |
| --- | --- | --- | --- | --- |
| Video (768-d) | 0.772 | 0.897 | 0.813 | primary |
| Audio (512-d) | 0.191 | 0.101 | 0.174 | secondary |
| Text (768-d) | 0.037 | 0.002 | 0.013 | structural null† |

Friedman / Wilcoxon reject equal shares. †Clip-level caption on this
extract ⇒ window-invariant; not a claim that language is uninformative
for retrieval.

### Cross-modal hop + mask (video+audio → text)

| Probe | Number |
| --- | --- |
| hop / hop+π_m MSE | 0.0898 / **0.0896** |
| π shares (live) | video≈0.71, audio≈0.26 |
| mask video patches (zero) ΔMSE | +0.0009 |
| mask audio patches (zero) ΔMSE | −0.0004 |
| Amazon top-PO mask zero/noise | +0.126 / +0.144 |

Reading: π_m slightly helps the cross-modal Ridge; video mask hurts more
than audio mask on this hop — consistent with video-primary L1 mass, mild
L3 coupling. Contaminations that would be “large” are still the opportunity
surface when hops get harder (injected audio concept, text that actually
varies with $W$, fusion models with attention soak).

---

## 5. Worked contamination playbook (ops)

1. **L0.** Run MMD + domain AUC + PO-risk on the hop. If all quiet → stop.
2. **L1.** Compute $\pi_m$ via `benchmark_feature_selection` →
   `modality_mass` (default RF-Domain; swap hook for MMD / PO).
3. **L2.** If one $m$ dominates, localize regions (Holm $\|d\|$, per-video
   shares). If shares flip across videos (e.g. video 13 audio-dominant on
   the board), emit **per-video** tickets — global π alone is wrong.
4. **L3.** Fit / score the cross-modal head; LOGO drop $m$ and PO-mask $m$;
   fill $C_{m\to m'}$. Cross-check against L1 (matrix in §3 O3).
5. **L4.** Emit one actuator:
   - π-only / mask-only / π+mask / EWMA λ / RefreshRef / Relearn
   - never “scale features by PO VIMP” for training
6. **Next hop.** Re-measure $C$; if contamination falls under π/mask, keep;
   if not, escalate Relearn on the contaminated head only.

---

## 6. Actionable items (multimodal contamination v1)

| Pri | Code | Item |
| --- | --- | --- |
| P0 | GateHop | No L1–L4 tickets if L0 quiet (MMD flat, domain AUC≈0.5, PO accept) |
| P1 | ModalityMass | Publish $\pi_m$ (RF / MMD / PO) every hop; flag structural nulls |
| P1 | ContamEdge | Compute $C_{m\to m'}$ for each donor→target edge on the deployed head |
| P1 | IntermCF | Intermediate CF: $\mathrm{Attr}_m=R(h^{\mathrm{obs}})-R(\mathrm{do}(h_m=\mu_{t-1}^m))$; also PO-mask inside $h_m$ for mechanism signature |
| P1 | CrossCheckMM | Require L1×L3×IntermCF cell on every reweight / mask ticket |
| P2 | PiHop | Apply modality-π hop weights when L1 mass moves and L3 mild |
| P2 | PoMask | Mask top-PO patches when L3 mask ΔMSE large; do not VIMP-reweight train |
| P2 | PerVideo | If Friedman / per-video dominant flips, ticket per video — not global π |
| P3 | EwmaFusion | Feed $(ĉ,δ̂,Δπ,ΔC)$ into `attribution_ewma_rate`; spike $ΔC$ alone → fusion adapter λ |
| P3 | NarrowSupportMM | Under mix shift, score only overlap support across modalities |
| P4 | RelearnHead | After π+mask trial fails under PO reject — relearn contaminated head only |
| P4 | StructNull | Report structural zeros explicitly (window-invariant text) — never “language unused” |

---

## 7. Why this is the opportunity *now*

- **Stack exists:** modality VIMP, LOGO, PO masks, modality-π hop, EWMA rate —
  already coded (`msrvtt_multimodal_attribution`, `po_hop_attribution`,
  `attribution_adapter`).
- **Gap is product geometry:** name contamination edges and layer→actuator
  routing the same way uplift monitoring named AUUC×IV tiers.
- **Cost asymmetry:** π and masks are hop-cheap; multimodal backbone retrain
  is not. Layer-wise attribution is how you stay on the cheap side.
- **Eval story:** tables-only contamination board (shares + $C_{m\to m'}$ +
  actuator) is publishable next to the existing MSR-VTT attribution tables.

---

## Run

```bash
PYTHONPATH=Python/src python3 scripts/run_msrvtt_multimodal_attribution.py
PYTHONPATH=Python/src python3 scripts/run_po_hop_attribution.py
PYTHONPATH=Python/src python3 scripts/run_attribution_adapter.py
```

Code: `Python/src/msrvtt_multimodal_attribution.py`,
`Python/src/po_hop_attribution.py`,
`Python/src/attribution_adapter.py`  
Companion LaTeX: `docs/method/Cross_Modal_Contamination_Attribution.tex`
