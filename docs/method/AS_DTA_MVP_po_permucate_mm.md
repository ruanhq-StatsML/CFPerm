# MVP for professor — Mag×Var PO, PermuCATE relation, multimodal attribution

## 0. One-screen pitch

> **MVP:** When mean OOD/PO rises, merge **Magnitude × Variance** of PO across
> bootstrap sub-environments to resolve **benign vs malignant shift** (2×2).  
> **Honest estimator note:** the AIGC Mag×Var MVP currently uses a **PO-risk
> proxy** $\mathrm{PO}_i=1-\hat p(Y_i\mid x_i)$, *not* a full PermuCATE fit.  
> **Extension:** same Mag×Var geometry plugs into your **PermuCATE** DR
> pseudo-outcome / R-risk, and into **multimodal attribution** (modality-block
> PO-risk VIMP) without changing the routing story.

## 1. How PO-risk is estimated *in this MVP* (important)

### What we use now (AS-DTA / Mag×Var boards)

| step | estimator |
| --- | --- |
| Detector | TF-IDF + logistic $\hat p(Y=1\mid x)$ on REF |
| Unit PO **proxy** | $\mathrm{PO}_i = 1-\hat p(Y_i\mid x_i)$ (or squared error twin) |
| Magnitude | $\mu=\mathbb{E}[\mathrm{PO}_i]$ |
| Variance | $s=\mathrm{Std}(r_b)$, $r_b=$ mean PO on bootstrap window $b$ |

**Not used in the AIGC Mag×Var path today:** cross-fit nuisance $(\hat\mu_t,\hat\pi)$,
DR pseudo-outcome $\varphi(Z)$, or PermuCATE feature-conditional VIMP.

Why proxy first: text LIVE + binary $Y$ already gives a labeled risk; Mag×Var
only needs a unit risk that rises under shift. Proxy is CPU-cheap for the MVP
board loop (HC3 / news20).

### What *your* PermuCATE estimates (canonical PO-risk)

From `Python/src/permuCATE_vimp.py`:

$$
\varphi(Z)=\frac{(Y-\mu_T(X))(T-\pi(X))}{\pi(X)(1-\pi(X))}+(\mu_1(X)-\mu_0(X)),
\qquad
R=\mathbb{E}\big[(\varphi(Z)-\hat\tau(X))^2\big].
$$

VIMP = permute / conditional fit importance on that R-risk.

FSDS / uplift monitors in-repo (`uplift_po_monitor`) use the same **family**:
batch $W$, reward $Y$, cross-fit $m,\hat e$, DR-style `po`, then
$R_{\mathrm{PO}}=\mathbb{E}[\hat\tau^2]$ — polarity aligned with CFPerm / PermuCATE,
still a compact RF stand-in rather than the full `permuCATE_fit_nuisance` stack.

### Bridge (justify the extension)

| layer | role |
| --- | --- |
| Mag×Var MVP | **shape router** on any unit risk field |
| PO proxy | fast when $Y$ is the supervised label (AIGC detect) |
| PermuCATE $\varphi$ / $R$ | when the scientific object is **treatment effect** ($T$=batch/hop, $Y$=reward) |

Swap unit field $\mathrm{PO}_i^{\mathrm{proxy}}\;\mapsto\;(\varphi_i-\hat\tau_i)^2$
(or $|\varphi_i|$) and keep the **same** Magnitude × Variance + 2×2 eval.
The MVP justifies the *geometry*; PermuCATE upgrades the *causal risk*.

**Professor line:** Mag×Var is estimator-agnostic. Current demo uses a supervised
PO proxy; production causal hops should plug PermuCATE / DR PO without changing
the benign↔malignant resolution logic.

## 2. Multimodal attribution — same intuition, justify

MSR-VTT path (`msrvtt_multimodal_attribution.py`):

1. Hop as treatment $W$ (batch / time), multimodal $X=(V,A,T)$.  
2. Build PO (DR / meta-learner family).  
3. **Attribution** = which modality block drives PO-risk:
   - modality mass of VIMP  
   - leave-one-group-out / **group-permute PO-risk** (shuffle one modality, keep $\tau$ fixed)

**Why intuitive with Mag×Var**

- Mag×Var asks: is LIVE risk *heterogeneous across sub-envs*? → benign vs malignant **shift type**.  
- Multimodal attribution asks: is that risk *heterogeneous across modality blocks*? → which **channel** to freeze / reweight.

Same “variance of PO mass” story: sparse malignant mass in text vs uniform video cov shift.
Landing: route actuator **by modality** (freeze text lexical; hold audio; conformal on ASR domain) once Mag×Var says malignant *and* PO-VIMP points at text.

## 3. MVP checklist to tell the professor

1. **Claim:** Magnitude × Variance of PO resolves benign vs malignant shift (2×2).  
2. **Evidence:** 5 boards, pooled contingency `[[5,0],[0,5]]`.  
3. **Estimator honesty:** demo uses $\mathrm{PO}=1-\hat p(Y\mid x)$; PermuCATE is the causal upgrade path (already in-repo).  
4. **Next natural paper beat:** Mag×Var on PermuCATE $\varphi$-risk + multimodal group-permute attribution as the channel drill-down.

## 4. Pointers

| artifact | path |
| --- | --- |
| Mag×Var logic | `docs/method/AS_DTA_magnitude_variance_logic.md` |
| Variance justify | `docs/method/AS_DTA_po_risk_variance_justify.md` |
| 2×2 eval | `docs/method/AS_DTA_mag_var_eval_metrics.md` |
| PermuCATE | `Python/src/permuCATE_vimp.py` |
| Uplift PO monitor | `Python/src/fsds_rl_closed_loop.py` → `uplift_po_monitor` |
| Multimodal attribution | `Python/src/msrvtt_multimodal_attribution.py` |
