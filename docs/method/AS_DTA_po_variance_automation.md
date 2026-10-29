# PO-risk variance — evaluate, visualize, automate (+ attack mixture logic)

## 1. Why variance (not mean alone)

Unit risk with labels: $\mathrm{PO}_i = 1 - \hat p(Y_i \mid x_i)$.

Across bootstrap **sub-environments** $b=1..K$ (windows of size $\approx n/5$):

$$
r_b = \frac{1}{|b|}\sum_{i\in b}\mathrm{PO}_i,\quad
\overline{r}=\mathrm{Mean}(r_b),\quad
s=\mathrm{Std}(r_b),\quad
\mathrm{PeakGap}=\max r_b - \mathrm{median}(r_b),\quad
\mathrm{CVaR}_{95}=\mathbb{E}[r_b \mid r_b \ge q_{95}].
$$

| Shape | Mean $\overline{r}$ | Sub-env Std $s$ / Peak / CVaR | Story |
| --- | --- | --- | --- |
| Quiet | ≈ REF | ≈ REF | hold |
| Benign OOD | ↑ | flat / mild | uniform calibration stress → conformal |
| Malignant | ↑↑ | spiked | concentrated failures → freeze-lexical |

**Evaluate:** report $(\Delta\mathrm{Mean},\, \Delta s,\, \Delta\mathrm{Peak},\, \Delta\mathrm{CVaR},\, \Delta\mathrm{Error})$ vs REF.  
**Visualize:** Mean×Std scatter; sub-env strip; unit-PO histogram; Δ bar panel.  
**Separate:** Mean can rise for both; **Std / Peak / Error** pull malignant away.

## 2. Automate (closed loop)

```
LIVE batch
  → unit PO_i (needs Y, or calibrated proxy)
  → bootstrap {r_b} → Mean, Std, Peak, CVaR
  → Δ vs rolling REF
  → gate:
       ΔMean < τ_m  and  ΔError < τ_e     → hold
       ΔMean ≥ τ_m  and  not spiked       → conformal_calibrate
       spiked (Δs / ΔPeak / ΔCVaR / ΔError) → freeze_lexical_adv_aug
  → after actuator: remeasure ΔError / ΔFN / ΔCost  (effectiveness with Y)
```

Spiked rule (MVP defaults, calibrate on a benign probe):

- $\Delta\mathrm{Error} \ge 0.12$ **or**
- $\Delta s \ge 0.018$ **or**
- $\Delta\mathrm{Var}_{\mathrm{unit}} \ge 0.03$ **or**
- $(\Delta\mathrm{CVaR}\ge 0.18 \land \Delta\mathrm{Error}\ge 0.08)$

**Automation tips**

1. Maintain a **benign envelope** from periodic synonym/BT probes; flag LIVE that exceeds envelope on $(s,\mathrm{Peak},\mathrm{Error})$.
2. Never hard-filter top-PO before computing $\{r_b\}$ — that deletes the Peak.
3. Never pseudo-label from PO; Y stays fixed for effectiveness.
4. Log stratum by channel; run the same gate per stratum (api vs open).

Code path: `Python/src/as_dta_benign_malignant_mvp.py` → `_shape_first_route`.

## 3. Why attack as a **concentrated subset** (mixture)

Malignant LIVE is built as:

$$
x_i' =
\begin{cases}
\mathrm{Attack}(x_i) & i \in S,\ |S|\approx n/3 \\
x_i & i \notin S
\end{cases}
$$

**Logic**

1. **Real attacks are sparse.** Evaders hit a minority of traffic (targeted spoof), not every row equally — unlike domain synonym drift.
2. **Heteroscedastic peaks.** Bootstrap windows that over-sample $S$ get high $r_b$; windows that miss $S$ stay near REF → large $s$ and PeakGap. Uniform benign edits raise all $r_b$ together → Mean↑, $s$ quiet.
3. **Separability.** Mixture is the statistical dual of “intentional injection”: mass concentrates on a PO minority (top-share / CVaR / FN on $S$).
4. **If you attack everyone equally,** geometry collapses toward “strong benign / domain shock” (Mean↑↑, flatter relative Peak) — still bad, but the **variance signature weakens**; then lean on $\Delta\mathrm{Error}/\mathrm{FN}$ with Y.

So subset concentration is not a hack — it is the **generative model of rule-breaking traffic** that Peak/CVaR are designed to catch.

## 4. Other attack / rule-break morphologies (yes, more)

| Kind | Mechanism | PO shape | Actuator lean |
| --- | --- | --- | --- |
| **A. Concentrated mixture** (MVP default) | Attack subset $S$ hard | Mean↑, Std/Peak/CVaR↑, FN on $S$ | freeze-lexical + adv-aug |
| **B. Uniform strong evasion** | Attack ≈ all rows lightly/hard | Mean↑↑, Peak flatter | still freeze / retrain head; variance gate weaker — use ΔError |
| **C. Feature-collision / trigger** | Inject opposite-class tokens | Mean↑, often Peak if sparse triggers | freeze lexical (zero surface) |
| **D. Watermark wash** | Paraphrase to kill WM bits | WM gray/neg; PO may stay mild | do **not** equate wash with A; need task Y |
| **E. Threshold gaming** | Stay just under $\tau$ | Mean mild, FN↑ at operating point | conformal + cost-sensitive $\tau$; PR/ROI panel |
| **F. Label / feedback poison** | Corrupt Y or reward | PO vs true Y explodes; monitor Y QC | human review + Y audit (not WM) |
| **G. Cross-domain mimic** | Looks like benign OOD but crafted | Mean like benign; Peak subtle | benign envelope + secondary probe |
| **H. Adaptive multi-hop** | Attacker reacts to your actuator | Peak migrates after freeze | re-diagnose each hop; IV/uplift gate |

**Rule-breaking ≠ only TextFooler.** Anything that **concentrates decision failures** while preserving labels is “malignant-shaped.” Anything that **shifts the whole mass mildly** is “benign-shaped.” Automation classifies **shape**, not attacker brand.

## 5. One-line claim

> Automate AS-DTA as Magnitude→Variance of PO sub-environments; use concentrated-subset attacks as the canonical malignant generator because real evasion is sparse and that sparsity is exactly what Peak/CVaR measure. Other morphologies exist — score them with the same $(\Delta\mathrm{Mean},\Delta s,\Delta\mathrm{Error})$ panel, not a single OOD score.
