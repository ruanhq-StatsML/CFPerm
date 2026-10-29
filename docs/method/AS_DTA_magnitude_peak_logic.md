# Logic: OOD Magnitude × OOD Variance (Peak) — joint read

Looking at **mean OOD / mean PO** alone is under-informative.  
Looking at **Peak / sub-env variance** alone can false-alarm on noise.  
**Together** they form a 2-D diagnosis that maps to actuators.

## 1. Two axes (same PO, different summaries)

Unit risk (with Y): $\mathrm{PO}_i = 1-\hat p(Y_i\mid x_i)$.

Bootstrap sub-env means $r_b$. Define:

| Axis | Estimator | Question |
| --- | --- | --- |
| **Magnitude** | $\mu = \mathbb{E}[\mathrm{PO}_i]$ or $\overline{r}$ | Is risk *elevated* vs REF? |
| **Variance / Peak** | $s=\mathrm{Std}(r_b)$, $\mathrm{PeakGap}=\max r_b-\mathrm{median}(r_b)$, $\mathrm{CVaR}_{95}$ | Is elevation *even* or *spiked*? |

Joint state $=(\Delta\mu,\;\Delta\mathrm{Peak}\text{ or }\Delta s)$ relative to REF.

## 2. Four quadrants (the logic)

```
                 ΔPeak / Δs →
            low (flat)              high (spiked)
         ┌────────────────────┬────────────────────┐
  Δμ     │                    │                    │
  low    │  Q1 QUIET          │  Q2 NOISE / TINY   │
         │  hold              │  monitor / human   │
         │  (no disease)      │  (Peak w/o Mean —  │
         │                    │   often sampling)  │
         ├────────────────────┼────────────────────┤
  Δμ     │  Q3 BENIGN OOD     │  Q4 MALIGNANT      │
  high   │  uniform shift     │  concentrated fail │
         │  → conformal τ     │  → freeze-lexical  │
         └────────────────────┴────────────────────┘
```

**Why joint is informative**

- **Q3 vs Q4** share “OOD high” on magnitude — Peak/variance is the *only* cheap separator before you spend on medicine.
- **Q1 vs Q2**: Peak alone without Mean is usually not worth freeze; avoid over-acting.
- **Peak is the readable face of variance** for ops: one number, same story as $s$ / CVaR when mixture attacks create heteroscedastic $r_b$.

## 3. Decision procedure (ops-ready)

1. Compute $(\Delta\mu,\;\Delta\mathrm{Peak},\;\Delta s,\;\Delta\mathrm{CVaR})$; optionally $\Delta\mathrm{Error}$ if Y exists.  
2. If $\Delta\mu$ quiet **and** Peak quiet → **hold**.  
3. If $\Delta\mu$ high **and** Peak/$s$ flat → **benign** → conformal.  
4. If Peak/$s$/CVaR spiked (esp. with $\Delta\mathrm{Error}$↑) → **malignant** → freeze-lexical.  
5. Remeasure after actuator (effectiveness) — joint panel should move back toward Q1.

Peak is sufficient as the *display* variance proxy; keep $s$ and CVaR as confirmatory in the same panel.

## 4. Information content (one sentence)

> Magnitude answers **whether** the batch is sick; Peak/variance answers **which disease** (uniform calibration vs concentrated evasion) — jointly they are the minimal informative statistic for AS-DTA routing.

## 5. Tie to dollars

Q3 medicine cuts FP (\$50); Q4 medicine cuts FN (\$20). Misreading Q3 as Q4 (or vice versa) is the wrong-drug tax that Mean×Peak exists to avoid.
