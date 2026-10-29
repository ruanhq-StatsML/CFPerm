# Magnitude × Variance — merge logic for benign vs malignant shift

## Scope（先说清楚）

**In scope:** raise *resolution* between two shift types that both look like
“OOD↑” on a single score:

| Shift | Generative story | Wanted actuator |
| --- | --- | --- |
| **Benign shift** | Uniform domain / synonym / paraphrase | `conformal_calibrate` (τ) |
| **Malignant shift** | Concentrated evasion / rule-break on a subset | `freeze_lexical_adv_aug` |

**Out of scope for this merge:** watermark privacy policy, full model SOTA chase,
“is there any OOD at all” as a binary alarm. Those use other tools. Here the
job is only: **same Magnitude band → still tell benign from malignant**.

Peak / CVaR may appear as *derived* summaries of the same bootstrap draws; they
are **not** the primary axis (Peak alone misleads). Primary second axis = **OOD /
PO-risk variance** across sub-environments.

## Why merge Magnitude and Variance

Either axis alone is coarse:

- **Magnitude only** ($\Delta\mu$ of PO / OOD): both shifts raise $\mu$ →
  resolution collapses; you cannot allocate FP-drug vs FN-drug.
- **Variance only** ($\Delta s=\Delta\mathrm{Std}(r_b)$): spikes without
  Magnitude can be sampling noise → false malignant.

**Merge** $=(\Delta\mu,\,\Delta s)$ (joint), or a scalar index for dashboards:

$$
J = \Delta\mu \cdot \bigl(1 + \kappa\cdot \Delta s_{+}\bigr)
\quad\text{with}\quad
\Delta s_{+}=\max(\Delta s,0)/\mathrm{Std}_{\mathrm{REF}}.
$$

- Benign: $\Delta\mu>0$, $\Delta s\approx 0$ → $J\approx\Delta\mu$ (mild).
- Malignant: $\Delta\mu>0$, $\Delta s\gg 0$ → $J$ amplified by variance.
- Quiet: both ≈0 → $J\approx 0$.

Justification: Magnitude supplies **presence** of shift; Variance supplies
**heterogeneity** of shift. Product / joint read is the minimal statistic that
*resolves* the two diseases that share an OOD alarm.

## Joint read (ops)

```
                 Δ Variance (sub-env Std of PO) →
            flat                         high
         ┌──────────────────┬──────────────────────┐
  Δμ     │ QUIET            │ noise / monitor      │
  low    │ hold             │ (var w/o magnitude)  │
         ├──────────────────┼──────────────────────┤
  Δμ     │ BENIGN SHIFT     │ MALIGNANT SHIFT      │
  high   │ conformal τ      │ freeze-lexical       │
         └──────────────────┴──────────────────────┘
```

Optional with Y: add $\Delta\mathrm{Error}/\mathrm{FN}$ as a confirmatory
third read — still variance-led for shape, not Peak-led.

## One-line claim

> Scope = higher-resolution separation of **benign vs malignant shift**.  
> Method = merge **OOD Magnitude × OOD Variance** (Peak is not the headline).  
> Value = right actuator (FP ledger vs FN ledger) when a single OOD score would conflate them.
