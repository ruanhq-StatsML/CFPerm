# PO-risk Var(u) — CF vs BCF vs bootstrap on existing scenes

Side-by-side CF / BCF / bootstrap Var(u) on existing synth + real scenes, plus 90% CI coverage of Δ=μ_LIVE−μ_REF (null Δ=0)

## Talk track

- Synth call agreement across methods: 5/5 boards (same triage call).
- Hit cf: synth 5/5 (1.00), real pairs 4/4 (1.00).
- Hit bcf: synth 5/5 (1.00), real pairs 4/4 (1.00).
- Hit boot: synth 5/5 (1.00), real pairs 4/4 (1.00).
- Coverage cf: quiet P(CI_Δ∋0)=0.00 (nominal 0.90, width≈0.118); malignant power P(CI_Δ∌0)=1.00; μ_LIVE point-cover=0.00.
- Coverage bcf: quiet P(CI_Δ∋0)=0.25 (nominal 0.90, width≈0.776); malignant power P(CI_Δ∌0)=1.00; μ_LIVE point-cover=0.12.
- Coverage boot: quiet P(CI_Δ∋0)=1.00 (nominal 0.90, width≈3.81); malignant power P(CI_Δ∌0)=1.00; μ_LIVE point-cover=0.50.
- Effect reading: Mag×Var + Y-flip dominate the triage call — CF/BCF/boot agree whenever the oracle is clear; CF V_debiased often clips near 0 (honesty: leaf bias removed) while boot Var(u) is heavy-tailed.
- Coverage reading (Δ=μ_LIVE−μ_REF, null Δ=0): quiet P(CI∋0) cf=0.00, bcf=0.25, boot=1.00 (nominal 0.90); malignant power cf=1.00, bcf=1.00, boot=1.00.

**Recommend:** Ship CF leaf-debiased as LIVE default for scene triage (same calls as BCF/boot on existing boards; bias_ratio honesty). For CI coverage of mean-risk shift Δ: bootstrap is the conservative cover (wide); BCF is the audit posterior middle; CF tree-CIs are sharp but undercover quiet — do not use CF draws alone as a Δ confidence statement. Note: coverage uses raw (unmerged) CF draws; leaf-merge stays on for product Var(u) stability.

## Synth boards

| board | oracle | cf | bcf | boot | agree? |
| --- | --- | --- | --- | --- | --- |
| quiet | quiet | quiet✓ | quiet✓ | quiet✓ | yes |
| benign_covariate | benign | quiet✓ | quiet✓ | quiet✓ | yes |
| malignant_sparse | malignant | malignant✓ | malignant✓ | malignant✓ | yes |
| malignant_dense | malignant | malignant✓ | malignant✓ | malignant✓ | yes |
| label_delay | malignant | malignant✓ | malignant✓ | malignant✓ | yes |

## Real-world pairs

| board | cf pair | bcf pair | boot pair | V_gap cf/bcf/boot |
| --- | --- | --- | --- | --- |
| hc3_finance | yes | yes | yes | 0/0.0154/-0.0401 |
| news20_sci | yes | yes | yes | 0/0.00554/-0.278 |
| ag_news | yes | yes | yes | 0/-0.0593/0.103 |
| imdb | yes | yes | yes | 0/0.128/-0.0047 |

## Coverage (90% CI for Δ = μ_LIVE − μ_REF)

| method | quiet P(CI∋0) | quiet width | malignant power | mal width | μ_LIVE point-cover |
| --- | --- | --- | --- | --- | --- |
| cf | 0.00 | 0.1182 | 1.00 | 11.35 | 0.00 |
| bcf | 0.25 | 0.7757 | 1.00 | 7.906 | 0.12 |
| boot | 1.00 | 3.805 | 1.00 | 8.99 | 0.50 |

Quiet: want P(CI_Δ∋0) near nominal 0.90. Malignant: want high power (CI excludes 0).
μ_LIVE point-cover is a sanity check that the μ_LIVE CI covers its point estimate.

```bash
po-risk-var compare --quick
```
