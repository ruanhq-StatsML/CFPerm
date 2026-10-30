# PO-risk Var(u) — business incremental value

Incremental business value of PO-risk Var(u) methods (CF / BCF / bootstrap) and five product applications

PO-risk Var(u) is an allocation rule: CF-debiased honesty for LIVE, BCF when auditors need posteriors, bootstrap as sanity — and the five apps turn that sensor into routing dollars ($22430 industrial), micro-loss canaries, clean flywheels, auditable freezes, and aborted agent cascades.

- Industrial LIVE routing save (measured): **$22430**
- Demo sim stack (canary+flywheel+compliance+agent): **$10534**
- Do not sum industrial + sim into one ROI — different ledgers. Pitch industrial $22430 as measured LIVE save; sim stack ($10534) as transparent unit-econ for canary/flywheel/compliance/agent on the demo panel.

## A. Estimator methods — what each one adds

### causal_forest_leaf_debiased  [product_default]

- Baseline: empirical Var(u) or CF V_naive (leaf-size biased)
- Adds: Debias leaf-size inflation in Var(tau) so quiet/benign LIVE does not look malignant; leaf-merge floor min n=30/50 keeps segment monitors stable.
- Ship as: **default LIVE monitor**
  - **false_freeze_avoidance:** Stops leaf-bias masquerading as integrity break → avoid wrong-drug freeze/retrain on synonym / mild covariate boards.
    - Evidence: Synth justify 5/5 + real 4/4 with CF-debiased; honesty bias_ratio flags when debias mattered.
    - $ proxy: Wrong freeze on a benign finance synonym cell alone leaves ~$3.4k/1k on the table (business PR: conformal was the right cheap drug).
  - **confirmed_malignant_routing:** Relative V_debiased↑ confirms Mag×Var / Y-flip → spend freeze+aug dollars only when shape+outcome agree.
    - Evidence: Estimator separation suite: cf_V_debiased / cf_Var_u_tree separate malignant vs benign (hit-rate 1.0).
    - $ proxy: Finance malignant cell: cost/1k $8404→$4067 (~$4.3k) when freeze+aug is correctly selected.
- Not enough when: Need posterior intervals / Bayesian stack → escalate to BCF; need a 30-second sanity check → bootstrap.

### bcf_posterior_rpy2  [bayesian_path]

- Baseline: CF point/debiased Var(u) without posterior draws
- Adds: Full posterior draws of tau via R bcf (rpy2 + pihat) → Var(u) as posterior spread; audit-grade uncertainty for regulated pipelines.
- Ship as: **Bayesian / audit lane**
  - **compliance_uncertainty:** Credible bands on unit risk for auditors: not only 'V is high' but 'posterior mass on integrity break'.
    - Evidence: bcf_Var_u separates malignant/benign (suite ratio ~29×); backend rpy2_bcf with pihat.
    - $ proxy: Pays for itself when a single avoided AI-Act / escalation dispute exceeds MCMC cost (minutes of CPU vs legal/ops review).
  - **bayesian_product_fit:** Drop-in when the org already ships BART/BCF; no second uncertainty language for risk teams.
    - Evidence: Same triage actuators; only Var(u) backend swaps.
    - $ proxy: Integration cost ↓ vs inventing a fourth estimator.
- Not enough when: Latency-critical LIVE default → stay on CF; tiny n without pihat quality → prefer CF + bootstrap sanity.

### bootstrap  [fallback_sanity]

- Baseline: ship CF/BCF with no cross-check
- Adds: Intuitive refit of PermuCATE tau → Winsorized Var(u); catches gross misspecification cheaply.
- Ship as: **fallback / CI sanity**
  - **engineering_confidence:** A/B against CF before a freeze: if boot and CF disagree wildly, hold + investigate instead of auto-escalate.
    - Evidence: boot_Var_u still separates (suite ratio ~15×) but is heavier-tailed — never the product default.
    - $ proxy: One blocked false escalation (GPU retrain / human on-call) typically dwarfs bootstrap CPU.
- Not enough when: Do not use as sole freeze sensor (fourth-moment noise).

## B. Product apps — incremental vs anti-patterns

### 1. Dynamic routing

- Baseline anti-patterns: always_heavy_model (pay GPU on every request), always_light_or_cache (miss integrity breaks → FN ledger), OOD_high_implies_retrain (wrong-drug tax)
- PO-risk adds: Var(u) + Mag×Var + Y-flip as a *confidence / disease* proxy: quiet/benign → light/cache/conformal; malignant → escalate.
- Headline incremental: **$22430**

### 2. Granular canary

- Baseline anti-patterns: global_rollback_on_any_alert (burn healthy segments), global_accuracy_only_monitor (miss segment integrity breaks)
- PO-risk adds: Segment Severity=J+Y-flip + CF Var(u): rollback only sick segments; keep canary elsewhere (micro-loss iteration).
- Headline incremental: **$600**

### 3. Data flywheel

- Baseline anti-patterns: ingest_all_online_feedback (data poisoning), drop_all_high_PO (throws away hard benign / synonym gold)
- PO-risk adds: Y-flip gate kills integrity poison; Y-flip=0 + mild J kept as distribution-expanding gold.
- Headline incremental: **$694**

### 4. Compliance guardrails

- Baseline anti-patterns: black_box_deny_or_silent_pass, global_accuracy_dashboard_only
- PO-risk adds: Segment risk-event narratives: window × segment × Y-break × J × V_debiased → local freeze — AI-Act-readable, not opaque.
- Headline incremental: **$9200**

### 5. Agent cascade block

- Baseline anti-patterns: run_workflow_to_end, only_check_final_answer
- PO-risk adds: Per-step Mag×Var / Var(u): abort + self-correct / human when a step's residual variance spikes — stop cascade amplification.
- Headline incremental: **$40**

## Justify before dollars

- Sensor hit-rates: synth **5/5**, real **4/4**
- Industrial anchor source: `results/as_dta_business_pr (HC3 finance/medicine LIVE)`

```bash
po-risk-var biz
po-risk-var products
po-risk-var justify --quick
```
