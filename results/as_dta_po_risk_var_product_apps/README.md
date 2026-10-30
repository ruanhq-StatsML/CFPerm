# PO-risk Var(u) — five product applications

Five PO-risk product applications on Mag×Var + Y-flip + CF-debiased Var(u)

Sensors: severity=J+Y-flip; Var(u)=CF leaf-debiased (min leaf merge n=30/50); BCF via rpy2; bootstrap fallback

Industrial LIVE cost saved (business PR): **$22430**

For per-method incremental $ / ops vs anti-patterns see `po-risk-var biz` → `AS_DTA_po_risk_var_biz_incremental.md`.

## 1. Dynamic routing & cost optimization

PO-risk variance as confidence proxy: benign → cheap path; malignant → escalate. Industrial LIVE boards saved $22430.

- Simulated save vs all-heavy: **$2660**
- Documented industrial save: **$22430**
- Incremental job: right-drug allocation (conformal vs freeze+aug) — not a bigger model

## 2. Granular canary & segment rollback

Segment-level Severity=J+Y-flip + Var(u): rollback only the sick segment; other segments keep canary (micro-loss iteration).

- Rolled back: **2** / kept canary: **3**
- Incremental job: vs global rollback — retain canary uplift on healthy segments

## 3. Data flywheel / active learning

Gatekeeper for fine-tune/RLHF: filter high Y-flip poison; keep Y-flip=0 + mild J as distribution-expanding gold data.

- Gold keep rate: **0.45**; drop (poison) rate: **0.55**
- Incremental job: vs ingest-all (poison) and vs drop-all-high-PO (lost gold)

## 4. Compliance / audit guardrails

Turn Mag×Var + Y-flip + Var(u) into segment risk-event reports readable by compliance (EU AI Act style).

- Risk events: **2**

- In 2026-W44, segment=medicine detected integrity break (Y-flip=0.32, J=0.90, V_debiased=0.08) → local freeze.
- In 2026-W44, segment=sci detected integrity break (Y-flip=0.28, J=0.70, V_debiased=0.05) → local freeze.

- Incremental job: vs black-box deny — auditable Y-break × segment freeze

## 5. Agentic cascade-error blocking

Embed Var(u)/Mag×Var at each agent step; abort + self-correct before early hallucination / attack cascades to the final act.

- Cascade blocked: **True** at step **3**
- Incremental job: vs run-to-end — abort before downstream tool/LLM waste

```bash
po-risk-var products
po-risk-var biz
```
