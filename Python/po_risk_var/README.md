# po-risk-var

**Product:** estimate $\mathrm{Var}(u)$ for PermuCATE unit PO-risk $u=(\varphi-\hat\tau)^2$, then triage benign vs malignant LIVE shift.

## Next-step logic (straightforward)

1. Fit $\varphi$, $\hat\tau$
2. $u=(\varphi-\hat\tau)^2$
3. Mag×Var $(\mu,s)$ — shape
4. **Var(u) = CF leaf-size-debiased + leaf merge (min n=30 / 50)** — default
5. $Y$-flip on text integrity boards
6. `triage()` → hold / conformal / VAI_freeze
7. Honesty: `bias_ratio = B_leaf / V_naive`

| Estimator | Role |
| --- | --- |
| Causal Forest debiased + leaf merge min n=30 (n&lt;200) / 50 (n≥200) | **default** |
| BCF posterior via **rpy2** + `pihat` (R `bcf`) | Bayesian path |
| Bootstrap | intuitive fallback / sanity check |

## Install (inside CFPerm checkout)

```bash
pip install -e Python/po_risk_var
po-risk-var logic
po-risk-var justify --quick
po-risk-var products
po-risk-var biz          # incremental business value per method / app
po-risk-var compare --quick   # CF vs BCF vs boot on scenes + coverage
```

Or without install:

```bash
PYTHONPATH=Python/po_risk_var:Python/src python3 -m po_risk_var justify --quick
PYTHONPATH=Python/po_risk_var:Python/src python3 -m po_risk_var products
PYTHONPATH=Python/po_risk_var:Python/src python3 -m po_risk_var biz
PYTHONPATH=Python/po_risk_var:Python/src python3 -m po_risk_var compare --quick
```

## API

```python
from po_risk_var import monitor_batch, triage, run_justify_suite, run_product_apps_suite

out = monitor_batch(X_ref, Y_ref, X_live, Y_live, method="cf", y_flip=0.0)
print(out["decision"])  # call / actuator / reason / bias_ratio

# Bayesian path (rpy2 → R bcf with pihat)
out_bcf = monitor_batch(X_ref, Y_ref, X_live, Y_live, method="bcf")
```

## Justify

- **Synthetic:** quiet / benign_covariate / malignant_sparse / malignant_dense / label_delay
- **Real-world:** HC3 finance·medicine, 20 Newsgroups, AG News, IMDB

Results: `results/as_dta_po_risk_var_apps/`

## Five product applications

`po-risk-var products` → `results/as_dta_po_risk_var_product_apps/`

1. **Dynamic routing & cost optimization** — PO-risk Var(u) as confidence proxy; benign → light/cache/conformal; malignant → escalate. Industrial LIVE boards saved **$22,430**.
2. **Granular canary & segment rollback** — Severity = J + Y-flip per segment; roll back only the sick segment.
3. **Data flywheel / active learning** — drop high Y-flip poison; keep Y-flip=0 + mild J as gold.
4. **Explainable compliance guardrails** — segment risk-event reports (AI Act style).
5. **Agentic cascade-error blocking** — abort multi-step agent when a step’s Mag×Var / Var(u) spikes.

## Business incremental value

`po-risk-var biz` elaborates **what each method adds vs its anti-pattern**:

| Layer | Method | Incremental job |
| --- | --- | --- |
| Estimator | CF leaf-debiased (default) | Avoid leaf-bias false freezes; confirm routing |
| Estimator | BCF rpy2 | Posterior Var(u) for audit / Bayesian stacks |
| Estimator | Bootstrap | Cheap sanity — not a freeze sensor |
| App | Routing | Right-drug $ vs always-heavy / always-light |
| App | Canary | Keep healthy-segment uplift vs global rollback |
| App | Flywheel | Drop poison + keep gold vs ingest-all / drop-all-PO |
| App | Compliance | Auditable freeze vs black-box deny |
| App | Agent | Abort cascade vs run-to-end |

Docs: `docs/method/AS_DTA_po_risk_var_biz_incremental.md`
