# OnlineRFPerm → T=0/T=1 PO re-fit → post-hoc √PO

## Claim

**Uniform should look slightly best overall.** Re-adjustment only on batches that OnlineRFPerm calls *significant* (small `p`, FDR reject). Weighting is **post-hoc**, after the test.

## Pipeline

```
batch_t arrives
    │
    ├─ OnlineRFPerm: T = MSE(f_ref,t) − E_ref
    │                 p = rank/EWMA vs history; online FDR
    │
    ├─ if NOT reject  →  w = 1  (uniform)  ✓ default
    │
    └─ if reject (small p)  →  POST-HOC PO re-fit
           T=0 : recent control window (batch just before the pair)
           T=1 : previous batch ∪ current batch
           fit μ0 on T=0  (optional μ1 on T=1)
           PO_i = |Y_i − μ0(X_i)|  (+ blend |μ1−μ0|)
           w_current ∝ √PO   on the *current* slice of T=1
           fit downstream RF with sample_weight=w_current
```

## Why T=0 / T=1 this way

| label | data | role |
|---|---|---|
| **T=0** | most recent *pre-pair* batch(es) | control / “what we just believed” |
| **T=1** | prev ∪ cur | treated / “the shifted regime we are adjusting on” |

PO-learner is **re-fit every time the gate opens**, so weights always use the latest cut — not a stale probe residual from the previous step.

## Modes in the runner

| mode | behavior |
|---|---|
| `uniform` | always w=1 |
| `sqrt` | always √PO from probe residual (no gate, no T0/T1 re-fit) |
| `sqrt_gated` | RFPerm gate + old probe residual √PO |
| `sqrt_gated_refit` | RFPerm gate + **T0/T1 PO re-fit** √PO ← main |
| `dre` | logistic density-ratio baseline |

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_rfperm_po_refit.py \
  --datasets metro_interstate beijing_pm25 stocks_AAPL waymo_proxy stocks_MSFT stocks_IWM \
  --batch-size 100 --n-batches 40 --n-burn 5 --n-control 1
```
