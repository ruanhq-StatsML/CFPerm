# PO-Boost MVP (minimal viable prototype)

> Status: **MVP frozen for advisor review.** Further productization follows advisor schedule.  
> Entry: `python3 scripts/po_boost_mvp.py`  
> Branch: `cursor/po-risk-posttrain-latex-abce` · PR #79

## What this MVP is

PO-risk as the **only sensor family** that reallocates **post-training update budget** under an Acc constraint.

Two orthogonal actuators (do not conflate):

| Dimension | Weight / knob | When |
|---|---|---|
| **Modality / tower** (feature column) | \(\alpha_m\) → freeze / step dump / LR | every window \(t\to t{+}1\) |
| **Sample row** | \(w_i\propto\sqrt{\mathrm{PO}_i}\) | **only** after reject event; calm \(w=1\) |

Short spike \(S=\Delta\mathrm{PO}\) → modality **step dump** (acceleration).  
Chronic \(L=\mathrm{EMA}(\mathrm{PO})\) → modality **freeze** (cut BWD / update variance).  
Reject → row IPTW. Not ε-greedy.

## Run (no Affec / no torch)

```bash
python3 scripts/po_boost_mvp.py
# or
python3 scripts/smoke_po_boost_synthetic.py
PYTHONPATH=. python3 -m pytest tests/test_po_boost_synthetic.py tests/test_reject_event.py tests/test_po_risk_fuse.py tests/test_po_roi_export.py -q
```

Expected headline (synthetic \(M{=}5\)): **`po_fuse` ships** — \(\mathrm{FLOPs}_{rel}\approx0.78\), \(T(\mathrm{Acc}^\star)=5\), Acc held.

## Code surface (keep small)

| Path | Role |
|---|---|
| `agod/po_risk_train.py` | sensors→α→actuators; `po_fuse`; `freeze_flops_rel`; continuous gains |
| `agod/reject_event.py` | reject: external RFPerm → hop OOS → proxy |
| `scripts/smoke_po_boost_synthetic.py` | no-Affec continuous-gains smoke |
| `scripts/po_boost_mvp.py` | one-command MVP entry |
| `configs/agod_po_schedule_cards.json` | \(M\ge3\) / \(M{=}2\) defaults |
| `docs/agod/po_posttrain_roi_map.sql` | ROI / ship gate (business reads SQL) |
| `docs/agod/AGOD_PO_Boost_*.tex` | writeups (scenarios, continuous-time, ε≠Softmax, OOD/FLOPs) |

## Explicitly out of MVP

- Affec cache fill (when data lands)  
- Real OnlineRFPerm wiring beyond `external_rejected=` hook  
- Drill / math-Eval / review-accel human-gate features  
- Inference-latency claims  

## Ship rule

\(\mathrm{FLOPs}<1\) **and** \(\Delta\mathrm{Acc}\ge-0.005\). Else no boost claim (\(M{=}2\) often honest noop).
