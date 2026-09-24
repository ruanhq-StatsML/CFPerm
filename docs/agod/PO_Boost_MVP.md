# PO-Boost MVP (minimal viable prototype)

> Status: **MVP frozen for advisor review.** Further productization follows advisor schedule.  
> Entry: `python3 scripts/po_boost_mvp.py`  
> Branch: `cursor/po-risk-posttrain-latex-abce` · PR #79

## What this MVP is

PO-risk reallocates **post-training update budget** under an Acc constraint.
Two layers (both useful; **sample-level is the concrete OOD→weight story**):

| Layer | What gets a weight | Concrete rule |
|---|---|---|
| **Modality** (optional / α-side) | which tower gets BWD / steps / LR | freeze low chronic \(L_m\); dump steps on spike \(S_m\) |
| **Sample** (core for OOD) | which rows in the batch get CE weight | if batch rejected: \(w_i\propto\sqrt{\mathrm{PO}_i}\); else \(w_i=1\) |

α Softmax is secondary. The advisor-facing claim is: **on OOD-large rows, up-weight the post-train loss**.

## Sample-level OOD → weights (concrete, not a pipeline poem)

**Step 1 — Is this batch an OOD event?** (batch gate, not per-row label)

```text
rejected = OnlineRFPerm/CFPerm flag
        OR hop_oos: e_t / e_{t-1} ≥ γ (default 1.5)
        OR proxy thresholds on mean PO / |ΔPO| / MMD
```

If `rejected=False` → **all** \(w_i=1\) (no sample reweight). Calm windows stay untouched.

**Step 2 — How large is OOD on each row?** (same PO family, observation resolution)

\[
\mathrm{PO}_i \;=\; |Y_i - \hat\mu(X_i)|
\quad\text{(residual under the control fit — not “intrinsic hard label”)}
\]

**Step 3 — Turn row-OOD into a loss weight** (only if rejected)

\[
w_i
=
\mathrm{clip}\!\left(
\frac{\sqrt{\mathrm{PO}_i}}{\mathrm{mean}_j\sqrt{\mathrm{PO}_j}},\;
[0.25,4]
\right)
\]

Code: `po_iptw_weights(po_i, mode="sqrt", rejected=True)`.  
Default **sqrt** (soft); `prop` / `cbrt` are ablations. Raw \(\propto\mathrm{PO}\) overfits the reject batch.

**Step 4 — Where it hits post-training**

Causal: decide at end of window \(t\), apply on Fit of window \(t{+}1\):

\[
\mathcal{L}_{t+1}
=
\frac{1}{n}\sum_i w_i\,\ell\!\big(f_{\theta}(X_i),Y_i\big)
(+\;\text{optional modality freeze/steps on the same Fit})
\]

In compare: `_weighted_ce(logits, y, w_row)` when `used_reject_w`.

## Why this justifies post-training acceleration

Post-train wall-clock is a **stream of mini-windows**. Uniform CE spends equal gradient mass on easy in-support rows and on the rows that actually carry the shift.

1. **Gate first** — only spend reweight budget when the batch is already flagged OOD (event time). Always-on \(\sqrt{\mathrm{PO}}\) regresses calm packs (measured in gated vs always-√ benches).  
2. **Weight the OOD-large rows** — high \(\mathrm{PO}_i\) = large residual under current \(\mu\); those rows are where the next update should bite if you want next-window MSE / Acc to recover.  
3. **Soft sqrt** — keeps mean \(w=1\) (optimizer scale stable) and caps extremes so one outlier cannot dominate the reject batch.  
4. **KPI** — report **next-MSE drop vs uniform \(w\)** on reject windows only; do not average into calm windows. Ship with Acc constraint; do not claim inference latency.

One-line justify:

> Post-training under shift is slow when updates chase easy rows; gated \(\sqrt{\mathrm{PO}}\) moves the same step budget onto high-residual (OOD-large) rows **only when the batch is rejected**, so \(T(\mathrm{Acc}^\star)\) / next-MSE improve without calm-pack tax.

Modality freeze/dump is the **orthogonal** column story (which head pays BWD). α ranking is optional scaffolding; sample weights do not require believing α.

## Run (no Affec / no torch)

```bash
python3 scripts/po_boost_mvp.py
PYTHONPATH=. python3 -m pytest tests/test_po_boost_synthetic.py tests/test_reject_event.py tests/test_po_risk_fuse.py tests/test_po_roi_export.py -q
```

Synthetic headline (\(M{=}5\)): **`po_fuse` ships** — FLOPs≈0.78, \(T(\mathrm{Acc}^\star)=5\) (modality side). Sample side is exercised whenever `reject_sources` shows hop_oos/proxy and `po_iptw_weights` is applied.

## Code surface (keep small)

| Path | Role |
|---|---|
| `agod/po_risk_train.py` | `row_po_residual`, `po_iptw_weights`, modality actuators |
| `agod/reject_event.py` | batch reject: external → hop_oos → proxy |
| `scripts/smoke_po_boost_synthetic.py` / `po_boost_mvp.py` | demo entry |
| `docs/agod/po_posttrain_roi_map.sql` | ROI_B = next_mse_drop on rejected windows |
| `docs/agod/AGOD_PO_Risk_PostTraining.tex` §Logic B | LaTeX writeup |

## Explicitly out of MVP

- Affec cache fill; live RFPerm beyond `external_rejected=`  
- Drill / math-Eval / review-accel human-gate features  
- Inference-latency claims; α-as-primary narrative  

## Ship rule

Modality: \(\mathrm{FLOPs}<1\) and \(\Delta\mathrm{Acc}\ge-0.005\).  
Sample: on rejected windows, next-MSE vs uniform \(w\) should improve; calm windows stay \(w=1\).
