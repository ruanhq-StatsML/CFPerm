# Method: online RFPerm + PO-risk

This is gated PO-risk reweighting on an online stream. Not AGOD.
Not DRE. Not always-on IPTW. Not PO-tail subset.

Every observation on the last two batches is scored:

\[
\texttt{po\_risk0}_i = |Y_i - \mu_0(X_i)|
\quad\text{(mix \(0.5\) with the batch gap)}
\]

\[
w_i =
\begin{cases}
1 & \text{quiet, or } T_i=0 \\
\sqrt{\texttt{po\_risk0}_i} & \text{fire and } T_i=1
\end{cases}
\quad\text{then mean 1, clip \((0.05,20)\).}
\]

That instance map is the quantification — not a batch scalar.
`quantify_last_two(...)` returns `po_risk0` and `w` for every row,
plus p10–p90. `run_online_rfperm(..., detail=True)` keeps the
vectors on each hop.

```python
from agod import run_online_rfperm, instance_po_risk, po_iptw_weights, quantify_last_two

rec = run_online_rfperm(stream, gate=1.5, learner="rf", detail=True)
hop = rec["history"][k]
# hop["po_risk0"][i], hop["w"][i]  — one number per observation
# hop["po_t1"]["p50"], hop["w_t1"]["p90"]  — hop summary
```

1. **Probe** — same shallow RF as the IPTW stream (`n_estimators=20`,
   `max_depth=4`) fit on 上一批 as \(T=0\).
2. **PO-risk** — `po_risk0_i = instance_po_risk(Y_i, μ0(X_i), ...)`
   (mix \(0.5\) with the batch gap).
3. **Gate** — consecutive OOS, no K-fold:
   \(e_{\mathrm{now}}/e_{\mathrm{prev}} \ge 1.5\). First hop is
   quiet (no previous OOS). In-sample \(e_1>e_0\) would fire every
   hop on trees — do not use it.
4. **Quiet** — last two batches, \(w=1\).
5. **Fire** — same rows. \(T=0\) stays \(1\);
   \(T=1\) gets `po_iptw_weights(po_risk0, mode="sqrt")` (mean 1).
6. **Anneal** — the next hop's \(e_{\mathrm{prev}}\) is that large
   error, so the ratio drops and the stream cools.

Drop-old (`run_resid_stream`) is a different, harder lever for a
one-shot concept cut. Method LaTeX:
`docs/po_risk/PO_risk_reweight_method.tex`. Boards:
`AGOD_overview.md`, `AGOD_po_refit_gated.md`,
`AGOD_po_refit_real.md`. Tables: `AGOD_performance_tables.tex`.
