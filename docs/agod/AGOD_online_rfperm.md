# Method: online RFPerm + PO-risk

This is the AGOD stream method. Not DRE. Not always-on IPTW. Not
PO-tail subset.

```python
from agod import run_online_rfperm, instance_po_risk, po_iptw_weights

rec = run_online_rfperm(stream, gate=1.5, learner="rf")
# rec["online_mse"]  next-batch score
# rec["fire_rate"]   fraction of hops that heated
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
one-shot concept cut. Boards: `AGOD_overview.md`,
`AGOD_po_refit_gated.md`, `AGOD_po_refit_real.md`.
