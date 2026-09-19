# LOCO PO-risk tip timing monitor

- tips J* = `[0, 1]`
- known shift @ batch **20**
- burn-in = 8, α = 0.05, L̄ = 0.7331831766521688, thr = 0.48682202434310673
- **hard-gate tip change-point t** = **20** (delay=0)
- tip OnlineRFPerm first reject t = 9
- full-PO first reject t = **11** (delay=-9)
- lead tip−PO = `9` (negative ⇒ tip earlier)

Score: `S_t = |L_tip(t) − L̄| + PO_t`; timing = first `S > mean(S_burn)+k·std(S_burn)`.

Artifacts: `results/loco_po_timing/summary.json`, `results/loco_po_timing/traj.csv`.
