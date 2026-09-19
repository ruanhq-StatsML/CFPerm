# FSDS multi-step overnight iterations

Locked scope: **图谱特征** → multi-step selection + **official FSDS** fuse.  
No community / ego / graph-local backends.

Official FSDS (stats method): `StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg`  
(fit W1 only; W2 = temporal holdout).

PO-risk (讲武德): `docs/tencent_gr/PO_RISK_FOR_DS.md` — period W; ranking prior only.

```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_fsds_multistep_iterate.py --iter-tag iterNN
```

Results: `results/tencent_gr_fsds_iterate/` (`ITERATION_LOG.md` + per-iter `FINDINGS.md`).

## Headlines
- **iter01**: π-stable / cmean+π slightly beat baseline on W2; hard corr@0.92 hurts
- **iter02**: fused official FSDS; soft-corr@0.98 best W2; J*-only prefilter hurts; rare-pos AP tiny; 3-fold π unstable
- **iter05**: PO-VIMP → FSDS mean W2 **0.722** > Z **0.720** > A **0.716** (seeds 0/1/2)
- **iter06**: `po_help_select` + rare-pos π; **P_po** still best mean; **Z_PO_rare** lowest σ (0.137); boot-π not free lunch
- **iter07**: α∈{0.3,0.5,0.7} **identical** on d=21; tight pool k hurts; P_po still best mean, rare-π best σ

Timer: every **20 minutes** until morning.
- **iter08**: k&lt;15 hurts; P_po best at **k=18** (0.724); τ̂²-rows ≈ P_po (no clear win)
