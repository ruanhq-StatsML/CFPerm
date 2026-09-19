# FSDS multi-step overnight iterations

Locked scope: **图谱特征** → multi-step selection + **official FSDS** fuse.  
No community / ego / graph-local backends.

Official FSDS: `StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg`  
(fit W1 only; W2 = temporal holdout).

**整理报告：** [`docs/tencent_gr/FSDS_PO_OVERNIGHT_REPORT.md`](../../docs/tencent_gr/FSDS_PO_OVERNIGHT_REPORT.md)  
**PO cookbook：** [`docs/tencent_gr/PO_RISK_FOR_DS.md`](../../docs/tencent_gr/PO_RISK_FOR_DS.md)  
**Summary：** `ITERATION_LOG.md` → OVERNIGHT_SUMMARY

```bash
# DS default
PYTHONPATH=. python3 scripts/tencent_gr/po_help_fsds.py --select-k 15 --seed 0

# Overnight harness
PYTHONPATH=. python3 scripts/tencent_gr/run_fsds_multistep_iterate.py --iter-tag iterNN
```

## Headlines
- **iter01**: π-stable / cmean+π 略好；硬 corr@0.92 伤
- **iter02–03**: 官方 FSDS fuse；收敛 cmean+π→FSDS
- **iter05**: PO-VIMP → FSDS mean W2 **0.722** > Z **0.720** > A **0.716**
- **iter06**: rare-π 最低 σ（0.137）；boot-π 非免费午餐
- **iter07**: α 无关；tight pool=k 伤
- **iter08**: k&lt;15 伤；P_po @k=18 → **0.724**；τ̂²-rows ≈ P_po
- **iter09**: seed-maj2 伤均值；LOO-pos maj σ **0.114**
- **iter10**: CLI + overnight wrap；定时器已停
