# FSDS multi-step overnight iterations

Locked scope: **图谱特征** → multi-step selection only.  
No community / ego / graph-local backends.

```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_fsds_multistep_iterate.py --iter-tag iterNN
```

Results: `results/tencent_gr_fsds_iterate/` (`ITERATION_LOG.md` + per-iter `FINDINGS.md`).

## iter01 headline
- π-stable F / cmean+π beat baseline slightly on W2 HGB (0.7735 vs 0.7679)
- Hard corr-prune@0.92 **hurts** (W2 0.64) — too aggressive on 21 feats
- Rare positives (≈3 in train) make 5-fold fragile — next iters: softer prune, fewer folds, seed sweeps
