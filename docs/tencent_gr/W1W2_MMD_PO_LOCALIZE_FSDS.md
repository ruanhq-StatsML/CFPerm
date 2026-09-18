# W1 vs W2: MMD + PO-risk + conditional-mean → subset → FSDS

换 dataset 跑 **FSDS 标准流**；最上面加 **Standardization**。无 network / 无 SMD 堆指标。

## Protocol
0. **StandardScaler** fit on W1（pipeline 最上面；localize + FSDS 共用）
1. `feature_engineer(t_start,t_end)` ×2，gap ≥ 30d
2. **Subset localization**（item）：rank-average of
   - conditional-mean ‖μ_W2−μ_W1‖（standardized space）
   - RBF-MMD²
   - PO-risk 聚集 mean(τ̂²)
3. Viz localized subset
4. **FSDS**：StandardScaler → VarianceThreshold → SelectKBest → HGB/LogReg → feature ranking  
   （W1-train fit；W2 = temporal holdout）

## Run
```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py \
  --root data/tencent_subset --max-users 20000 --gap-days 30 --localize-k 200
```

Outputs: `results/tencent_gr_w1w2_mmd_po_fsds/`
