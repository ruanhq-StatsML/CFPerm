# 滑动窗候选特征库（算力够就挂上去）

结论和预期一样：算力够的时候，在滑动窗上维护一库统计量，当监控的 **candidate 特征**（以后也可当 clever covariate）。online PCA 这类推理很快。特征都在同一口 \(D_{\mathrm{ref}}\) 钟上，看板会直观。

MMD 仍是 \(\mathrm{MMD}^2(X_{\mathrm{window}},X_{\mathrm{ref}})\)，不是对历史 pairwise。窗是为了滚动统计，不是换参考分布。

---

## 表 1 · 特征库（都对冻结的 \(D_{\mathrm{ref}}\)）

| 特征 | 算什么 | 成本 | 直观读法 |
|---|---|---|---|
| `mean_l2` | \(\|\bar x_W-\bar x_{\mathrm{ref}}\|\) | \(O(p)\) | X 的位置 |
| `pca_recon_excess` | 窗在冻结 \(\mathrm{PCA}_{\mathrm{ref}}\) 上的重构残差 − 基准。滑动窗用 \(C\leftarrow C\pm X^\top X\) 增量更新就够 | \(O(npk)\) / 增量 \(O(p^2)\) | \(P(X)\) 有没有离开参考子空间 |
| `pca_recon_group[g]` | 每个模态自己一块 PCA recon | 同样便宜 | 哪一座塔的 X 在走 |
| `pca_subspace_gap` | 滑动 \(C=X^\top X\) 的 top-k vs \(\mathrm{PCA}_{\mathrm{ref}}\) | 便宜 | 几何变了没有 |
| `mmd_vs_ref` | 核 MMD vs 参考窗 | 中 | 和看板同一口径 |
| `rfperm_T` / `brier_excess` | 冻住 RF / 现模型 | 便宜 | **服务误差**、WHEN |
| `po_risk` / LOGO \(\pi\) | 可选 | 重 | **份额**；算力够再每窗打 |

`SlidingWindowBank.vector()` 就是 clever covariate：recon、T、Brier、分组 recon 拼成一条，进 shadow / 监控模型。和之前 talk 里的 rolling statistics / clever covariate 是同一件事。

---

## 表 2 · 挂上之后该一眼看到的（符合预期）

| DGP | T / Brier | PCA recon / 分组 PCA | 读法 |
|---|---|---|---|
| 渐进 concept（\(P(X)\) 不动，video \(\beta\) 转） | 后期抬（T 0.07→0.24，Brier excess →0.22） | 全局 recon 不主导；分组 PCA 乱、不锁塔 | 服务误差在走，X 子空间没走 → concept |
| 渐进 covariate（audio 均值走） | T/Brier **不抬** | 全局 recon 0.03→0.29；**PCA audio 0.08→0.59**，argmax=audio | X 在 audio 上离开参考子空间，房租还在 → 别按 concept 冻顶 |

分组 PCA 就是「算力够时」不必每窗打 LOGO 也能先定位 \(P(X)\)。份额（PO-LOGO）仍然留给 \(P(Y\mid X)\)；两套指针不要混。

---

## 表 3 · 和四格、latency 怎么接

| 特征库看到 | 对应四格 | 下一窗 |
|---|---|---|
| T/Brier 走、分组 PCA 安静 | 服务响、X 份额安静 | watch / 准备更新；**不要**按 PCA 去 `train_stem` |
| 分组 PCA 走、T/Brier 安静 | X 指针响、房租还在 | `train_stem` 候选在那一组；先不冻顶 |
| 两头都走 | both | 服务确认 + 指针指向 |
| 都安静 | quiet | 全量 |

online PCA 每窗毫秒级，不占 LOGO 那种 500ms。算力够：窗上 **先** 打这库，LOGO/PO 降频。算力紧：只留 T、Brier、`pca_recon_group`。

图：`results/sliding_window_bank/bank_walk.png`。代码：`Python/src/sliding_window_bank.py`。
