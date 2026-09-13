# Continuous soft-LR · sample size · evaluation

## Point (fewer switches)

不需要一堆 L3 gate 变体。核心是：

1. **连续 next-stage 步长**（L2 soft LR）：`LR_m = lr0 · (β + (1-β)·α_m·|M|)`
2. **要多少 adapt samples** 才能相对 equal-LR 出现 Acc lift（`N*`)
3. **不必 dense equal adapter**：用 attribution 加权的连续步长；可选一个 sparse gate 做对照

## How to characterize continuous step sizes

| 量 | 含义 |
|---|---|
| `lr_ratio` = max LR / min LR | 模态步长拉开程度（=1 ⇒ 退化成 equal） |
| `lr_std` / `lr_cv` | 步长离散度 |
| `H(α)` | routing 熵（equal ⇒ log\|M\|；越小越尖） |
| Acc lift(N) | 效用随 sample size 的曲线 |
| `flops_rel` | soft_lr=1；soft_sparse<1（只省 adapt BWD） |

## How to evaluate

- **相对 equal**：同一 N、同一 window 流，比 mean Acc lift
- **N\***：升序 N 曲线上，**第一个** mean Acc lift > 0 的 mean `N_adapt`
- **不是越尖越好**：看 lift 与 `lr_ratio` / `H(α)` 是否同向；尖但 lift≤0 ⇒ 过冲
- **sparse 不是必须**：soft_lr 已是非 equal；soft_sparse 只在要省 adapt FLOPs 时加

## Smoke board

| Dataset | policy | N_cur | N_adapt | Acc lift | flops_rel | lr_ratio | H(α) |
|---|---|---:|---:|---:|---:|---:|---:|
| amazon | equal | 40 | 24 | +0.000 | 1.000 | 1.00 | 0.345 |
| amazon | equal | 80 | 45 | +0.060 | 1.000 | 1.00 | 0.593 |
| amazon | equal | 120 | 67 | +0.050 | 1.000 | 1.00 | 0.563 |
| amazon | equal | 160 | 83 | +0.023 | 1.000 | 1.00 | 0.525 |
| amazon | equal | 200 | 107 | +0.002 | 1.000 | 1.00 | 0.570 |
| amazon | soft_lr | 40 | 24 | +0.021 | 1.000 | 4.95 | 0.397 |
| amazon | soft_lr | 80 | 45 | +0.046 | 1.000 | 2.73 | 0.556 |
| amazon | soft_lr | 120 | 67 | +0.026 | 1.000 | 3.23 | 0.518 |
| amazon | soft_lr | 160 | 83 | +0.031 | 1.000 | 2.51 | 0.578 |
| amazon | soft_lr | 200 | 107 | +0.010 | 1.000 | 3.04 | 0.529 |
| amazon | soft_sparse | 40 | 24 | +0.000 | 0.733 | 4.62 | 0.420 |
| amazon | soft_sparse | 80 | 45 | +0.046 | 0.733 | 3.01 | 0.532 |
| amazon | soft_sparse | 120 | 67 | +0.001 | 0.733 | 3.15 | 0.528 |
| amazon | soft_sparse | 160 | 83 | +0.021 | 0.733 | 3.20 | 0.517 |
| amazon | soft_sparse | 200 | 107 | +0.002 | 0.733 | 2.79 | 0.552 |
| msrvtt | equal | 80 | 52 | -0.024 | 1.000 | 1.00 | nan |
| msrvtt | equal | 160 | 104 | +0.063 | 1.000 | 1.00 | nan |
| msrvtt | equal | 240 | 156 | -0.002 | 1.000 | 1.00 | nan |
| msrvtt | equal | 320 | 208 | +0.057 | 1.000 | 1.00 | nan |
| msrvtt | equal | 480 | 312 | +0.042 | 1.000 | 1.00 | nan |
| msrvtt | soft_lr | 80 | 52 | -0.012 | 1.000 | 2.23 | nan |
| msrvtt | soft_lr | 160 | 104 | +0.059 | 1.000 | 1.75 | nan |
| msrvtt | soft_lr | 240 | 156 | +0.010 | 1.000 | 2.83 | nan |
| msrvtt | soft_lr | 320 | 208 | +0.058 | 1.000 | 1.77 | nan |
| msrvtt | soft_lr | 480 | 312 | +0.037 | 1.000 | 3.14 | nan |
| msrvtt | soft_sparse | 80 | 52 | +0.006 | 0.810 | 2.23 | nan |
| msrvtt | soft_sparse | 160 | 104 | +0.060 | 0.810 | 1.75 | nan |
| msrvtt | soft_sparse | 240 | 156 | +0.004 | 0.810 | 2.83 | nan |
| msrvtt | soft_sparse | 320 | 208 | +0.085 | 0.810 | 1.77 | nan |
| msrvtt | soft_sparse | 480 | 312 | +0.038 | 0.810 | 3.14 | nan |

### N* (first mean Acc lift > 0)

| Dataset | policy | N* | max lift |
|---|---|---:|---:|
| amazon | equal | 45 | +0.060 |
| amazon | soft_lr | 24 | +0.046 |
| amazon | soft_sparse | 45 | +0.046 |
| msrvtt | equal | 104 | +0.063 |
| msrvtt | soft_lr | 104 | +0.059 |
| msrvtt | soft_sparse | 52 | +0.085 |

```bash
PYTHONPATH=. python3 scripts/run_agod_soft_lr_nsize.py
```
