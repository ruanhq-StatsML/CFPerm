# Gradient-cosine × modality LR

## Point

在 online-window 流上，用 **梯度余弦相似度** 刻画模态更新几何，并可选地折进 L2 步长：

| 量 | 含义 |
|---|---|
| `cos(g_m, g_shared)` | 模态投影梯度与 head 的对齐 |
| `pair_cos` | 模态两两梯度夹角（<0 ⇒ 冲突） |
| `temporal_cos` | 跨窗 `cos(g_t, g_{t-1})`（步长是否抖） |
| `align_gain` | `½(1+cos_to_shared)` ∈ [0,1] |
| `soft_gradcos` | `LR_m = soft(α_m) · [(1-λ)+λ·align_gain_m]` |

不必 dense equal：α 负责 **谁该大步**；grad-cos 负责 **谁的更新方向不打架**。

## Smoke

| Dataset | scheduler | Acc lift | pair_cos | conflict | temporal_cos | lr_ratio |
|---|---|---:|---:|---:|---:|---:|
| amazon | equal | -0.009 | +0.688 | 0.00 | -0.205 | 1.00 |
| amazon | soft | -0.007 | +0.669 | 0.00 | -0.301 | 3.10 |
| amazon | soft_gradcos | +0.000 | +0.710 | 0.00 | -0.241 | 3.42 |
| msrvtt | equal | +0.010 | +0.190 | 0.06 | +0.302 | 1.00 |
| msrvtt | soft | +0.006 | +0.216 | 0.06 | +0.161 | 2.83 |
| msrvtt | soft_gradcos | -0.006 | +0.215 | 0.06 | +0.413 | 2.82 |

### Readout
- **amazon**: best `soft_gradcos` lift +0.000 (vs equal -0.009, Δ=+0.009); pair_cos=+0.710, conflict=0.00
- **msrvtt**: best `soft` lift +0.006 (vs equal +0.010, Δ=-0.004); pair_cos=+0.216, conflict=0.06

### Other opportunities
- **Conflict damp**: if `cos(g_m,g_shared)<0`, force damp schedule on that mod
- **Temporal gate**: low `temporal_cos` → raise β (flatten) until grads stabilize
- **Budget × align**: redistribute fixed excess-LR by `α_m · align_gain_m`
- **Sample-grain grads**: cosine on per-example grads → replay weights (different actuator)

```bash
PYTHONPATH=. python3 scripts/run_agod_gradcos_lr.py
```
