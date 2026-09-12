# 模态梯度相关如何指导下一步学习率

同一件几何事实 `R_ij = cos(g_i, g_j)`，用 **ensemble learning** 和 **统计方法** 各说一遍，然后合成下一步的 `LR_m`.

## 设定

窗口 `t` 上，模态 `m ∈ M` 给出（common-dim）梯度签名 `g_m`。下一步融合更新是加权 ensemble：

```
Δθ_{t+1} = −η_{t+1}  Σ_m π_m g_m ,     π ≥ 0,  1ᵀπ = 1
```

`α_m` 是 attribution 路由质量（谁该动），`R` 是梯度相关（他们的票相不相干）。下一步 LR 必须同时吃这两口：

```
π_gls ∝ relu(R̃^{-1} α)
γ     = clip(1 − N_eff/|M|, 0, 1)      # independent → 0, collinear → 1
π*    = (1−γ) α + γ π_gls              # empirical-Bayes stacking mix
LR_m  = η_{t+1} · (β + (1−β) |M| π*_m) · c_m
η_{t+1} = √( (1/|M|) / (π*ᵀ R̃ π*) ) · φ_t
```

`R̃` 是 PSD + equicorrelation shrinkage 之后的 Gram；`φ_t = temporal_gain(cos(g_t, g_{t-1}))`；`c_m` 是 `cos(g_m, g_shared)<0` 时的 conflict damp。FWD 始终开，只改下一窗的 adapt 步长。

---

## Ensemble 视角：模态是选民，不是独立基学习器

把每个模态的 `g_m` 当成一个 voter / base-learner 的下降方向。

### 1. Error–ambiguity（Krogh–Vedelsby）

Ensemble 风险 `E_ens = Ē − Ā`。ambiguity `Ā` 随 pairwise 相关上升而消失：

```
Ā ≈ (1 − ρ̄₊)_+
```

`ρ̄ → 1` ⇒ 选民共线 ⇒ 再给每个模态相同的 LR 是在重复买同一票。有效选民数不是 `|M|`，而是

```
N_eff^{Kish} = |M| / (1 + (|M|−1) ρ̄) = |M|² / (1ᵀ R 1)
N_eff^{GLS}  = 1ᵀ R^{-1} 1     (≥ Kish；相关非均匀时 GLS 更有效)
```

**指导下一步：** `η_{t+1}` 按 `√(N_eff / |M|)` 收缩，使融合步的方差稳住在独立等权基线 `Var = 1/|M|`。这就是 variance-stabilising step。

### 2. 多样性–精度分解（连续版 leader / diversifier / redundant）

| 角色 | 离散（`soft_decorr`） | 连续（本文 `π*`） |
|---|---|---|
| leader | unique mass 最大，占住共享方向 | `π*_m` 大：高精度票 |
| diversifier | 相对 leader 的残差 `1−cos₊` 高 | partial uniqueness `1/P_{mm}` 高 |
| redundant | 与 leader 共线 | `π*` 被 clip 到 0 / 地板 |

等权 + 等相关时，**统计量只缩全局 `η`、不改分配**（选民可交换）；**离散角色会任意点一个 leader**。`α` 一旦不均，两种观点合流：高相关把质量集中到 `argmax α`。

### 3. Stacking / 最小方差组合

线性 ensemble 的最优权（非负约束前）就是最小方差投资组合 / stacking：

```
min_π  πᵀ R π    s.t. πᵀ α 最大化   ⇒   π* ∝ R^{-1} α
```

这是 **SNR 匹配滤波器**：`(πᵀα)² / (πᵀ R π)`。Negative correlation learning 的合法部分落在 `R_ij<0` 且仍与 `g_shared` 对齐——那是真多样性。`R_ij<0` 且 `cos(g_m, g_shared)<0` 是 **冲突**，不是多样性：conflict damp，而不是加权。

---

## 统计视角：相关分数、BLUE、收缩、检验

### 1. 相关 Gram 与收缩

`R_ij = cos(g_i, g_j)`（签名已按 common-dim 对齐；中心化后即 Pearson）。`|M|∈{2,3}` 时样本 Gram 极噪：

```
R_psd = Higham(R)                 # clip 负特征值（冲突子空间）
R̃    = (1−λ) R_psd + λ T         # T = equicorrelation((1−ρ̄)I + ρ̄ 11ᵀ)
λ     = clip( log(cond) / log(κ*), 0, λ_max )
```

Equicorrelation 先验 = 单因子模型（模态共享一个 common descent factor）。`λ` 随条件数升高——这就是 Ledoit–Wolf / Schäfer–Strimmer 在小 `M` 上的可用形式。

### 2. Gauss–Markov / GLS

观察 `g_m = α_m μ + ε_m`, `Cov(ε)=σ² R`。`μ` 的 BLUE / 匹配滤波权就是 `R^{-1} α`。`soft_gls` 直接用这组权；`soft_stat` 再用 `γ=1−N_eff/|M|` 把它们往独立选民先验 `α` 上收缩，并乘方差稳定 `η` 和时序/冲突增益。`|M|∈{2,3}` 时纯 GLS 会过早把弱模态 clip 到 0，混合是必要的。

`R̃ ≈ I` 时 `γ≈0`、`π* ≈ α`、`η ≈ 1` → **退回 soft LR**。`R̃ ≈ 11ᵀ` 时 `N_eff → 1`、`γ→1`、`η → 1/√|M|`，且 `π*` 集中到高于平均的 `α`。

### 3. 偏相关 / 精度矩阵

`P = R̃^{-1}`。模态 `m` 对他人回归后的残差方差是 `1/P_{mm}`（uniqueness）。偏相关

```
ρ_{ij|rest} = −P_{ij} / √(P_{ii} P_{jj})
```

uniqueness 高 = 真正的 diversifier；低 = 冗余票，不该再领大步长。

### 4. Fisher z：相关够不够大，才配改 LR

```
z = artanh(ρ̄),   SE = 1/√(n−3)
```

`n` 取签名长度（单窗 Gram；梯度签名 dim 往往 ≫ 窗数）或窗数（轨迹 ρ̄）。6 个 online window 的 Fisher z **检验力不够**：Amazon `ρ̄≈0.69` 的轨迹级 p 约 0.07，点估计和 `N_eff≈1.2` 已经够用来改 LR，但还不能称为“显著结构相关”。MSR-VTT `ρ̄≈0.21` 则明确更接近独立选民。

只有点估计足够大（`N_eff` 明显小于 `|M|`）时，decorr / GLS 集中才值得开；否则应接近 equal/soft。这避免把一次 noisy cosine 当成结构相关。

### 5. 时序 AR(1)

`φ_t = mean_m cos(g_m,t, g_m,t−1)`：

- `φ → +1`：方向稳，放大 trust region（`η` 上至 cap）
- `φ → −1`：振荡，收缩步长（Polyak 式）

### 6. 冲突 vs 多样性

`R` 非 PSD 的负特征值 = 冲突子空间，先投影再求逆。`cos(g_m, g_shared)<0` 的模态乘 `c_m ∈ [floor, 1]`，避免 PCGrad 意义上的对头更新领大 LR。

---

## 合成扫描（等权 α）

| ρ | N_eff | η | ambiguity | LR ratio | decorr |
|---|---:|---:|---:|---:|---:|
| -0.15 | 4.29 | 1.20 | 1.00 | 1.00 | n |
| +0.12 | 2.40 | 0.89 | 0.88 | 1.00 | n |
| +0.40 | 1.67 | 0.75 | 0.60 | 1.00 | n |
| +0.67 | 1.28 | 0.65 | 0.32 | 1.00 | Y |
| +0.95 | 1.03 | 0.59 | 0.05 | 1.00 | Y |

## 合成扫描（peaked α = 0.70 / 0.20 / 0.10）

| ρ | π video | π audio | γ | LR_stat v | LR_stat a | LR_soft v | roles v/a |
|---|---:|---:|---:|---:|---:|---:|---:|
| -0.15 | 0.70 | 0.10 | 0.00 | 1.67 | 0.31 | 1.99 | d/d |
| +0.12 | 0.73 | 0.08 | 0.20 | 1.51 | 0.23 | 1.99 | d/d |
| +0.40 | 0.83 | 0.06 | 0.44 | 1.49 | 0.16 | 1.99 | d/d |
| +0.67 | 0.87 | 0.04 | 0.57 | 1.47 | 0.13 | 1.99 | l/d |
| +0.95 | 0.90 | 0.03 | 0.66 | 1.46 | 0.11 | 1.99 | l/r |

ρ 低：γ≈0，`π ≈ α`，stat ≈ soft（独立选民）。ρ 高：γ↑，`π` 向 GLS 集中到 video，audio 的下一步 LR 被统计权+全局 `η` 同时压下。

## 因子模型

- **A**（匹配滤波）：`g_m = α_m μ + ε_m`，指标是融合方向与真 μ 的 **cosine**（simplex 权不解偏尺度，所以不用生 MSE）。
- **B**（异方差投票）：`g_m = μ + ε_m`, `Var(ε_m)∝1/α_m`，指标是 BLUE 组合相对 μ 的 MSE。

| ρ | N_eff | cos equal | cos stat | Δcos GLS−eq | MSE eq (B) | MSE BLUE (B) |
|---|---:|---:|---:|---:|---:|---:|
| 0.00 | 3.00 | 0.335 | 0.412 | +0.077 | 1.825 | 0.999 |
| 0.17 | 2.25 | 0.299 | 0.398 | +0.103 | 2.362 | 1.235 |
| 0.33 | 1.80 | 0.276 | 0.392 | +0.125 | 2.832 | 1.376 |
| 0.50 | 1.50 | 0.242 | 0.384 | +0.160 | 3.287 | 1.420 |
| 0.67 | 1.28 | 0.236 | 0.387 | +0.172 | 3.814 | 1.412 |
| 0.84 | 1.12 | 0.219 | 0.382 | +0.188 | 4.288 | 1.402 |

低相关时 cosine 上 GLS/stat 优于等权。ρ 升高时等权方向因多样性崩塌而变差，GLS 仍盯住高 α 票，Δcos 拉大。Model B 的 BLUE 方差始终低于等权（异方差选民的 Gauss–Markov）。

## Amazon / MSR-VTT 轨迹 replay

用已有 grad-cos 窗口的 `pair_cos / α / temporal_cos` 反推 `N_eff, η, π*`（不重训，刻画 *guidance* 本身）。

| traj | ρ̄ | N_eff | η | γ | ambiguity | Fisher p | ρ>0 sig |
|---|---:|---:|---:|---:|---:|---:|---:|
| amazon:equal | +0.688 | 1.19 | 0.63 | 0.41 | 0.31 | 0.072 | n |
| amazon:soft | +0.669 | 1.20 | 0.61 | 0.40 | 0.33 | 0.081 | n |
| amazon:soft_gradcos | +0.710 | 1.17 | 0.62 | 0.41 | 0.29 | 0.062 | n |
| msrvtt:equal | +0.190 | 2.27 | 0.81 | 0.25 | 0.81 | 0.370 | n |
| msrvtt:soft | +0.216 | 2.19 | 0.76 | 0.27 | 0.78 | 0.352 | n |
| msrvtt:soft_gradcos | +0.215 | 2.17 | 0.83 | 0.28 | 0.78 | 0.352 | n |

经验对照：Amazon 窗间 `ρ̄` 高、`N_eff≈1.2` → 统计量会 **降全局步长并（在 γ 升高后）集中 π***；MSR-VTT `ρ̄` 低、`N_eff≈2.2`（接近 `|M|=3`）→ 接近独立选民，γ 小，soft 即可，强行 GLS 集中没有统计理由。轨迹级 Fisher z 因 T=6 检验力不足，应读 `N_eff` / `ρ̄` 而不是 p 值。

## 调度器

| name | 做什么 |
|---|---|
| `equal` | 稠密等 LR |
| `soft` | `α → LR_m` |
| `soft_decorr` | 离散角色增益（高相关才激活） |
| `soft_gls` | 只用 `π* ∝ R̃^{-1} α` |
| `soft_stat` | `π* = (1-γ)α + γ GLS` × 方差稳定 `η` × 时序 × 冲突 |

```bash
PYTHONPATH=. python3 scripts/run_agod_grad_corr_stat.py
PYTHONPATH=. python3 -m tests.test_agod_grad_corr_stat
```
