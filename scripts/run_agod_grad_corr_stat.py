#!/usr/bin/env python3
"""Characterize how modality-gradient correlation guides next-stage LR.

Ensemble view: modalities are voters; high ρ collapses diversity → GLS stacking
π* ∝ R^{-1} α (continuous leader / diversifier / redundant).

Statistical view: Gauss–Markov / matched filter, N_eff, variance-stabilising
η, Fisher-z on ρ, partial uniqueness, AR(1) temporal gain.

  PYTHONPATH=. python3 scripts/run_agod_grad_corr_stat.py
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.ensemble_decorr import soft_decorr_lr
from agod.grad_corr_stat import (
    characterize_stat_traj,
    shrink_correlation,
    correlation_from_pairs,
    effective_ensemble_size,
    stat_corr_lr,
)
from agod.lr_controller import alpha_to_lr

OUT = ROOT / "results" / "agod_grad_corr_stat"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_grad_corr_stat")
GRADCOS = ROOT / "results" / "agod_gradcos_lr" / "agod_gradcos_lr.json"

MODS3 = ("video", "text", "audio")
SEED = 2026


def _iso(c: float, mods=MODS3) -> dict[str, float]:
    return {
        f"{a}|{b}": float(c)
        for i, a in enumerate(mods)
        for b in mods[i + 1 :]
    }


def isotropic_sweep(mods=MODS3, n: int = 25):
    """ρ-grid: N_eff, η, GLS π, LR vs soft / decorr (equal + peaked α)."""
    equal = {m: 1.0 / len(mods) for m in mods}
    peaked = {m: float(v) for m, v in zip(mods, (0.70, 0.20, 0.10))}
    rows = []
    for rho in np.linspace(-0.15, 0.95, n):
        pair = _iso(rho, mods)
        for name, alpha in (("equal", equal), ("peaked", peaked)):
            packed = stat_corr_lr(
                alpha, pair, mods, use_temporal=False, use_conflict=False
            )
            gls_only = stat_corr_lr(
                alpha, pair, mods, use_temporal=False, use_conflict=False, gls_only=True
            )
            soft = alpha_to_lr(alpha, mods, beta=0.10)
            decorr = soft_decorr_lr(alpha, pair, mods)
            st = packed["stats"]
            rows.append(
                {
                    "rho": float(rho),
                    "alpha_kind": name,
                    "mean_rho": st["mean_rho"],
                    "n_eff": st["n_eff"]["n_eff"],
                    "n_eff_kish": st["n_eff"]["n_eff_kish"],
                    "eta_global": packed["eta_global"],
                    "ambiguity": st["ambiguity"],
                    "cond": st["shrink"]["cond"],
                    "lam": st["shrink"]["lam"],
                    "gamma_gls": packed["gamma_gls"],
                    "pi": packed["pi"],
                    "lr_stat": {m: packed["lr"][m] for m in mods},
                    "lr_gls": {m: gls_only["lr"][m] for m in mods},
                    "lr_soft": {m: soft[m] for m in mods},
                    "lr_decorr": {m: decorr["lr"][m] for m in mods},
                    "roles": decorr["decomp"]["roles"],
                    "leader": decorr["decomp"]["leader"],
                    "decorr_active": decorr["decomp"]["decorr_active"],
                    "fisher_p": st["fisher_z"]["p_value"],
                    "lr_ratio_stat": float(
                        max(packed["lr"][m] for m in mods)
                        / max(min(packed["lr"][m] for m in mods), 1e-9)
                    ),
                }
            )
    return rows


def factor_mse_sweep(mods=MODS3, d: int = 64, n_trials: int = 400):
    """Fused-direction quality under two generative models.

    A. Matched filter: g_m = α_m μ + ε_m, Corr(ε)=R.
       Metric = cosine(πᵀg, μ)  (scale-free; simplex weights are not de-biased).
    B. Heteroscedastic votes: g_m = μ + ε_m, Var(ε_m)∝1/α_m, Corr=R.
       Metric = MSE of BLUE-like combo vs μ.
    """
    rng = np.random.default_rng(SEED)
    peaked = np.array([0.70, 0.20, 0.10], float)
    peaked = peaked / peaked.sum()
    rows = []
    for rho in np.linspace(0.0, 0.92, 12):
        pair = _iso(rho, mods)
        r = shrink_correlation(correlation_from_pairs(pair, mods), lam=0.08)["R_shrink"]
        alpha = {m: float(peaked[i]) for i, m in enumerate(mods)}
        packed = stat_corr_lr(alpha, pair, mods, use_temporal=False, use_conflict=False)
        pi = np.array([packed["pi"][m] for m in mods], float)
        pi_gls = np.array([packed["pi_gls"][m] for m in mods], float)
        eq = np.full(len(mods), 1.0 / len(mods))
        evals, evecs = np.linalg.eigh(r)
        evals = np.clip(evals, 1e-8, None)
        L_r = evecs * np.sqrt(evals)
        # Model B covariance: Σ = D^{1/2} R D^{1/2}, D=diag(1/α)
        dstd = 1.0 / np.sqrt(np.clip(peaked, 1e-6, None))
        sigma = (dstd[:, None] * r) * dstd[None, :]
        evals_s, evecs_s = np.linalg.eigh(0.5 * (sigma + sigma.T))
        evals_s = np.clip(evals_s, 1e-8, None)
        L_s = evecs_s * np.sqrt(evals_s)
        # BLUE of common mean: π ∝ Σ^{-1} 1
        prec_s = evecs_s * (1.0 / evals_s)
        prec_s = prec_s @ evecs_s.T
        w_blue = prec_s @ np.ones(len(mods))
        w_blue = np.clip(w_blue, 0.0, None)
        w_blue = w_blue / max(float(w_blue.sum()), 1e-12)

        cos_eq, cos_stat, cos_gls = [], [], []
        mse_eq, mse_blue, mse_stat = [], [], []
        for _ in range(n_trials):
            mu = rng.normal(size=d)
            mu = mu / max(float(np.linalg.norm(mu)), 1e-12)
            # A: signal-scaled + correlated noise
            z = rng.normal(size=(len(mods), d))
            g_a = peaked[:, None] * mu[None, :] + 0.20 * (L_r @ z)
            c = lambda v: float(
                np.dot(v, mu) / max(np.linalg.norm(v) * np.linalg.norm(mu), 1e-12)
            )
            cos_eq.append(c(eq @ g_a))
            cos_stat.append(c(pi @ g_a))
            cos_gls.append(c(pi_gls @ g_a))
            # B: common μ, heteroscedastic correlated votes
            z2 = rng.normal(size=(len(mods), d))
            g_b = mu[None, :] + L_s @ z2
            mse_eq.append(float(np.mean((eq @ g_b - mu) ** 2)))
            mse_blue.append(float(np.mean((w_blue @ g_b - mu) ** 2)))
            mse_stat.append(float(np.mean((pi @ g_b - mu) ** 2)))
        rows.append(
            {
                "rho": float(rho),
                "n_eff": effective_ensemble_size(r)["n_eff"],
                "gamma": packed["gamma_gls"],
                "cos_equal": float(np.mean(cos_eq)),
                "cos_stat": float(np.mean(cos_stat)),
                "cos_gls": float(np.mean(cos_gls)),
                "mse_equal": float(np.mean(mse_eq)),
                "mse_blue": float(np.mean(mse_blue)),
                "mse_stat": float(np.mean(mse_stat)),
                "pi": {m: float(pi[i]) for i, m in enumerate(mods)},
            }
        )
    return rows


def replay_gradcos(path: Path) -> dict:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    traj = payload.get("trajectory") or {}
    out = {}
    for key, rows in traj.items():
        if ":" not in key or not rows:
            continue
        ds, sched = key.split(":", 1)
        mods = list(rows[0].get("alpha") or {})
        if not mods:
            continue
        packed_rows = []
        for r in rows:
            pair = r.get("pair_cos") or {}
            if not pair:
                mpc = float(r.get("mean_pair_cos", 0.0))
                pair = _iso(mpc, mods)
            alpha = r.get("alpha") or {m: 1.0 / len(mods) for m in mods}
            packed = stat_corr_lr(
                alpha,
                pair,
                mods,
                temporal_cos=r.get("temporal_cos"),
                cos_to_shared=r.get("cos_to_shared"),
            )
            packed_rows.append(
                {
                    "t": r.get("t"),
                    "mean_pair_cos": r.get("mean_pair_cos"),
                    "acc_lift": r.get("acc_lift"),
                    "pi": packed["pi"],
                    "eta_global": packed["eta_global"],
                    "n_eff": packed["stats"]["n_eff"]["n_eff"],
                    "ambiguity": packed["stats"]["ambiguity"],
                    "roles": packed["roles"],
                    "lr_stat": {m: packed["lr"][m] for m in mods},
                    "lr_logged": r.get("lr_mult"),
                }
            )
        summary = characterize_stat_traj(rows, mods)
        summary["dataset"] = ds
        summary["source_scheduler"] = sched
        summary["windows"] = packed_rows
        lifts = [float(x["acc_lift"]) for x in rows if x.get("acc_lift") is not None]
        if lifts:
            summary["mean_acc_lift_source"] = float(np.mean(lifts))
        # correlate logged pair_cos with would-be η (guidance tightness)
        rhos = [float(x.get("mean_pair_cos") or 0.0) for x in rows]
        etas = [float(x["eta_global"]) for x in packed_rows]
        if len(rhos) >= 3 and np.std(rhos) > 1e-8 and np.std(etas) > 1e-8:
            summary["corr_rho_eta"] = float(np.corrcoef(rhos, etas)[0, 1])
        out[key] = summary
    return out


def plot_board(sweep, mse_rows, replay, path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.2), facecolor="#f7f5f1")
    fig.suptitle(
        "Modality-grad correlation → next LR  (ensemble GLS + statistical scale)",
        fontsize=13,
        fontweight="bold",
    )

    eq = [r for r in sweep if r["alpha_kind"] == "equal"]
    pk = [r for r in sweep if r["alpha_kind"] == "peaked"]
    xs = [r["rho"] for r in eq]

    ax = axes[0, 0]
    ax.plot(xs, [r["n_eff"] for r in eq], lw=2.2, label=r"$N_{\mathrm{eff}}$ (GLS)")
    ax.plot(xs, [r["n_eff_kish"] for r in eq], lw=1.4, ls="--", label=r"$N_{\mathrm{eff}}$ (Kish)")
    ax.plot(xs, [r["eta_global"] for r in eq], lw=2.0, label=r"$\eta_{t+1}$ (var-stab)")
    ax.plot(xs, [r["ambiguity"] for r in eq], lw=1.4, ls=":", label="ambiguity $1-\\rho_+$")
    ax.axhline(1.0, color="#999", ls="--", lw=0.7)
    ax.set_title("Independent voters → collinear ensemble")
    ax.set_xlabel(r"pairwise $\cos(g_m, g_{m'})$")
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0, 3.2)

    ax = axes[0, 1]
    ax.plot(xs, [r["lr_stat"]["video"] for r in pk], lw=2.2, label="stat LR video (α=0.70)")
    ax.plot(xs, [r["lr_stat"]["audio"] for r in pk], lw=2.2, label="stat LR audio (α=0.10)")
    ax.plot(xs, [r["lr_soft"]["video"] for r in pk], lw=1.2, ls="--", label="soft LR video")
    ax.plot(xs, [r["lr_decorr"]["video"] for r in pk], lw=1.2, ls=":", label="decorr LR video")
    ax.set_title("Peaked α: GLS concentrates as ρ↑")
    ax.set_xlabel(r"pairwise $\cos(g_m, g_{m'})$")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    mx = [r["rho"] for r in mse_rows]
    ax.plot(mx, [r["cos_equal"] for r in mse_rows], lw=2.0, label="equal-weight cos")
    ax.plot(mx, [r["cos_stat"] for r in mse_rows], lw=2.0, label="stat π* cos")
    ax.plot(mx, [r["cos_gls"] for r in mse_rows], lw=1.4, ls="--", label="pure GLS cos")
    ax.set_title("Factor model A: fused direction vs true μ (cosine)")
    ax.set_xlabel(r"noise correlation $\rho$")
    ax.set_ylabel("cosine")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    if replay:
        keys = sorted(replay)
        labels = [k.replace(":", "\n") for k in keys]
        x = np.arange(len(keys))
        w = 0.35
        ne = [replay[k]["mean_n_eff"] for k in keys]
        et = [replay[k]["mean_eta_global"] for k in keys]
        ax.bar(x - w / 2, ne, w, label=r"mean $N_{\mathrm{eff}}$")
        ax.bar(x + w / 2, et, w, label=r"mean $\eta_{t+1}$")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_title("Replay Amazon / MSR-VTT windows")
        ax.legend(frameon=False, fontsize=8)
        ax.axhline(1.0, color="#999", ls="--", lw=0.7)
    else:
        ax.text(0.5, 0.5, "no gradcos trajectory to replay", ha="center", va="center")
        ax.set_axis_off()

    fig.text(
        0.5,
        0.015,
        "pi* = (1-gamma) alpha + gamma GLS   gamma = 1 - N_eff/|M|   "
        "eta next ~ 1/sqrt(pi' R pi)   FWD always on",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _md_table(headers, rows):
    head = "| " + " | ".join(headers) + " |"
    sep = "|" + "|".join("---" if h != headers[0] else "---" for h in headers) + "|"
    # numeric right-align: keep simple
    sep = "|" + "|".join("---:" if i else "---" for i, _ in enumerate(headers)) + "|"
    body = "\n".join("| " + " | ".join(r) + " |" for r in rows)
    return "\n".join([head, sep, body])


def write_docs(sweep, mse_rows, replay, path: Path):
    eq = [r for r in sweep if r["alpha_kind"] == "equal"]
    pk = [r for r in sweep if r["alpha_kind"] == "peaked"]
    pick = [0, 6, 12, 18, 24]
    pick = [i for i in pick if i < len(eq)]

    def _r4(x):
        return f"{x:.2f}"

    sweep_rows = []
    for i in pick:
        r = eq[i]
        sweep_rows.append(
            [
                f"{r['rho']:+.2f}",
                _r4(r["n_eff"]),
                _r4(r["eta_global"]),
                _r4(r["ambiguity"]),
                _r4(r["lr_ratio_stat"]),
                "Y" if r["decorr_active"] else "n",
            ]
        )
    peak_rows = []
    for i in pick:
        r = pk[i]
        peak_rows.append(
            [
                f"{r['rho']:+.2f}",
                _r4(r["pi"]["video"]),
                _r4(r["pi"]["audio"]),
                _r4(r["gamma_gls"]),
                _r4(r["lr_stat"]["video"]),
                _r4(r["lr_stat"]["audio"]),
                _r4(r["lr_soft"]["video"]),
                r["roles"]["video"][0] + "/" + r["roles"]["audio"][0],
            ]
        )
    mse_pick = mse_rows[::2] or mse_rows
    mse_tbl = []
    for r in mse_pick:
        mse_tbl.append(
            [
                _r4(r["rho"]),
                _r4(r["n_eff"]),
                f"{r['cos_equal']:.3f}",
                f"{r['cos_stat']:.3f}",
                f"{r['cos_gls'] - r['cos_equal']:+.3f}",
                f"{r['mse_equal']:.3f}",
                f"{r['mse_blue']:.3f}",
            ]
        )
    replay_tbl = []
    for k, s in sorted(replay.items()):
        z = s.get("fisher_z_pooled") or {}
        replay_tbl.append(
            [
                k,
                f"{s['mean_rho']:+.3f}",
                _r4(s["mean_n_eff"]),
                _r4(s["mean_eta_global"]),
                _r4(s.get("mean_gamma", float("nan"))),
                _r4(s["mean_ambiguity"]),
                f"{z.get('p_value', float('nan')):.3f}",
                "Y" if z.get("significant_pos") else "n",
            ]
        )

    md = f"""# 模态梯度相关如何指导下一步学习率

同一件几何事实 `R_ij = cos(g_i, g_j)`，用 **ensemble learning** 和 **统计方法** 各说一遍，然后合成下一步的 `LR_m`.

## 设定

窗口 `t` 上，模态 `m ∈ M` 给出（common-dim）梯度签名 `g_m`。下一步融合更新是加权 ensemble：

```
Δθ_{{t+1}} = −η_{{t+1}}  Σ_m π_m g_m ,     π ≥ 0,  1ᵀπ = 1
```

`α_m` 是 attribution 路由质量（谁该动），`R` 是梯度相关（他们的票相不相干）。下一步 LR 必须同时吃这两口：

```
π_gls ∝ relu(R̃^{{-1}} α)
γ     = clip(1 − N_eff/|M|, 0, 1)      # independent → 0, collinear → 1
π*    = (1−γ) α + γ π_gls              # empirical-Bayes stacking mix
LR_m  = η_{{t+1}} · (β + (1−β) |M| π*_m) · c_m
η_{{t+1}} = √( (1/|M|) / (π*ᵀ R̃ π*) ) · φ_t
```

`R̃` 是 PSD + equicorrelation shrinkage 之后的 Gram；`φ_t = temporal_gain(cos(g_t, g_{{t-1}}))`；`c_m` 是 `cos(g_m, g_shared)<0` 时的 conflict damp。FWD 始终开，只改下一窗的 adapt 步长。

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
N_eff^{{Kish}} = |M| / (1 + (|M|−1) ρ̄) = |M|² / (1ᵀ R 1)
N_eff^{{GLS}}  = 1ᵀ R^{{-1}} 1     (≥ Kish；相关非均匀时 GLS 更有效)
```

**指导下一步：** `η_{{t+1}}` 按 `√(N_eff / |M|)` 收缩，使融合步的方差稳住在独立等权基线 `Var = 1/|M|`。这就是 variance-stabilising step。

### 2. 多样性–精度分解（连续版 leader / diversifier / redundant）

| 角色 | 离散（`soft_decorr`） | 连续（本文 `π*`） |
|---|---|---|
| leader | unique mass 最大，占住共享方向 | `π*_m` 大：高精度票 |
| diversifier | 相对 leader 的残差 `1−cos₊` 高 | partial uniqueness `1/P_{{mm}}` 高 |
| redundant | 与 leader 共线 | `π*` 被 clip 到 0 / 地板 |

等权 + 等相关时，**统计量只缩全局 `η`、不改分配**（选民可交换）；**离散角色会任意点一个 leader**。`α` 一旦不均，两种观点合流：高相关把质量集中到 `argmax α`。

### 3. Stacking / 最小方差组合

线性 ensemble 的最优权（非负约束前）就是最小方差投资组合 / stacking：

```
min_π  πᵀ R π    s.t. πᵀ α 最大化   ⇒   π* ∝ R^{{-1}} α
```

这是 **SNR 匹配滤波器**：`(πᵀα)² / (πᵀ R π)`。Negative correlation learning 的合法部分落在 `R_ij<0` 且仍与 `g_shared` 对齐——那是真多样性。`R_ij<0` 且 `cos(g_m, g_shared)<0` 是 **冲突**，不是多样性：conflict damp，而不是加权。

---

## 统计视角：相关分数、BLUE、收缩、检验

### 1. 相关 Gram 与收缩

`R_ij = cos(g_i, g_j)`（签名已按 common-dim 对齐；中心化后即 Pearson）。`|M|∈{{2,3}}` 时样本 Gram 极噪：

```
R_psd = Higham(R)                 # clip 负特征值（冲突子空间）
R̃    = (1−λ) R_psd + λ T         # T = equicorrelation((1−ρ̄)I + ρ̄ 11ᵀ)
λ     = clip( log(cond) / log(κ*), 0, λ_max )
```

Equicorrelation 先验 = 单因子模型（模态共享一个 common descent factor）。`λ` 随条件数升高——这就是 Ledoit–Wolf / Schäfer–Strimmer 在小 `M` 上的可用形式。

### 2. Gauss–Markov / GLS

观察 `g_m = α_m μ + ε_m`, `Cov(ε)=σ² R`。`μ` 的 BLUE / 匹配滤波权就是 `R^{{-1}} α`。`soft_gls` 直接用这组权；`soft_stat` 再用 `γ=1−N_eff/|M|` 把它们往独立选民先验 `α` 上收缩，并乘方差稳定 `η` 和时序/冲突增益。`|M|∈{{2,3}}` 时纯 GLS 会过早把弱模态 clip 到 0，混合是必要的。

`R̃ ≈ I` 时 `γ≈0`、`π* ≈ α`、`η ≈ 1` → **退回 soft LR**。`R̃ ≈ 11ᵀ` 时 `N_eff → 1`、`γ→1`、`η → 1/√|M|`，且 `π*` 集中到高于平均的 `α`。

### 3. 偏相关 / 精度矩阵

`P = R̃^{{-1}}`。模态 `m` 对他人回归后的残差方差是 `1/P_{{mm}}`（uniqueness）。偏相关

```
ρ_{{ij|rest}} = −P_{{ij}} / √(P_{{ii}} P_{{jj}})
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

{_md_table(["ρ", "N_eff", "η", "ambiguity", "LR ratio", "decorr"], sweep_rows)}

## 合成扫描（peaked α = 0.70 / 0.20 / 0.10）

{_md_table(["ρ", "π video", "π audio", "γ", "LR_stat v", "LR_stat a", "LR_soft v", "roles v/a"], peak_rows)}

ρ 低：γ≈0，`π ≈ α`，stat ≈ soft（独立选民）。ρ 高：γ↑，`π` 向 GLS 集中到 video，audio 的下一步 LR 被统计权+全局 `η` 同时压下。

## 因子模型

- **A**（匹配滤波）：`g_m = α_m μ + ε_m`，指标是融合方向与真 μ 的 **cosine**（simplex 权不解偏尺度，所以不用生 MSE）。
- **B**（异方差投票）：`g_m = μ + ε_m`, `Var(ε_m)∝1/α_m`，指标是 BLUE 组合相对 μ 的 MSE。

{_md_table(["ρ", "N_eff", "cos equal", "cos stat", "Δcos GLS−eq", "MSE eq (B)", "MSE BLUE (B)"], mse_tbl)}

低相关时 cosine 上 GLS/stat 优于等权。ρ 升高时等权方向因多样性崩塌而变差，GLS 仍盯住高 α 票，Δcos 拉大。Model B 的 BLUE 方差始终低于等权（异方差选民的 Gauss–Markov）。

## Amazon / MSR-VTT 轨迹 replay

用已有 grad-cos 窗口的 `pair_cos / α / temporal_cos` 反推 `N_eff, η, π*`（不重训，刻画 *guidance* 本身）。

{_md_table(["traj", "ρ̄", "N_eff", "η", "γ", "ambiguity", "Fisher p", "ρ>0 sig"], replay_tbl) if replay_tbl else "_no gradcos json_"}

经验对照：Amazon 窗间 `ρ̄` 高、`N_eff≈1.2` → 统计量会 **降全局步长并（在 γ 升高后）集中 π***；MSR-VTT `ρ̄` 低、`N_eff≈2.2`（接近 `|M|=3`）→ 接近独立选民，γ 小，soft 即可，强行 GLS 集中没有统计理由。轨迹级 Fisher z 因 T=6 检验力不足，应读 `N_eff` / `ρ̄` 而不是 p 值。

## 调度器

| name | 做什么 |
|---|---|
| `equal` | 稠密等 LR |
| `soft` | `α → LR_m` |
| `soft_decorr` | 离散角色增益（高相关才激活） |
| `soft_gls` | 只用 `π* ∝ R̃^{{-1}} α` |
| `soft_stat` | `π* = (1-γ)α + γ GLS` × 方差稳定 `η` × 时序 × 冲突 |

```bash
PYTHONPATH=. python3 scripts/run_agod_grad_corr_stat.py
PYTHONPATH=. python3 -m tests.test_agod_grad_corr_stat
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(sweep, mse_rows, replay, path: Path):
    eq = [r for r in sweep if r["alpha_kind"] == "equal"]
    pick = [i for i in (0, 6, 12, 18, 24) if i < len(eq)]
    lines = []
    for i in pick:
        r = eq[i]
        lines.append(
            f"{r['rho']:+.2f} & {r['n_eff']:.2f} & {r['eta_global']:.2f} & "
            f"{r['ambiguity']:.2f} & {r['lr_ratio_stat']:.2f} \\\\"
        )
    replay_lines = []
    for k, s in sorted(replay.items()):
        replay_lines.append(
            f"{k.replace('_', '\\_')} & {s['mean_rho']:+.3f} & {s['mean_n_eff']:.2f} & "
            f"{s['mean_eta_global']:.2f} & {s['mean_ambiguity']:.2f} \\\\"
        )
    tex = (
        "% Modality-grad correlation → next LR\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Effective ensemble size and variance-stabilising next-step "
        "scale as pairwise gradient cosine varies (equal $\\alpha$).}\n"
        "\\label{tab:agod-grad-corr-stat}\n"
        "\\begin{tabular}{rrrrr}\\toprule\n"
        "$\\rho$ & $N_{\\mathrm{eff}}$ & $\\eta_{t+1}$ & ambiguity & LR ratio \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    if replay_lines:
        tex += (
            "\n\\begin{table}[t]\\centering\n"
            "\\caption{Replay of Amazon / MSR-VTT grad-cos windows: "
            "correlation-guided next LR diagnostics.}\n"
            "\\label{tab:agod-grad-corr-replay}\n"
            "\\begin{tabular}{lrrrr}\\toprule\n"
            "trajectory & $\\bar\\rho$ & $N_{\\mathrm{eff}}$ & $\\eta$ & ambiguity \\\\\n"
            "\\midrule\n"
            + "\n".join(replay_lines)
            + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    print("isotropic sweep…", flush=True)
    sweep = isotropic_sweep()
    print("factor-model MSE…", flush=True)
    mse_rows = factor_mse_sweep()
    print("replay gradcos trajectories…", flush=True)
    replay = replay_gradcos(GRADCOS)

    payload = {
        "agod_version": "0.1.0",
        "focus": "modality gradient correlation → next LR (ensemble + statistics)",
        "control_law": {
            "pi": "relu(R_shrink^{-1} alpha) / 1^T relu(.)",
            "eta": "sqrt((1/M) / (pi^T R pi)) * temporal_gain(phi)",
            "LR_m": "eta * (beta + (1-beta) |M| pi_m) * conflict_gain_m",
        },
        "sweep": sweep,
        "factor_mse": mse_rows,
        "replay": {
            k: {kk: vv for kk, vv in s.items() if kk != "windows"}
            | {"n_windows_logged": len(s.get("windows") or [])}
            for k, s in replay.items()
        },
        "replay_windows": {k: s.get("windows") for k, s in replay.items()},
    }
    (OUT / "agod_grad_corr_stat.json").write_text(json.dumps(payload, indent=2))
    plot_board(sweep, mse_rows, replay, OUT / "AGOD_Grad_Corr_Stat_Board.png")
    write_docs(sweep, mse_rows, replay, OUT / "README.md")
    write_latex(sweep, mse_rows, replay, OUT / "AGOD_grad_corr_stat_tables_only.tex")
    write_docs(sweep, mse_rows, replay, DOCS / "AGOD_grad_corr_stat.md")
    write_latex(sweep, mse_rows, replay, DOCS / "AGOD_grad_corr_stat_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Grad_Corr_Stat_Board.png", ART / "AGOD_Grad_Corr_Stat_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")

    print("\n=== equal-α sweep (subset) ===", flush=True)
    for r in sweep:
        if r["alpha_kind"] != "equal":
            continue
        if abs(r["rho"] - round(r["rho"] * 4) / 4) > 0.03:
            continue
        print(
            f"  ρ={r['rho']:+.2f}  N_eff={r['n_eff']:.2f}  η={r['eta_global']:.2f}  "
            f"amb={r['ambiguity']:.2f}  LRratio={r['lr_ratio_stat']:.2f}",
            flush=True,
        )
    print("\n=== replay ===", flush=True)
    for k, s in sorted(replay.items()):
        print(
            f"  {k:24s} ρ̄={s['mean_rho']:+.3f}  N_eff={s['mean_n_eff']:.2f}  "
            f"η={s['mean_eta_global']:.2f}  amb={s['mean_ambiguity']:.2f}",
            flush=True,
        )
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
