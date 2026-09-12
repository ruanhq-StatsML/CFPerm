#!/usr/bin/env python3
"""Finite-step and other stacking logics — explicitly not MoE.

  PYTHONPATH=. python3 scripts/run_agod_stack_logics.py
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

from agod.stack_logics import (
    cauchy_step,
    entropy_regularized,
    eta_path,
    frank_wolfe_linear,
    mgda_min_norm,
    mixloss_weights,
    n_path,
    pseudo_bma_weights,
    quadratic_finite_step,
    rayleigh_gls,
    tau_path,
    taylor_remainder_scale,
)

OUT = ROOT / "results" / "agod_stack_logics"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_stack_logics")

EXP = ("a", "b", "c")
SCORES = {"a": 0.85, "b": 0.45, "c": 0.10}
LOSSES = {"a": 0.15, "b": 0.55, "c": 0.90}


def gram():
    r = np.eye(3)
    r[0, 1] = r[1, 0] = 0.80
    r[0, 2] = r[2, 0] = 0.10
    r[1, 2] = r[2, 1] = 0.10
    return r


def complementary():
    """Unit grads e0, e1, e2 matching (e0+e1)/√2 — linear gain is a vertex."""
    s = {"a": 0.5**0.5, "b": 0.5**0.5, "c": 0.0}
    r = np.eye(3)
    return {
        "scores": s,
        "vertex": quadratic_finite_step(s, r, EXP, eta=0.0),
        "train": quadratic_finite_step(s, r, EXP, eta=3e-3),
        "eta1": quadratic_finite_step(s, r, EXP, eta=1.0),
        "joint": rayleigh_gls(s, r, EXP)["pi"],
        "fw": frank_wolfe_linear(s, EXP),
    }


def plot_board(etas, taus, snapshots, taylor, path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.2), facecolor="#f7f5f1")
    fig.suptitle(
        "Not MoE: finite-step / entropy / MGDA stacking logics",
        fontsize=13,
        fontweight="bold",
    )
    xs = [r["eta"] for r in etas]
    ax = axes[0, 0]
    for e in EXP:
        ax.semilogx(xs, [r["pi"][e] for r in etas], lw=2.0, label=f"π_{e}")
    ax.set_title("Finite-step η-path  (vertex → min-var)")
    ax.set_xlabel("virtual step η")
    ax.legend(frameon=False, fontsize=8)
    ax.axvline(3e-3, color="#999", ls="--", lw=0.8)
    ax.axvline(snapshots["joint_eta"], color="#c45", ls=":", lw=1.0)

    xs = [r["tau"] for r in taus]
    ax = axes[0, 1]
    for e in EXP:
        ax.semilogx(xs, [r["pi"][e] for r in taus], lw=2.0, label=f"π_{e}")
    ax.set_title("Entropy τ-path  (vertex → equal, ignores R)")
    ax.set_xlabel("temperature τ")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    names = ["vertex", "train η", "η=1", "joint GLS", "entropy", "MGDA", "BMA n=50"]
    keys = ["vertex", "train", "eta1", "joint", "ent", "mgda", "bma50"]
    x = np.arange(len(names))
    w = 0.22
    for i, e in enumerate(EXP):
        ax.bar(x + (i - 1) * w, [snapshots[k][e] for k in keys], w, label=e)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=22, ha="right", fontsize=8)
    ax.set_title("Same (s, R): sequential η vs joint vs other logics")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    et = [t["eta"] for t in taylor]
    ax.loglog(et, [t["linear"] for t in taylor], lw=2.0, label="linear ~ η")
    ax.loglog(et, [t["quadratic"] for t in taylor], lw=2.0, label="quadratic ~ η²")
    ax.loglog(et, [t["cubic"] for t in taylor], lw=1.4, ls="--", label="O(η³)")
    ax.axvline(3e-3, color="#999", ls="--", lw=0.8)
    ax.set_title("Taylor: adapt LR is linear-regime")
    ax.set_xlabel("η")
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.015,
        "π is global on votes, not π(x)  ·  joint (η,π) = GLS  ·  train LR ≈ vertex  ·  BMA ≠ stacking",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(snapshots, etas, taylor_lr, ns, comp, path: Path):
    def fmt(pi):
        return ", ".join(f"{e}={pi[e]:.2f}" for e in EXP)

    md = f"""# 不是 MoE：有限步长和其他 stacking 逻辑

结论先说：**这不是 Mixture of Experts。** MoE（Jacobs et al. 1991；Switch Transformer 也是）是 **输入条件门** `π(x)=softmax(g(x))`，expert 按 x 的区域分工、和门一起训练。这里的 π 是 **全局（或慢变）的票组合权**，用 holdout / one-step-ahead 损失更新——Super Learner / forecast combination，门不看 x。

`softmax(s/τ)` 长得像 MoE 的门，但它的 logit 是 **expert 分数 s**，不是 `Wx`。没有 x，就不是 MoE。

---

## 有限步长（线性 gain 的下一项）

沿堆叠方向走真实一步：

```
L(θ − η Gπ) = L − η πᵀs + (η²/2) πᵀ H π + O(η³)
```

Gauss–Newton 下 `H ≈ GᵀG = R`。丢掉余项，单纯形上的 stacking 是

```
min_{{π∈Δ}}  −πᵀs + (η/2) πᵀ R π
```

| η | 行为 |
|---|---|
| → 0 | 线性项主导 → **顶点**（离散 OSL / Frank–Wolfe） |
| 中等 | 精度–方差权衡，内点组合 |
| → ∞ | 二次项主导 → **min-var** `min πᵀRπ`（趋近 MGDA） |
| 训练 LR ~ 3×10⁻³ | `quad/lin = η/2 ≪ 1`，真实 SGD 步几乎就是线性 gain |

所以：**不能指望实际学习率带来的曲率帮你组合。** 虚步长 η 是 meta-loss 上「敢不敢离开 winner-take-all」的旋钮，和训练 LR 不是同一个东西。

### 顺序 vs 联合：有限步长 *就是* GLS

把 (η, π) 一起放进二阶 Taylor：

1. **顺序**：先钉死训练 LR `η≈3e-3`，再选 π → 线性项主导 → **顶点**。
2. **线搜索（π 固定）**：Cauchy 步 `η* = (πᵀs)/(πᵀRπ)`。
3. **联合**：把 η* 代回去得到 `min_π −(πᵀs)² / (2 πᵀRπ)`，即最大化 Rayleigh / SNR。无约束解 `π ∝ R⁻¹ s`，就是方向匹配 / 静态 GLS。

有限步长不是另一种估计器，它是 **离开顶点的机制**。联合优化步长和权重，回到已经写过的 GLS；只用训练 LR 走一步，永远走不出离散 OSL。

Taylor 在 η=3e-3：linear={taylor_lr['linear']:.2e}，quad/lin={taylor_lr['quad_over_lin']:.2e}（线性区）；η=1 才离开线性区。本例联合 Cauchy `η*={snapshots['joint_eta']:.2f}`（图中红点线）。

η-path（s_a=0.85, s_b=0.45, s_c=0.10，a–b 相关 0.8）：

| η | π | N_eff | top |
|---|---|---:|---|
""" + "\n".join(
        f"| {r['eta']:.3g} | {fmt(r['pi'])} | {r['n_eff']:.2f} | {r['top']} |"
        for r in etas[::3]
    ) + f"""

互补单位梯度（holdout ∥ (e₀+e₁)/√2，R=I）：线性 / 训练 LR / Frank–Wolfe 都塌到一个轴；η=1 和联合 GLS 才是 ½–½。

| logic | π |
|---|---|
| vertex / train / FW | {fmt(comp['vertex'])} |
| finite-step η=1 | {fmt(comp['eta1'])} |
| joint GLS | {fmt(comp['joint'])} |

---

## 其他逻辑（同样不是 MoE）

| 逻辑 | 公式 | 用 holdout 靶？ | 用 R？ | π 的路径 |
|---|---|---|---|---|
| 线性 gain | `max πᵀs` | 是（s） | 否 | 只有顶点 |
| Frank–Wolfe | 线性目标的 FW 步 | 是 | 否 | 每步仍是顶点 |
| **有限步长二次（顺序）** | `−πᵀs+(η/2)πᵀRπ`，η 给定 | 是 | 是 | 顶点 → min-var |
| **联合 (η,π)** | `max (πᵀs)²/(πᵀRπ)` | 是 | 是 | = GLS / 方向匹配 |
| 方向匹配 | `‖Gπ−ĝ_hold‖²` | 是 | 是 | 同上，η 无关的内点 |
| **熵正则** | `softmax(s/τ)` | 是 | **否** | 顶点 → **equal**（不拆 clone） |
| **MGDA** | `min ‖Gπ‖²` | **否** | 是 | 冲突几何的 min-norm，不是 stacking |
| mixloss / 聚合 | `−log Σ π_m e^{{-L_m}}` | 损失 | 否 | Vovk aggregating；proper scoring of the *mixture* |
| pseudo-BMA | `π ∝ e^{{-n L}}` | 边缘似然/损失 | 否 | `n=1` 温和；`n=|window|` 在 M-open 塌缩 |
| Hedge+share | 乘性 × 均匀 | 是 | 否 | tracking，不是门 |
| GLS-EWMA | `Σ̂^{{-1}}1` | 误差 | 是 | 在线 Bates–Granger |

同一组 (s, R) 上的快照：

| logic | π |
|---|---|
| vertex η=0 | {fmt(snapshots['vertex'])} |
| train η=3e-3 | {fmt(snapshots['train'])} |
| finite-step η=1 | {fmt(snapshots['eta1'])} |
| finite-step η=8 | {fmt(snapshots['eta8'])} |
| joint GLS (η*={snapshots['joint_eta']:.2f}) | {fmt(snapshots['joint'])} |
| entropy τ=1 | {fmt(snapshots['ent'])} |
| MGDA（无靶） | {fmt(snapshots['mgda'])} |
| pseudo-BMA n=1 | {fmt(snapshots['bma1'])} |
| pseudo-BMA n=50 | {fmt(snapshots['bma50'])} |

读图：

- **训练 LR ≈ 顶点**：`train` 和 `vertex` 几乎一样。真实一步不会组合。
- **η=8** 把质量从共线的 a/b 分给独立的 c（二次项惩罚相关）。
- **联合 GLS** 在虚步长 η* 处取 Rayleigh 最优，不是「再把训练 LR 调大一点」。
- **熵正则不管 R**：τ 变大只是摊成 1/3，clone trap 解不了。
- **MGDA 没有 s**：更像「别打架」而不是「跟 holdout 对齐」。Sener & Koltun NeurIPS 2018 是多任务 Pareto，不是 stacking。
- **pseudo-BMA**：n=1 时损失差 0.4 nats 几乎还在混；n=50（ELPD 尺度）塌向 a。Yao, Vehtari, Simpson, Gelman (*Bayesian Analysis* 2018) 才主张 stacking 而不是 BMA——因为 stacking 优化的是预测风险，softmax **不乘 n**。

Armijo / 线搜索是有限步长的自适应 η：在 holdout 上沿 `d=Gπ` 找使 `L(θ−η d)` 下降的 η。训练 LR 太小，线搜索给出的虚步长才会进入组合区；把这条虚步长代回 π 的二次型，就是联合 GLS。

### 还相邻、但不要混进来的

| 名字 | 为什么不是这套 stacking |
|---|---|
| MoE / Switch Transformer | `π(x)=softmax(Wx)`，门看样本 |
| Gradient Blending (Wang et al. CVPR 2020) | 用过拟合间隙调模态权重，不是 OOF 组合 |
| PCGrad / OGM-GE | 梯度投影 / 冲突手术，不在单纯形上估 π |
| NCL (Liu & Yao) | 训练专家时罚相关，改的是专家本身不是 meta π |
| Mixup / 输入混合 | 混合的是 x，不是专家票 |

---

## 和 MoE 的对照（避免混）

| | MoE | 这里的 online stacking |
|---|---|---|
| π | `π_m(x)`，看样本 | `π_m(t)`，看 expert 票 / 损失 |
| 训练 | 门和 expert 联合（EM 或联合梯度） | expert 当黑盒，meta 只用 OOF |
| 目标 | 混合似然 `log Σ π_m(x) p_m(y\\|x)` | 预测 / 方向风险 `L(π·v, y)` 或 `‖Gπ−g_hold‖²` |
| 多样性 | specialization by region of x | error–ambiguity / `N_eff` of votes |
| 文献 | Jacobs 1991；Switch Transformer | Wolpert；Super Learner；OSL；Bates–Granger |

熵正则的 softmax **不是** MoE 门，只是把分数温度化。若把 s 换成 `Wx`，那才滑向 MoE——不要滑。

```bash
PYTHONPATH=. python3 scripts/run_agod_stack_logics.py
PYTHONPATH=. python3 -m tests.test_agod_stack_logics
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(snapshots, path: Path):
    lines = []
    labels = [
        ("vertex", "vertex $\\eta=0$"),
        ("train", "train $\\eta=3\\times10^{-3}$"),
        ("eta1", "finite-step $\\eta=1$"),
        ("eta8", "finite-step $\\eta=8$"),
        ("joint", "joint GLS"),
        ("ent", "entropy $\\tau=1$"),
        ("mgda", "MGDA (no target)"),
        ("bma1", "pseudo-BMA $n=1$"),
        ("bma50", "pseudo-BMA $n=50$"),
    ]
    for k, lab in labels:
        p = snapshots[k]
        lines.append(f"{lab} & {p['a']:.2f} & {p['b']:.2f} & {p['c']:.2f} \\\\")
    tex = (
        "% Finite-step and other stacking logics (not MoE)\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Same scores and Gram, different combination logics. "
        "Not an input-conditional MoE gate. Sequential tiny $\\eta$ stays "
        "at a vertex; joint $(\\eta,\\pi)$ recovers GLS.}\n"
        "\\label{tab:agod-stack-logics}\n"
        "\\begin{tabular}{lccc}\\toprule\n"
        "logic & $\\pi_a$ & $\\pi_b$ & $\\pi_c$ \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    r = gram()
    etas = eta_path(SCORES, r, EXP, np.logspace(-2, 1.1, 16))
    taus = tau_path(SCORES, EXP, np.logspace(-2, 1.1, 16))
    ns = n_path(LOSSES, EXP, np.logspace(0, 2, 12))
    joint = rayleigh_gls(SCORES, r, EXP)
    equal = {e: 1.0 / 3.0 for e in EXP}
    snapshots = {
        "vertex": quadratic_finite_step(SCORES, r, EXP, eta=0.0),
        "train": quadratic_finite_step(SCORES, r, EXP, eta=3e-3),
        "eta1": quadratic_finite_step(SCORES, r, EXP, eta=1.0),
        "eta8": quadratic_finite_step(SCORES, r, EXP, eta=8.0),
        "joint": joint["pi"],
        "joint_eta": float(joint["eta_star"]),
        "cauchy_equal": cauchy_step(equal, SCORES, r, EXP),
        "fw": frank_wolfe_linear(SCORES, EXP),
        "ent": entropy_regularized(SCORES, EXP, tau=1.0),
        "mgda": mgda_min_norm(r, EXP),
        "bma1": pseudo_bma_weights(LOSSES, EXP, n_eff=1.0),
        "bma50": pseudo_bma_weights(LOSSES, EXP, n_eff=50.0),
        "mix": mixloss_weights(LOSSES, EXP, eta=2.0),
    }
    taylor = [taylor_remainder_scale(float(e)) for e in np.logspace(-3.5, 0.5, 18)]
    taylor_lr = taylor_remainder_scale(3e-3)
    comp = complementary()

    payload = {
        "not_moe": True,
        "scores": SCORES,
        "snapshots": snapshots,
        "eta_path": etas,
        "tau_path": taus,
        "n_path": ns,
        "taylor_lr": taylor_lr,
        "complementary": comp,
        "joint_snr": joint["snr"],
    }
    (OUT / "agod_stack_logics.json").write_text(json.dumps(_jsonable(payload), indent=2))
    plot_board(etas, taus, snapshots, taylor, OUT / "AGOD_Stack_Logics_Board.png")
    write_docs(snapshots, etas, taylor_lr, ns, comp, OUT / "README.md")
    write_latex(snapshots, OUT / "AGOD_stack_logics_tables_only.tex")
    write_docs(snapshots, etas, taylor_lr, ns, comp, DOCS / "AGOD_stack_logics.md")
    write_latex(snapshots, DOCS / "AGOD_stack_logics_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Stack_Logics_Board.png", ART / "AGOD_Stack_Logics_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    shutil.copy2(
        OUT / "AGOD_Stack_Logics_Board.png",
        Path("/opt/cursor/artifacts/agod_stack_logics_board.png"),
    )
    pi_keys = ("vertex", "train", "eta1", "joint", "ent", "mgda", "bma1", "bma50")
    shown = {k: {e: round(snapshots[k][e], 2) for e in EXP} for k in pi_keys}
    print("snapshots:", shown, flush=True)
    print(f"joint η*={snapshots['joint_eta']:.3f}  snr={joint['snr']:.3f}", flush=True)
    print(
        f"taylor η=3e-3 linear_regime={taylor_lr['linear_regime']} "
        f"quad/lin={taylor_lr['quad_over_lin']:.3e}",
        flush=True,
    )
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
