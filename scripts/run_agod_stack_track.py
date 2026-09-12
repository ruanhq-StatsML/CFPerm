#!/usr/bin/env python3
"""Half-stream switch, OOF features, and π variation.

  PYTHONPATH=. python3 scripts/run_agod_stack_track.py
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

from agod.stack_track import (
    alignment_leak,
    oof_vs_leaky_stack,
    path_tv,
    run_switch_methods,
    sticky_path,
)

OUT = ROOT / "results" / "agod_stack_track"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_stack_track")
SEED = 2026


def stream_switch(t: int = 80, seed: int = SEED):
    rng = np.random.default_rng(seed)
    y = np.zeros(t)
    mid = t // 2
    rows = []
    for i in range(t):
        if i < mid:
            a, b = rng.normal(0, 0.12), rng.normal(0, 1.05)
        else:
            a, b = rng.normal(0, 1.05), rng.normal(0, 0.12)
        rows.append((a, b, rng.normal(0, 0.65)))
    return np.array(rows), y, ("a", "b", "c"), mid


def stream_stationary(t: int = 80, seed: int = SEED + 3):
    rng = np.random.default_rng(seed)
    y = np.zeros(t)
    rows = np.column_stack(
        [rng.normal(0, 0.25, t), rng.normal(0, 0.28, t), rng.normal(0, 0.30, t)]
    )
    return rows, y, ("a", "b", "c")


def oof_demo(seed: int = 6) -> dict:
    rng = np.random.default_rng(seed)
    n, d = 120, 60
    xs = rng.normal(size=(n, 1))
    y = xs[:, 0] + 0.18 * rng.normal(size=n)
    xn = rng.normal(size=(n, d))
    return oof_vs_leaky_stack(xs, xn, y, n_folds=4, lam=1e-3)


def grad_leak_demo(seed: int = SEED + 9) -> dict:
    rng = np.random.default_rng(seed)
    g = rng.normal(size=64)
    return alignment_leak(
        g + 0.55 * rng.normal(size=64),
        g + 0.55 * rng.normal(size=64),
        g,
    )


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if obj is None:
        return None
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def plot_board(switch, packed, sticky, oof, path: Path):
    names = ("a", "b", "c")
    mid = switch["mid"]
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.2), facecolor="#f7f5f1")
    fig.suptitle(
        "Half-stream switch / OOF features / π variation",
        fontsize=13,
        fontweight="bold",
    )

    ax = axes[0, 0]
    for key, lab, ls in (
        ("hedge_share0", "Hedge share=0", "-"),
        ("hedge_share", "Hedge+share", "-"),
        ("osl_disc", "discrete OSL", "--"),
        ("osl_window", "window OSL", "-"),
    ):
        ax.plot(
            [p["b"] for p in packed[key]["traj"]],
            lw=2.0,
            ls=ls,
            label=lab,
        )
    ax.axvline(mid, color="#999", ls=":", lw=1.0)
    ax.set_title("π_b(t) after a→b switch")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    for key, lab in (
        ("hedge_share0", "share=0"),
        ("hedge_share", "share=0.08"),
        ("osl_disc", "osl_disc"),
        ("osl_window", "window"),
    ):
        traj = packed[key]["traj"]
        tvs = [0.0] + [
            0.5
            * sum(abs(traj[i][n] - traj[i + 1][n]) for n in names)
            for i in range(len(traj) - 1)
        ]
        ax.plot(tvs, lw=1.6, label=lab)
    ax.axvline(mid, color="#999", ls=":", lw=1.0)
    ax.set_title("TV(π_t, π_{t+1})  — spend variation at the switch")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    keys = ["hedge_share0", "hedge_share", "osl_disc", "osl_window", "osl_sgd"]
    labels = ["share=0", "share", "OSL disc", "window", "OSL SGD"]
    delays = [
        packed[k]["delay_from_switch"] if packed[k]["delay_from_switch"] is not None else 40
        for k in keys
    ]
    tvs = [packed[k]["path_tv"] for k in keys]
    tvs.append(path_tv(sticky))
    delays.append(
        next(
            (
                i - mid
                for i, p in enumerate(sticky)
                if i >= mid and p["b"] >= 0.45
            ),
            40,
        )
    )
    labels.append("sticky λ=0.12")
    ax.scatter(tvs, delays, s=48)
    for x, y, lab in zip(tvs, delays, labels):
        ax.annotate(lab, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax.set_xlabel("path TV")
    ax.set_ylabel("delay after switch")
    ax.set_title("Tracking delay vs π variation")

    ax = axes[1, 1]
    labs = ["in-sample", "holdout"]
    x = np.arange(len(labs))
    ax.bar(x - 0.15, [oof["in_sample_oof"], oof["holdout_oof"]], 0.3, label="OOF π")
    ax.bar(x + 0.15, [oof["in_sample_leaky"], oof["holdout_leaky"]], 0.3, label="leaky π")
    ax.set_xticks(x)
    ax.set_xticklabels(labs)
    ax.set_title("OOF vs leaky level-1 features")
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.015,
        "share=0 / discrete OSL lock on first-half winner  ·  OOF votes only  ·  sticky π for the actuator",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(packed, oof, leak, stat_tv, switch_tv, path: Path):
    def fmt_pi(pi):
        return ", ".join(f"{k}={pi[k]:.2f}" for k in pi)

    rows = []
    for key in ("hedge_share0", "hedge_share", "osl_disc", "osl_window", "osl_sgd"):
        r = packed[key]
        d = r["delay_from_switch"]
        rows.append(
            f"| `{key}` | {r['preq_mean']:.3f} | {r['path_tv']:.2f} | "
            f"{r['tv_after_frac']:.2f} | {d if d is not None else '∞'} | "
            f"{r.get('final_leader', '?')} | {fmt_pi(r['final_pi'])} |"
        )

    md = f"""# 半程换人、OOF 特征、π 的变动

上一张讲的是 **选哪一种 meta-loss**（线性顶点 / 有限步长 / GLS）。这一张是 **协议里还没讲透的三件事**：最好的 expert 会换人、meta 只能看 OOF 特征、π 在单纯形上怎么走。还是不是 MoE：门不看 x。

---

## 半程换人（tracking，不是静态遗憾）

最好 expert 在 `t=T/2` 从 a 换成 b。对 **全程最好的固定 expert** 的遗憾是错的 oracle——那个固定赢家往往是一直中等的 c。对的 oracle 是 **最多换 k 次的比较器**（Herbster & Warmuth 1998 switching regret）。

三种锁死 / 三种花钱方式：

| 机制 | 会怎样 |
|---|---|
| 离散 OSL（累加损失） | 选出全程最好的**固定**专家（本例是中等的 c），看不见 a→b |
| Hedge `share=0` | 乘性更新从 ~0 质量爬升是对数慢的，半程经常来不及 |
| Hedge **Fixed-Share** `π ← (1-α)π⊙e^{{-ηL}} + α/\\|E\\|` | 每个死专家每步至少 `α/\\|E\\|`，能复活 |
| 滑窗 OSL（只看最近 W 步） | 忘掉前半程，W 步后跳到 b |
| 睡眠专家 / specialists | 前半程没投票的人不当成「很差」，只当没出场 |

share 不是学习率。它是 **允许 π 花掉的复活预算**：α=0 跟踪不动；α 太大则全程抖动，定常流上也在烧 TV。

本例 `T=80`、a→b 在 t=40：

| method | preq | path TV | TV after switch | delay | leader | final π |
|---|---:|---:|---:|---:|---|---|
""" + "\n".join(rows) + f"""

读表：离散 OSL 的累加损失选出的是全程最好的**固定**专家 c（一直中等），π 停在 c 上，a→b 的切换它看不见。Fixed-Share delay=5 把质量搬到 b，preq 从 0.40 掉到 0.06。`share=0` 最后也到了 b，但晚了 35 步，整段风险已经付过了。滑窗 OSL 是一次迟到的跳跃（delay=10）。Hedge `share=0` 从 ~0 质量爬升是对数慢的。

**负 regret vs best fixed expert 是预期，不是 bug。** 跟踪赢的是换人的比较器；离散 OSL 对齐的是那个「一直中等」的固定赢家。

---

## OOF 特征（Wolpert 的 level-1 矩阵）

Super Learner 的 meta 特征是 `Z_{{i,m}} = f_m^{{(-i)}}(x_i)`：第 m 个专家在 **没见过 i** 时对 i 的预测。用 `f_m(x_i)`（in-sample）当 stacking 特征，就是 Wolpert 警告过的 leaky stacking。流上的翻译：

```
P_t probe  → v_{{m,t}}     # level-1 票（OOF）
H_t holdout → y_t / g_hold # 靶，和 P_t 不相交
score with frozen π_t
train bases on A_t
update π_{{t+1}}
```

三种泄漏，严重程度不一样：

| 泄漏 | 发生了什么 | 看起来 |
|---|---|---|
| in-sample 票 | 专家在同一点上训过再投票 | meta 以为噪声专家也好 |
| 同一窗既当票又当靶 | `g` 同时是 vote 和 `g_hold` | 方向余弦被抬成 1 |
| 先 update 再 score | `leaky_step` | 低估 prequential risk |

本例：信号专家看 1 维 x，噪声专家看 60 维纯噪声（in-sample 能背下来）。

| | π_good | π_noise | corr(Z_noise, y) | in-sample MSE | holdout MSE | optimism |
|---|---:|---:|---:|---:|---:|---:|
| OOF | {oof['pi_oof']['good']:.2f} | {oof['pi_oof']['noise']:.2f} | {oof['corr_oof_noise']:.2f} | {oof['in_sample_oof']:.3f} | {oof['holdout_oof']:.3f} | {oof['optimism_oof']:.3f} |
| leaky | {oof['pi_leak']['good']:.2f} | {oof['pi_leak']['noise']:.2f} | {oof['corr_leak_noise']:.2f} | {oof['in_sample_leaky']:.3f} | {oof['holdout_leaky']:.3f} | {oof['optimism_leaky']:.3f} |

诊断是 **corr(Z_noise, y)**：leaky 票是噪声专家在同一点上背下来的，会跟 y 相关；OOF 票几乎不相关。π 跟着走——OOF 把质量放在 good 上，leaky 给噪声专家更多票。holdout MSE / optimism 是后果，小样本会抖；相关才是特征有没有漏。

梯度协议里同一 batch 的对齐余弦 = {leak['leaky_cos']:.2f}，probe/holdout 不相交才是 {leak['honest_cos']:.2f}。票和靶必须拆开。

能进 level-1 的列：OOF 预测、OOF 残差、OOF 方向分数 `⟨ĝ_m, ĝ_hold⟩`。不能塞：本窗 in-sample 损失、x 本身、时间下标 t——把 x 或 t 塞进门，就滑向 MoE / 时变门，不是 stacking。

---

## π 的变动（要花在刀刃上）

`TV(π_t, π_{{t+1}}) = ½‖π_{{t+1}}−π_t‖_1`。`path TV = Σ_t TV_t` 是这条轨迹的复杂度（和 switching 次数是一家）。

| 流 | 想要的 π |
|---|---|
| 半程换人 | 切换后 **局部** TV 尖峰，把质量从 a 搬到 b |
| 定常 | path TV 小；大 TV = 在拟合上一窗噪声 |
| clone | TV 可以有，但应拆冗余票，不是来回抖 |

本例切换流 Hedge+share path TV = {switch_tv:.2f}。同样方法在**定常**流上是 {stat_tv:.2f}——三个差不多好的专家时，Fixed-Share 每步掺 `α/K`，会在没事的时候烧 TV（复活税）。sticky λ=0.12 把执行器 TV 压下来。跟踪算法的 TV 应该花在切换点，而不是当成默认抖动。

离散 OSL：TV 是 **一次跳到顶点**（不平滑）。Hedge：连续、乘性。滑窗：窗口一满就跳。Fixed-Share：每步至少掺 `α/K`，底噪 TV 换复活。

**Two-timescale（给执行器）：** meta 可以用快的 π（跟踪），LR 用慢的 `π_slow ← (1-λ)π_slow + λ π_fast`。FWD 一直开着的时候，π 一抖 LR 就抖。λ=0.12 的 sticky 路径 path TV 更小、延迟稍长——这是执行器和跟踪之间的预算。

share α、窗宽 W、sticky λ 是 **三个花 TV 的旋钮**，不是三个新模型。α/W 管「半程换人复不复活」，λ 管「下一步 LR 跟不跟得上」。

---

## 和前面几张的关系

| 已经说清 | 这一张补的 |
|---|---|
| 线性 gain → 顶点 | 半程换人时顶点会 **锁死在旧人** |
| 有限步长 / 联合 GLS | 那是一张窗里怎么组合；这里是 π **跨窗怎么走** |
| honesty gap | 把 gap 拆成 OOF 特征泄漏 vs 先 update 再 score |
| Hedge+share 能跟踪 | share、窗、sticky 各自花掉多少 TV、换来多少 delay |

```bash
PYTHONPATH=. python3 scripts/run_agod_stack_track.py
PYTHONPATH=. python3 -m tests.test_agod_stack_track
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(packed, oof, path: Path):
    lines = []
    for key in ("hedge_share0", "hedge_share", "osl_disc", "osl_window", "osl_sgd"):
        r = packed[key]
        d = r["delay_from_switch"]
        ds = "inf" if d is None else str(d)
        pi = r["final_pi"]
        lines.append(
            f"{key.replace('_', '\\_')} & {r['preq_mean']:.3f} & {r['path_tv']:.2f} & "
            f"{ds} & {pi.get('a', 0):.2f} & {pi.get('b', 0):.2f} \\\\"
        )
    tex = (
        "% Half-stream switch / OOF / pi variation\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Tracking delay and path total variation after a mid-stream "
        "expert switch. Discrete OSL locks; Fixed-Share spends TV after the switch.}\n"
        "\\label{tab:agod-stack-track}\n"
        "\\begin{tabular}{lccccc}\\toprule\n"
        "method & preq & path TV & delay & $\\pi_a$ & $\\pi_b$ \\\\\n"
        "\\midrule\n"
        + "\n".join(lines)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    votes, y, names, mid = stream_switch()
    packed = run_switch_methods(votes, y, names, t_switch=mid)
    sticky = sticky_path(packed["hedge_share"]["traj"], names, lam=0.12)
    oof = oof_demo()
    leak = grad_leak_demo()
    sv, sy, sn = stream_stationary()
    stat = run_switch_methods(sv, sy, sn, t_switch=len(sy) // 2)
    switch_tv = packed["hedge_share"]["path_tv"]
    stat_tv = stat["hedge_share"]["path_tv"]

    slim = {}
    for k, v in packed.items():
        slim[k] = {kk: vv for kk, vv in v.items() if kk != "traj"}
    payload = {
        "switch": slim,
        "oof": oof,
        "grad_leak": leak,
        "stationary_hedge_tv": stat_tv,
        "switch_hedge_tv": switch_tv,
        "sticky_tv": path_tv(sticky),
        "mid": mid,
    }
    (OUT / "agod_stack_track.json").write_text(json.dumps(_jsonable(payload), indent=2))
    plot_board({"mid": mid}, packed, sticky, oof, OUT / "AGOD_Stack_Track_Board.png")
    write_docs(packed, oof, leak, stat_tv, switch_tv, OUT / "README.md")
    write_latex(packed, oof, OUT / "AGOD_stack_track_tables_only.tex")
    write_docs(packed, oof, leak, stat_tv, switch_tv, DOCS / "AGOD_stack_track.md")
    write_latex(packed, oof, DOCS / "AGOD_stack_track_tables_only.tex")
    ART.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUT / "AGOD_Stack_Track_Board.png", ART / "AGOD_Stack_Track_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    shutil.copy2(
        OUT / "AGOD_Stack_Track_Board.png",
        Path("/opt/cursor/artifacts/agod_stack_track_board.png"),
    )
    print("switch:", {k: {kk: slim[k][kk] for kk in ("preq_mean", "path_tv", "delay_from_switch", "locked")} for k in slim}, flush=True)
    print(
        "oof:",
        {k: oof[k] for k in ("pi_oof", "pi_leak", "corr_oof_noise", "corr_leak_noise", "optimism_oof", "optimism_leaky")},
        flush=True,
    )
    print(f"TV switch={switch_tv:.3f} stationary={stat_tv:.3f} sticky={path_tv(sticky):.3f}", flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
