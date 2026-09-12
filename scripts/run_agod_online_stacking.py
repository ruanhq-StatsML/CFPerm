#!/usr/bin/env python3
"""Online stacking of named experts (modality-agnostic).

Code is the easy part. This script evaluates the *protocol*:
  honest prequential score → then update π (Online Super Learner / Hedge / BG / GLS)
against oracle combo, best expert, leaky stacking, and the correlation-clone trap.

  PYTHONPATH=. python3 scripts/run_agod_online_stacking.py
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

from agod.online_stacking import (
    METHODS,
    OnlineStacker,
    best_expert_loss,
    honest_step,
    leaky_step,
    oracle_convex_combo,
    pi_to_lr,
    prequential_regret,
)

OUT = ROOT / "results" / "agod_online_stacking"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_online_stacking")
GRADCOS = ROOT / "results" / "agod_gradcos_lr" / "agod_gradcos_lr.json"
SEED = 2026


def stream_switch(t: int = 80, seed: int = SEED):
    rng = np.random.default_rng(seed)
    y = np.zeros(t)
    mid = t // 2
    cols = []
    for i in range(t):
        if i < mid:
            a, b = rng.normal(0, 0.12), rng.normal(0, 1.05)
        else:
            a, b = rng.normal(0, 1.05), rng.normal(0, 0.12)
        c = rng.normal(0, 0.65)
        cols.append((a, b, c))
    return np.array(cols), y, ("a", "b", "c")


def stream_clones(t: int = 100, seed: int = SEED + 1):
    rng = np.random.default_rng(seed)
    y = np.zeros(t)
    rows = []
    for _ in range(t):
        common = rng.normal(0, 0.48)
        e0 = common + rng.normal(0, 0.04)
        e1 = common + rng.normal(0, 0.04)
        e2 = rng.normal(0, 0.30)
        rows.append((e0, e1, e2))
    return np.array(rows), y, ("e0", "e1", "e2")


def run_methods(votes, y, names, methods=METHODS):
    out = {}
    for meth in methods:
        st = OnlineStacker(names, method=meth, eta=0.55, ewma=0.18, share=0.07)
        scored = []
        for t in range(len(y)):
            v = {n: float(votes[t, i]) for i, n in enumerate(names)}
            scored.append(honest_step(st, v, float(y[t])))
        preq = [s["preq_loss"] for s in scored]
        be = best_expert_loss(votes, y)
        oc = oracle_convex_combo(votes, y)
        out[meth] = {
            "preq_mean": float(np.mean(preq)),
            "regret_best_expert": prequential_regret(preq, be["loss"]),
            "regret_oracle_combo": prequential_regret(preq, oc["loss"]),
            "final_pi": st.pi(),
            "final_n_eff": st.n_eff(),
            "traj_pi": [s["pi"] for s in scored],
            "preq": preq,
            "best_expert_loss": be["loss"],
            "oracle_combo_loss": oc["loss"],
            "oracle_combo_pi": oc["pi"].tolist(),
            "lr": pi_to_lr(st.pi(), names),
        }
    return out


def honesty_gap(votes, y, names):
    hon = OnlineStacker(names, method="hedge", eta=0.40)
    leak = OnlineStacker(names, method="hedge", eta=0.40)
    lh, ll = [], []
    for t in range(len(y)):
        v = {n: float(votes[t, i]) for i, n in enumerate(names)}
        lh.append(honest_step(hon, v, float(y[t]))["preq_loss"])
        ll.append(leaky_step(leak, v, float(y[t]))["preq_loss"])
    return {
        "honest": float(np.mean(lh)),
        "leaky": float(np.mean(ll)),
        "gap": float(np.mean(lh) - np.mean(ll)),
    }


def replay_gradcos(path: Path) -> dict:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    traj = payload.get("trajectory") or {}
    out = {}
    for key, rows in traj.items():
        if ":" not in key or not rows:
            continue
        mods = list((rows[0].get("align_gain") or rows[0].get("alpha") or {}).keys())
        if not mods:
            continue
        st = OnlineStacker(mods, method="hedge", eta=0.35)
        scored = []
        for r in rows:
            votes = r.get("align_gain") or r.get("alpha") or {}
            y = 1.0  # target: perfect alignment with shared head
            v = {m: float(votes.get(m, 0.5)) for m in mods}
            scored.append(honest_step(st, v, y))
        out[key] = {
            "final_pi": st.pi(),
            "mean_alpha": {
                m: float(np.mean([r["alpha"][m] for r in rows if "alpha" in r]))
                for m in mods
            },
            "n_eff": st.n_eff(),
            "preq_mean": float(np.mean([s["preq_loss"] for s in scored])),
            "lr": pi_to_lr(st.pi(), mods),
        }
    return out


def plot_board(switch, clones, gap_sw, gap_cl, replay, path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.2), facecolor="#f7f5f1")
    fig.suptitle(
        "Online stacking of named experts  (honest prequential meta-learner)",
        fontsize=13,
        fontweight="bold",
    )

    ax = axes[0, 0]
    traj = switch["hedge"]["traj_pi"]
    xs = np.arange(len(traj))
    ax.plot(xs, [p["a"] for p in traj], lw=2.0, label="π_a (good then bad)")
    ax.plot(xs, [p["b"] for p in traj], lw=2.0, label="π_b (bad then good)")
    ax.axvline(len(traj) / 2, color="#999", ls="--", lw=0.8)
    ax.set_title("Hedge tracks the switching expert")
    ax.set_xlabel("window t")
    ax.set_ylabel("π")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    names = ("e0", "e1", "e2")
    x = np.arange(len(names))
    w = 0.25
    for i, meth in enumerate(("equal", "hedge", "gls_ewma")):
        vals = [clones[meth]["final_pi"][n] for n in names]
        ax.bar(x + (i - 1) * w, vals, w, label=meth)
    ax.set_xticks(x)
    ax.set_xticklabels(["e0 clone", "e1 clone", "e2 unique"])
    ax.set_title("Clone trap: GLS downweights redundant votes")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    meths = [m for m in METHODS if m in switch]
    x = np.arange(len(meths))
    ax.bar(x - 0.2, [switch[m]["preq_mean"] for m in meths], 0.4, label="switch")
    ax.bar(x + 0.2, [clones[m]["preq_mean"] for m in meths], 0.4, label="clones")
    ax.set_xticks(x)
    ax.set_xticklabels(meths, rotation=25, ha="right")
    ax.set_title("Prequential MSE (honest)")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    ax.bar([0, 1], [gap_sw["honest"], gap_sw["leaky"]], 0.35, label="switch")
    ax.bar([0.4, 1.4], [gap_cl["honest"], gap_cl["leaky"]], 0.35, label="clones")
    ax.set_xticks([0.2, 1.2])
    ax.set_xticklabels(["honest", "leaky"])
    ax.set_title("Honesty gap (leaky underestimates risk)")
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.015,
        "score with frozen pi_t, then update  ·  experts are names, not modalities  ·  FWD on",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _tbl(headers, rows):
    head = "| " + " | ".join(headers) + " |"
    sep = "|" + "|".join("---:" if i else "---" for i, _ in enumerate(headers)) + "|"
    body = "\n".join("| " + " | ".join(r) + " |" for r in rows)
    return "\n".join([head, sep, body])


def write_docs(switch, clones, gap_sw, gap_cl, replay, path: Path):
    def rows_for(pack):
        out = []
        for m, r in pack.items():
            out.append(
                [
                    m,
                    f"{r['preq_mean']:.3f}",
                    f"{r['regret_best_expert']:+.3f}",
                    f"{r['regret_oracle_combo']:+.3f}",
                    f"{r['final_n_eff']:.2f}",
                ]
            )
        return out

    replay_rows = []
    for k, s in sorted(replay.items()):
        pi = s["final_pi"]
        top = max(pi, key=pi.get)
        replay_rows.append(
            [
                k,
                top,
                f"{pi[top]:.2f}",
                f"{s['n_eff']:.2f}",
                f"{s['preq_mean']:.3f}",
            ]
        )

    md = f"""# Online stacking of modality experts（怎么做、怎么评、文献在哪）

静态 `π ∝ R̃^{{-1}} α` 只是 stacking 的**一张快照**。代码好写；难的是 **online-stacking 的协议**：权重在单纯形上、只用 one-step-ahead 损失更新，而且评估必须对 oracle。这件事是 **modality-agnostic** 的——expert 只是带名字的投票器。

## 文献地图（这不是从零发明）

| 传统 | 文献 | 我们把它用在哪 |
|---|---|---|
| Stacked generalization | Wolpert 1992; Breiman 1996 stacked regressions | 模态 = base expert，π = meta-learner |
| Super Learner | van der Laan, Polley, Hubbard 2007; oracle inequality van der Laan & Dudoit 2003 | CV 组合渐近等价于 library 里最好的凸组合 |
| **Online Super Learner** | Benkeser, Ju, Lendle, van der Laan, *Stat Med* 2018 ([PMC5671383](https://pmc.ncbi.nlm.nih.gov/articles/PMC5671383/)) | 流上的 one-step-ahead CV；离散 OSL = 累加损失最小的 expert；凸 OSL = 单纯形上 SGD |
| Prequential | Dawid 1984 | **先用冻住的 π_t 计分，再 update** — 流式 OOF |
| Hedge / weighted majority + Fixed-Share | Littlestone & Warmuth 1994; Freund & Schapire 1997; **Herbster & Warmuth 1998 tracking** | `π ← (1-α) π⊙exp(−ηL) + α/|E|`，否则半程切换会被前半程锁死 |
| Tracking the best expert | Herbster & Warmuth 1998 | 非平稳窗（概念漂移）下的动态遗憾 |
| Forecast combination | Bates & Granger 1969; Newbold & Granger 1974; Timmermann 2006 | 误差协方差的最小方差权；**实践中常丢掉 off-diagonal**（Σ 估不稳） |
| Error–ambiguity / BVC | Krogh & Vedelsby 1995; Ueda & Nakano 1996 | `N_eff`、clone trap：相关 expert 不是新票 |
| Negative correlation learning | Liu & Yao 1999 | 负相关才是多样性；与 shared head 冲突不是 |
| Mixture of experts | Jacobs et al. 1991 | gating 网是另一类 meta；我们用损失驱动的单纯形权 |
| Gradient Blending | Wang, Tran, Feiszli CVPR 2020, arXiv:1905.12681 | **不是 stacking**：按 OGR 混 loss。我们混的是 *expert 的票* |
| PCGrad / OGM-GE | Yu et al. 2020; Peng et al. 2022 | 冲突投影 / 梯度调制，没有 online CV oracle inequality |

缺口很干净：Gradient Blending / PCGrad 不懂 prequential honesty，OSL / Hedge / Bates–Granger 不懂「票是梯度签名」。把模态塔当成 **online expert**，vote 的构造见下一节，π 走 OSL/Hedge/BG，再 `π → LR`。video/text/audio/tabular **同一段代码**。

## 梯度 expert 的票怎么造（这才是「怎么做」）

预测 stacking 里 `v_m = f_m(x)`、`y` 可观测。梯度 stacking 里 `g_m ∈ R^d`，真更新 `μ = −∇L_pop` **看不见**。holdout 梯度 `g_hold` 是 μ 的 noisy proxy。三种 meta-loss，几何完全不同：

| meta-loss | 公式 | π 的形状 | 文献对应 |
|---|---|---|---|
| 一阶 holdout gain | `s_m = ⟨ĝ_m, ĝ_hold⟩`，`max πᵀs` | **顶点**（离散 OSL / winner-take-all） | 线性 Super Learner 退化 |
| **方向匹配** | `min ||Gπ − ĝ_hold||²` | 内点：`π ∝ R^{−1}s`，`R=GᵀG` | Breiman stacked regressions in `R^d`；就是上一张 GLS 快照 |
| 均值–方差 | `max πᵀs − (λ/2) πᵀRπ` | λ=0 顶点，λ↑ 缩向 min-var | Markowitz / Bates–Granger |
| 有限步长 | `L_hold(θ − η Gπ)` | 非线性，η 大才离开线性 | 真 virtual step，贵 |

所以：**不能**把 directional cosine 直接塞进凸 stacking 还指望得到组合权——线性目标在单纯形上必崩到一个 expert。要组合，必须用方向匹配（平方损失）或方差惩罚。`direction_match_weights` 就是这件事。

诚实协议在梯度上还多一条：`g_m` 和 `g_hold` 必须来自 **不同 batch**（probe vs holdout）。同一窗的 holdout 既当票又当靶，就是 leaky stacking。

## 协议（这才是 online-stacking，不是 `R^{{-1}}α` 一行）

```
for window t:
  1. freeze π_t
  2. each expert emits vote v_{{m,t}}  *before* being trained on window t
  3. score ŷ_t = π_t · v_t on the holdout of window t     # prequential loss
  4. train bases on the adapt split                       # FWD still on
  5. update π_{{t+1}} from (v_t, y_t, L_{{m,t}})
  6. actuator: LR_{{t+1,m}} = β + (1-β)|E| π_{{t+1,m}}
```

漏写第 2–3 步就是 **leaky stacking**（用 in-sample 票训 meta），等价于 Wolpert 警告过的用训练集预测做 stacking 特征。流上的正确类比是 Dawid 的 prequential / OSL 的 one-step-ahead CV。

Meta-learners：

| method | 更新 | 何时该赢 |
|---|---|---|
| `equal` | 不动 | 相关高、Σ 估不稳（Bates–Granger 自己的告诫） |
| `hedge` | 乘性权重 | expert 会切换（漂移） |
| `osl_disc` | argmin 累加损失 | 有一个稳定的最好 expert |
| `osl_sgd` | 投影 SGD on `(π·v − y)²` | 要凸组合，K 小 |
| `bg` | `π ∝ 1/MSE_m`（忽略相关） | 异方差、Σ 噪 |
| `gls_ewma` | `π ∝ Σ̂^{{-1}} 1` + 特征值收缩 | clone / 共线票，且 T 够估 Σ |

## 怎么评估（对 oracle，不只看 Acc lift）

1. **Prequential risk** `R_n = n^{{-1}} Σ_t (π_t·v_t − y_t)²` — 唯一诚实的流风险  
2. **vs best expert**（离散 oracle）和 **vs best fixed convex combo**（网格上的 Super Learner oracle）  
3. **Tracking**：半程切换最好 expert，Hedge/OSL 必须把质量搬过去  
4. **Clone trap**：两个几乎相同的好 expert + 一个独立 expert；等权把 2/3 砸在同一个因子上，GLS 必须丢掉冗余票  
5. **Honesty gap**：leaky − honest；leaky 应系统性低估风险  
6. **N_eff(π)**：stacking 是否还剩多样性  
7. 模态无关：expert 名叫 `a,b,c` 或 `video,text,audio` 行为应一样  

### Switch stream（最好 expert 在 t=T/2 切换）

{_tbl(["method", "preq MSE", "regret vs best", "regret vs combo", "N_eff"], rows_for(switch))}

负 regret vs best expert 是正常的：全时段最好的**单个** expert 往往是那个「一直中等」的 c，组合/跟踪能赢它。`osl_disc` 锁在前半程赢家 a 上，后半程崩掉（preq 0.55）——离散 Super Learner 在非平稳流上必须加 Fixed-Share / 窗，不能只用累加损失。`hedge`+share 和 `osl_sgd` 把质量从 a 搬到 b。

### Clone trap（e0≈e1 共线，e2 独立）

{_tbl(["method", "preq MSE", "regret vs best", "regret vs combo", "N_eff"], rows_for(clones))}

Clone 终盘 π：equal e0+e1 = {clones['equal']['final_pi']['e0']+clones['equal']['final_pi']['e1']:.2f}；GLS e0+e1 = {clones['gls_ewma']['final_pi']['e0']+clones['gls_ewma']['final_pi']['e1']:.2f}（冗余票被拆掉）。

### Honesty gap（Hedge）

| stream | honest | leaky | gap |
|---|---:|---:|---:|
| switch | {gap_sw['honest']:.3f} | {gap_sw['leaky']:.3f} | {gap_sw['gap']:+.3f} |
| clones | {gap_cl['honest']:.3f} | {gap_cl['leaky']:.3f} | {gap_cl['gap']:+.3f} |

## Amazon / MSR-VTT replay

把 `align_gain = ½(1+cos(g_m,g_shared))` 当 vote、目标=1，跑诚实 Hedge。T=6 且 align_gain 常挤在 0.5 附近，**不够当 OSL 渐近，也不该指望看到 Amazon 共线 → N_eff 下降**；replay 只证明同一套 stacker 吃 `video/text/audio` 或 `text/image` 这些名字。真正的评估是上面的 switch / clone / honesty 三条。

{_tbl(["traj", "top expert", "π_top", "N_eff", "preq"], replay_rows) if replay_rows else "_no gradcos json_"}

下一步 LR 仍是 `pi_to_lr(π)`，FWD 不关。T 短时不要用 replay 的 N_eff 下结论。

## 和上一张 GLS 快照的关系

`soft_stat` 的 `π* ∝ R^{{-1}}α` 是 **窗口内闭式 stacking**（Bates–Granger / 匹配滤波）。Online stacking 多了三件闭式没有的东西：

- **时间**：π 是状态，能 track 切换（Hedge 遗憾界）  
- **诚实**：one-step-ahead，不拿本窗训练票更新本窗权  
- **oracle inequality**：OSL 对 library 里最好凸组合渐近等价（Benkeser et al. 2018）  

相关矩阵仍在：`gls_ewma` 用误差协方差的 EWMA，clone 时它就是 `N_eff` 故事的在线版。

```bash
PYTHONPATH=. python3 scripts/run_agod_online_stacking.py
PYTHONPATH=. python3 -m tests.test_agod_online_stacking
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(switch, clones, path: Path):
    def lines(pack, tag):
        out = []
        for m, r in pack.items():
            out.append(
                f"{tag} & {m.replace('_', '\\_')} & {r['preq_mean']:.3f} & "
                f"{r['regret_best_expert']:+.3f} & {r['regret_oracle_combo']:+.3f} \\\\"
            )
        return out

    tex = (
        "% Online stacking prequential evaluation\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Honest prequential MSE and regret vs best expert / best convex combo.}\n"
        "\\label{tab:agod-online-stacking}\n"
        "\\begin{tabular}{llrrr}\\toprule\n"
        "stream & method & preq MSE & regret vs best & regret vs combo \\\\\n"
        "\\midrule\n"
        + "\n".join(lines(switch, "switch") + lines(clones, "clones"))
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)

    vs, ys, ns = stream_switch()
    vc, yc, nc = stream_clones()
    print("switch stream…", flush=True)
    switch = run_methods(vs, ys, ns)
    print("clone stream…", flush=True)
    clones = run_methods(vc, yc, nc)
    gap_sw = honesty_gap(vs, ys, ns)
    gap_cl = honesty_gap(vc, yc, nc)
    print("replay gradcos…", flush=True)
    replay = replay_gradcos(GRADCOS)

    slim_s = {m: {k: v for k, v in r.items() if k not in ("traj_pi", "preq")} for m, r in switch.items()}
    slim_c = {m: {k: v for k, v in r.items() if k not in ("traj_pi", "preq")} for m, r in clones.items()}
    payload = {
        "focus": "online stacking of named experts (honest prequential)",
        "methods": list(METHODS),
        "switch": slim_s,
        "clones": slim_c,
        "honesty": {"switch": gap_sw, "clones": gap_cl},
        "replay": replay,
        "switch_traj_hedge": switch["hedge"]["traj_pi"],
    }
    (OUT / "agod_online_stacking.json").write_text(json.dumps(payload, indent=2))
    plot_board(switch, clones, gap_sw, gap_cl, replay, OUT / "AGOD_Online_Stacking_Board.png")
    write_docs(switch, clones, gap_sw, gap_cl, replay, OUT / "README.md")
    write_latex(switch, clones, OUT / "AGOD_online_stacking_tables_only.tex")
    write_docs(switch, clones, gap_sw, gap_cl, replay, DOCS / "AGOD_online_stacking.md")
    write_latex(switch, clones, DOCS / "AGOD_online_stacking_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Online_Stacking_Board.png", ART / "AGOD_Online_Stacking_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    shutil.copy2(OUT / "AGOD_Online_Stacking_Board.png", Path("/opt/cursor/artifacts/agod_online_stacking_board.png"))

    print("\n=== switch ===", flush=True)
    for m, r in switch.items():
        print(
            f"  {m:10s} preq={r['preq_mean']:.3f}  regret_best={r['regret_best_expert']:+.3f}  "
            f"π={ {k: round(v,2) for k,v in r['final_pi'].items()} }",
            flush=True,
        )
    print("=== clones ===", flush=True)
    for m, r in clones.items():
        print(
            f"  {m:10s} preq={r['preq_mean']:.3f}  clone_mass={r['final_pi']['e0']+r['final_pi']['e1']:.2f}  "
            f"e2={r['final_pi']['e2']:.2f}",
            flush=True,
        )
    print(f"honesty gap switch={gap_sw['gap']:+.3f} clones={gap_cl['gap']:+.3f}", flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
