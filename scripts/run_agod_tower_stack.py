#!/usr/bin/env python3
"""Three-tower votes → online π(t). Snapshot GLS vs Hedge vs leak vs MMoE.

  PYTHONPATH=. python3 scripts/run_agod_tower_stack.py
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

from agod.three_tower import CATS, NegPool
from agod.tower_stack import run_tower_stack, switch_events, switch_report

OUT = ROOT / "results" / "agod_tower_stack"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_tower_stack")
SEED = 2026
N_EACH = 18


def catalog(n: int = 48, seed: int = SEED) -> NegPool:
    rng = np.random.default_rng(seed)
    cats = list(CATS[:-1])
    items = []
    for k in range(n):
        items.append(
            {
                "parent_asin": f"i{k}",
                "category": cats[int(rng.integers(0, len(cats)))],
                "price": float(rng.uniform(12, 36)),
            }
        )
    return NegPool(items)


def plot_board(hon, leak, disc, t_switch: int, path: Path):
    traj = hon["traj"]
    ts = np.arange(len(traj))
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.1), facecolor="#f7f5f1")
    fig.suptitle(
        "Three-tower online stack: cat switch → π(t), not a GLS snapshot",
        fontsize=13,
        fontweight="bold",
    )

    ax = axes[0, 0]
    ax.axvline(t_switch, color="#999", ls="--", lw=0.9)
    ax.plot(ts, [r["snap"]["cat_hist"].get("Tools", 0) for r in traj], lw=2.0, label="cat_hist Tools")
    ax.plot(ts, [r["snap"]["cat_hist"].get("Sports", 0) for r in traj], lw=2.0, label="cat_hist Sports")
    ax.plot(ts, [r["votes"]["user"] for r in traj], lw=1.8, ls=":", label="user vote ⟨hist, pos⟩")
    ax.plot(ts, [r["votes"]["item"] for r in traj], lw=1.8, ls="--", label="item vote (price)")
    ax.set_title("Frozen portrait + stack votes")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[0, 1]
    ax.axvline(t_switch, color="#999", ls="--", lw=0.9)
    for name, rows, ls in (
        ("hedge π_user", [r["pi"]["user"] for r in traj], "-"),
        ("hedge π_item", [r["pi"]["item"] for r in traj], "-"),
        ("GLS snap π_user", [r["pi_gls"]["user"] for r in traj], ":"),
        ("osl_disc π_user", [r["pi"]["user"] for r in disc["traj"]], "--"),
    ):
        ax.plot(ts, rows, lw=2.0, ls=ls, label=name)
        ax.set_title("π(t): Hedge leaves user; GLS snapshot does not")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1, 0]
    ax.axvline(t_switch, color="#999", ls="--", lw=0.9)
    hon_loss = np.array([r["preq_loss"] for r in traj], float)
    leak_loss = np.array([r["preq_loss"] for r in leak["traj"]], float)
    k = 5
    kernel = np.ones(k) / k
    ax.plot(ts, np.convolve(hon_loss, kernel, mode="same"), lw=2.0, label="honest preq")
    ax.plot(ts, np.convolve(leak_loss, kernel, mode="same"), lw=2.0, label="leaky stacker preq")
    ax.set_title("Prequential loss (leaky is optimistic)")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    ax.axvline(t_switch, color="#999", ls="--", lw=0.9)
    for tw in ("user", "item", "neg"):
        ax.plot(ts, [r["lr"][tw] for r in traj], lw=2.0, label=tw)
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("Sticky LR×  (conf-damped, FWD on)")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.015,
        "user vote = cat_hist match  ·  item vote = price affinity  ·  leaky π scores in-sample  ·  not MMoE",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(hon, leak, disc, t_switch: int, path: Path):
    rh = switch_report(hon, t_switch)
    rl = switch_report(leak, t_switch)
    rd = switch_report(disc, t_switch)
    last = hon["traj"][-1]

    def fmt(pi):
        return ", ".join(f"{k}={pi[k]:.2f}" for k in ("user", "item", "neg"))

    md = f"""# 三塔票的在线 stacking：类目切换、泄漏、快照 GLS

上一张三塔把 `(cat_match, s_pos, s_neg)` **每步独立** GLS 一遍。y 全是 1 时，oracle 轴上的快照 `π` 几乎不动。缺的是 stacking 协议本身：

```
冻住画像 且 冻住 π_t
votes ← (snap, item, pool)     # 不准看 y
honest_step(stacker, votes, y) # 先计分再 update
sticky π → conf-damped LR      # FWD 开着
然后才 Welford / seen
```

这不是 MMoE。MMoE 是 `π(x)=softmax(W mean(U,I,N))`，门看输入。这里 `π_m(t)` 看的是 expert 的 one-step-ahead 损失。

---

## 票怎么造（类目频率 / 价格，不是原始塔向量）

| expert | 票（冻住的） | 类目一切换会怎样 |
|---|---|---|
| **user** | `⟨cat_hist, onehot(pos.cat)⟩` | 掉到 0，直到画像追上 |
| **item** | 价格亲和（和类目无关的 item 坐标） | 同价位时还活着 |
| **neg** | `1 − cos(U,N)` | 负样本池硬度 |

原始 `{{U,I,N}}` 不能当 GLS 库（会在 item 顶点塌掉）。在线 stacking 吃的是这些标量票，预测 y。

本流：前 {t_switch} 步 Tools，之后 Sports，价格带相同。user 票 {rh['vote_user_pre']:.2f} → {rh['vote_user_post']:.2f}。

---

## 谁在跟踪切换

| 方法 | π_item 前 → 后 | path TV | 延迟（π_item≥0.4） |
|---|---|---|---|
| Hedge + Fixed-Share（诚实） | {rh['pi_item_pre']:.2f} → {rh['pi_item_post']:.2f} | {rh['path_tv']:.2f} | {rh['delay_item']} |
| 离散 OSL | {rd['pi_item_pre']:.2f} → {rd['pi_item_post']:.2f} | {rd['path_tv']:.2f} | {rd['delay_item']} |
| 每步 GLS 快照 | TV={hon['path_tv_gls']:.3f} | 几乎不走 | — |

半程换人花的是单纯形 TV。快照 GLS 的 `π_user` 钉在约 2/3，**不是** `π(t)`。离散 OSL 是顶点（乱跳）。Hedge+share 在切换处把质量从 user 挪走。

---

## 两种泄漏

1. **画像泄漏**（他们 `CausalPosNeg` 的反面）：先 `update` 再计分。冷启动第一步 `cat_match` 从 0 变成 1。类目切换时泄漏只有 `1/n`，不是主戏。
2. **stacker 泄漏**：先 `π ← update` 再计分。in-sample 乐观。本流诚实 preq **{hon['mean_preq']:.3f}**，leaky stacker **{leak['mean_preq']:.3f}**（切换后 {rh['preq_post']:.3f} vs {rl['preq_post']:.3f}）。

终步诚实 Hedge `π`：{fmt(last['pi'])}。sticky actuator TV **{hon['path_tv_act']:.2f}** < 快 π TV **{hon['path_tv']:.2f}**（两套时间尺度，FWD 仍开）。

```bash
PYTHONPATH=. python3 -m tests.test_agod_tower_stack
PYTHONPATH=. python3 scripts/run_agod_tower_stack.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(hon, leak, disc, t_switch: int, path: Path):
    rh = switch_report(hon, t_switch)
    rd = switch_report(disc, t_switch)
    tex = (
        "% Three-tower online stack vs snapshot GLS\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Category switch: Hedge $\\pi(t)$ moves mass to the price "
        "item vote; per-step GLS snapshot almost does not. "
        "Leaky stacker prequential loss is optimistic.}\n"
        "\\label{tab:agod-tower-stack-switch}\n"
        "\\begin{tabular}{lcccc}\\toprule\n"
        "method & $\\pi_{\\mathrm{item}}$ pre & post & path TV & mean preq \\\\\n"
        "\\midrule\n"
        f"Hedge honest & {rh['pi_item_pre']:.2f} & {rh['pi_item_post']:.2f} & "
        f"{rh['path_tv']:.2f} & {hon['mean_preq']:.3f} \\\\\n"
        f"OSL discrete & {rd['pi_item_pre']:.2f} & {rd['pi_item_post']:.2f} & "
        f"{rd['path_tv']:.2f} & {disc['mean_preq']:.3f} \\\\\n"
        f"GLS snapshot TV & --- & --- & {hon['path_tv_gls']:.3f} & --- \\\\\n"
        f"leaky stacker & --- & --- & {leak['path_tv']:.2f} & {leak['mean_preq']:.3f} \\\\\n"
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
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
    if isinstance(obj, set):
        return list(obj)
    return obj


def _slim(payload):
    rows = []
    for r in payload["traj"]:
        rows.append(
            {
                "category": r["category"],
                "snap": r["snap"],
                "votes": r["votes"],
                "pi": r["pi"],
                "pi_act": r["pi_act"],
                "pi_gls": r["pi_gls"],
                "pi_mmoe": r["pi_mmoe"],
                "lr": r["lr"],
                "preq_loss": r["preq_loss"],
            }
        )
    return {
        "traj": rows,
        "final_snap": payload["final_snap"],
        "method": payload["method"],
        "leak": payload["leak"],
        "path_tv": payload["path_tv"],
        "path_tv_act": payload["path_tv_act"],
        "path_tv_gls": payload["path_tv_gls"],
        "mean_preq": payload["mean_preq"],
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    events, t_switch = switch_events(N_EACH)
    pool = catalog()
    hon = run_tower_stack(events, pool, method="hedge", leak="none", seed=SEED, share=0.10, eta=0.9)
    leak = run_tower_stack(events, pool, method="hedge", leak="stacker", seed=SEED, share=0.10, eta=0.9)
    disc = run_tower_stack(events, pool, method="osl_disc", leak="none", seed=SEED)
    out = {
        "honest": _slim(hon),
        "leaky_stacker": _slim(leak),
        "osl_disc": _slim(disc),
        "t_switch": t_switch,
        "switch_honest": switch_report(hon, t_switch),
        "switch_leaky": switch_report(leak, t_switch),
        "switch_disc": switch_report(disc, t_switch),
    }
    (OUT / "agod_tower_stack.json").write_text(json.dumps(_jsonable(out), indent=2))
    plot_board(hon, leak, disc, t_switch, OUT / "AGOD_Tower_Stack_Board.png")
    write_docs(hon, leak, disc, t_switch, OUT / "README.md")
    write_latex(hon, leak, disc, t_switch, OUT / "AGOD_tower_stack_tables_only.tex")
    write_docs(hon, leak, disc, t_switch, DOCS / "AGOD_tower_stack.md")
    write_latex(hon, leak, disc, t_switch, DOCS / "AGOD_tower_stack_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Tower_Stack_Board.png", ART / "AGOD_Tower_Stack_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    shutil.copy2(
        OUT / "AGOD_Tower_Stack_Board.png",
        Path("/opt/cursor/artifacts/agod_tower_stack_board.png"),
    )
    rh = out["switch_honest"]
    print("hist", hon["final_snap"]["cat_hist"], "n", hon["final_snap"]["review_cnt"], flush=True)
    print("user vote", round(rh["vote_user_pre"], 3), "->", round(rh["vote_user_post"], 3), flush=True)
    print("π_item", round(rh["pi_item_pre"], 3), "->", round(rh["pi_item_post"], 3), "tv", round(hon["path_tv"], 3), "gls tv", round(hon["path_tv_gls"], 3), flush=True)
    print("preq honest", round(hon["mean_preq"], 3), "leaky", round(leak["mean_preq"], 3), flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
