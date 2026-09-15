#!/usr/bin/env python3
"""Three-tower prototype: causal cat_hist / review_cnt, neg pool, AGOD π→LR.

  PYTHONPATH=. python3 scripts/run_agod_three_tower.py
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

from agod.online_portrait import OnlinePortrait
from agod.three_tower import (
    CATS,
    NegPool,
    honest_three_tower_step,
    pool_eval,
    score_triplet,
    user_tower,
)

OUT = ROOT / "results" / "agod_three_tower"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_three_tower")
SEED = 2026


def catalog(n: int = 40, seed: int = SEED) -> list[dict]:
    rng = np.random.default_rng(seed)
    cats = list(CATS[:-1])
    items = []
    for k in range(n):
        cat = cats[int(rng.integers(0, len(cats)))]
        items.append(
            {
                "parent_asin": f"i{k}",
                "category": cat,
                "price": float(rng.uniform(8, 90)),
            }
        )
    return items


def user_events():
    """Mostly Tools, then a Tools pos and a Fashion pos."""
    hist = [{"category": "Tools", "price": 18.0 + i, "parent_asin": f"h{i}", "y": 1} for i in range(8)]
    hist += [{"category": "Sports", "price": 40.0, "parent_asin": "hS", "y": 1}]
    return hist


def run_stream(items: list[dict], seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    pool = NegPool(items)
    portrait = OnlinePortrait()
    seen: set[str] = set()
    traj = []
    for ev in user_events():
        y = float(ev["y"])
        pos = {k: ev[k] for k in ("parent_asin", "category", "price")}
        step = honest_three_tower_step(
            portrait, pos, pool, seen, y, rng, k_neg=4, hard=False
        )
        traj.append({**step, "category": ev["category"]})
    snap = portrait.snapshot()
    pos = {"parent_asin": "probe", "category": "Tools", "price": 20.0}
    ev = pool_eval(snap, pos, pool, seen, rng, k_neg=6)
    return {"traj": traj, "final_snap": snap, "pool": ev, "seen": len(seen)}


def plot_board(payload: dict, path: Path):
    traj = payload["traj"]
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.0), facecolor="#f7f5f1")
    fig.suptitle(
        "Three-tower: cat_hist / review_cnt / neg-pool → AGOD π→LR",
        fontsize=13,
        fontweight="bold",
    )

    ax = axes[0, 0]
    ts = np.arange(len(traj))
    ax.plot(ts, [r["snap"]["cat_hist"].get("Tools", 0) for r in traj], lw=2.0, label="cat_hist Tools")
    ax.plot(ts, [r["snap"]["cat_hist"].get("Sports", 0) for r in traj], lw=2.0, label="cat_hist Sports")
    ax.plot(ts, [r["trip"]["cat_match_pos"] for r in traj], lw=1.8, ls=":", label="⟨hist, pos.cat⟩")
    ax.plot(ts, [r["trip"]["conf"] for r in traj], lw=2.0, ls="--", label="review_cnt confidence")
    ax.set_title("Frozen portrait (emit then update)")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    names = ["random_global", "seen_excluded", "same_cat_hard"]
    gaps = [payload["pool"][n]["gap"] for n in names]
    cont = [payload["pool"][n]["contam"] for n in names]
    x = np.arange(len(names))
    ax.bar(x - 0.18, gaps, 0.36, label="pos−neg gap")
    ax.bar(x + 0.18, cont, 0.36, label="in-seen contamination")
    ax.set_xticks(x)
    ax.set_xticklabels(["global", "seen-excl", "same-cat hard"], fontsize=8)
    ax.set_title("Negative pool")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    for tw in ("user", "item", "neg"):
        ax.plot(ts, [r["scores"][tw] for r in traj], lw=2.0, label=rf"$s_{{{tw}}}$")
        ax.plot(ts, [r["pi"][tw] for r in traj], lw=1.4, ls="--", label=rf"$\pi_{{{tw}}}$")
    ax.set_title("AGOD votes s and alignment π")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=7, ncol=2)

    ax = axes[1, 1]
    for tw in ("user", "item", "neg"):
        ax.plot(ts, [r["lr"][tw] for r in traj], lw=2.0, label=tw)
    ax.axhline(1.0, color="#999", ls="--", lw=0.8)
    ax.set_title("LR×  (damped by review_cnt confidence)")
    ax.set_xlabel("t")
    ax.legend(frameon=False, fontsize=8)

    fig.text(
        0.5,
        0.015,
        "U = shrunk cat_hist  ·  match = ⟨hist, item.cat⟩  ·  n→0 ⇒ uniform + LR→1  ·  not MMoE",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(payload: dict, path: Path):
    pool = payload["pool"]
    last = payload["traj"][-1]
    snap = payload["final_snap"]

    def fmt_pi(pi):
        return ", ".join(f"{k}={pi[k]:.2f}" for k in ("user", "item", "neg"))

    md = f"""# 三塔原型：类目频率、评论次数、负样本池 → AGOD

不做特征系统。在线状态只有因果画像：`cat_hist`、`review_cnt`、价格 μ/σ。三塔是

| 塔 | 输入（冻住的） | 向量 |
|---|---|---|
| **U user** | `cat_hist` + `review_cnt` | 类目单纯形向均匀先验收缩，收缩量由 n 决定 |
| **I item** | 当前商品类目 + 价格 | one-hot + 价格，目录不是评分 |
| **N neg** | 同一 item 塔，seen-excluded 池 | 他们 `CausalPosNeg` 的负样本 |

分数 `s = cos(U,I) − mean_k cos(U,N_k)`。AGOD **不**把原始 {{U,I,N}} 当 GLS 库（holdout 在它们张成的空间里会在 item 顶点塌掉）。票是正交打包的标量，holdout 只看 y，再 `π→LR`。顺序仍是 **先冻住画像再计分，最后才 update Welford**。

---

## 类目频率怎么刻画

`cat_hist` 是用户塔在类目单纯形上的**方向**。和当前商品的对齐就是

```
match = ⟨cat_hist, onehot(item.cat)⟩
```

也就是「过去有多大比例的评论落在这个类目」。正样本 match 应高于负样本；same-cat hard 负样本会把这个差压小——这就是负样本池硬度。

终盘画像：`{snap['cat_hist']}`，`review_cnt={snap['review_cnt']}`。

---

## 评论次数怎么刻画

`review_cnt` 不是又一个稀疏 ID，是用户塔的**置信度**：

```
conf = 1 − exp(−n / n0)
U = (1−conf)·uniform + conf·cat_hist
LR_damped = 1 + conf · (LR − 1)
```

n→0：U 塌成均匀，AGOD 不许猛调 LR（冷启动）。n 大：相信类目方向，才允许 π 离开 1/3、LR 离开 1。这和 stacking 里 `N_eff` / Kish 是同一类量。

---

## 负样本池

| pool | pos−neg gap | in-seen | cat_match neg |
|---|---:|---:|---:|
| random global | {pool['random_global']['gap']:.3f} | {pool['random_global']['contam']:.2f} | {pool['random_global']['cat_match_neg']:.2f} |
| seen-excluded | {pool['seen_excluded']['gap']:.3f} | {pool['seen_excluded']['contam']:.2f} | {pool['seen_excluded']['cat_match_neg']:.2f} |
| same-cat hard | {pool['same_cat_hard']['gap']:.3f} | {pool['same_cat_hard']['contam']:.2f} | {pool['same_cat_hard']['cat_match_neg']:.2f} |

他们流水线是 seen-excluded + 全局 10k 池（上架时间不过滤）。hard 负样本 gap 更小，才是在考 `cat_hist` 而不是考「类目不同」。

---

## 回到 AGOD 对齐

三塔向量负责 **打分**。AGOD 的 expert 票是打到正交轴上的标量（否则 GLS 会在 item 顶点塌掉）。holdout 只看 y：

```
π = GLS / direction-match(vecs, g_hold)
LR = π_to_lr(π)   # 再乘 conf
```

holdout 是 **y 决定的 oracle 轴**（正样本要 user/item 轴，负样本要 neg 轴），不是把 {{U,I,N}} 再投影回去——那会在 item 上塌成顶点。expert 票是正交打包的标量：user=`⟨cat_hist, pos.cat⟩`，item=`cos(U,I)`，neg=`cos(U,N)`。所以 **s_user 跟着类目频率走**，π 是 GLS 组合，LR 再乘 `conf(review_cnt)`。终步 `π`：{fmt_pi(last['pi'])}。这不是 MMoE：门不看 x，看 holdout。MMoE 的 `softmax(W mean(experts))` 只对照，不进这条环。

```bash
PYTHONPATH=. python3 -m tests.test_agod_three_tower
PYTHONPATH=. python3 scripts/run_agod_three_tower.py
```
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(payload: dict, path: Path):
    pool = payload["pool"]
    lines = []
    for name in ("random_global", "seen_excluded", "same_cat_hard"):
        r = pool[name]
        lines.append(
            f"{name.replace('_', '\\_')} & {r['gap']:.3f} & {r['contam']:.2f} & "
            f"{r['cat_match_neg']:.2f} \\\\"
        )
    tex = (
        "% Three-tower neg pool\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Negative pool: seen-excluded vs same-category hard negs. "
        "cat\\_hist match on negatives is the hardness diagnostic.}\n"
        "\\label{tab:agod-three-tower-neg}\n"
        "\\begin{tabular}{lccc}\\toprule\n"
        "pool & pos--neg gap & in-seen & cat match (neg) \\\\\n"
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
    if isinstance(obj, set):
        return list(obj)
    return obj


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    items = catalog()
    payload = run_stream(items)
    slim_traj = []
    for r in payload["traj"]:
        slim_traj.append(
            {
                "category": r["category"],
                "snap": r["snap"],
                "trip": r["trip"],
                "pi": r["pi"],
                "scores": r["scores"],
                "lr": r["lr"],
                "contam": r["contam"],
            }
        )
    out = {
        "traj": slim_traj,
        "final_snap": payload["final_snap"],
        "pool": payload["pool"],
        "seen": payload["seen"],
    }
    (OUT / "agod_three_tower.json").write_text(json.dumps(_jsonable(out), indent=2))
    plot_board(payload, OUT / "AGOD_Three_Tower_Board.png")
    write_docs(payload, OUT / "README.md")
    write_latex(payload, OUT / "AGOD_three_tower_tables_only.tex")
    write_docs(payload, DOCS / "AGOD_three_tower.md")
    write_latex(payload, DOCS / "AGOD_three_tower_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Three_Tower_Board.png", ART / "AGOD_Three_Tower_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    shutil.copy2(
        OUT / "AGOD_Three_Tower_Board.png",
        Path("/opt/cursor/artifacts/agod_three_tower_board.png"),
    )
    print("final cat_hist", payload["final_snap"]["cat_hist"], "n", payload["final_snap"]["review_cnt"], flush=True)
    print("pool", {k: {kk: round(vv, 3) for kk, vv in v.items() if kk != "n"} for k, v in payload["pool"].items()}, flush=True)
    print("last π", {k: round(v, 3) for k, v in payload["traj"][-1]["pi"].items()}, flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
