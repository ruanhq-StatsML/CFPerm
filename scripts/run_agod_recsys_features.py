#!/usr/bin/env python3
"""Causal recsys features (Amazon MMoE pipeline) vs leaky current-review.

  PYTHONPATH=. python3 scripts/run_agod_recsys_features.py
"""
from __future__ import annotations

import json
import shutil
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from agod.recsys_features import (
    agod_concat_blob,
    build_item_text,
    build_user_text,
    emit_causal_rows,
    leak_report,
    leaky_user_feat,
    shortcut_gap,
)

OUT = ROOT / "results" / "agod_recsys_features"
DOCS = ROOT / "docs" / "agod"
ART = Path("/opt/cursor/artifacts/agod_recsys_features")
SEED = 2026


def _ts(d: date) -> int:
    return int(datetime(d.year, d.month, d.day, tzinfo=timezone.utc).timestamp() * 1000)


def toy_users(n_users: int = 24, seed: int = SEED) -> list[dict]:
    rng = np.random.default_rng(seed)
    cats = ["Tools", "Sports", "Fashion", "Electronics"]
    rows = []
    start = date(2023, 1, 5)
    for u in range(n_users):
        n = int(rng.integers(3, 8))
        taste = rng.random()
        for k in range(n):
            d = start + timedelta(days=int(rng.integers(0, 280)))
            good = int(rng.random() < 0.35 + 0.5 * taste)
            text = (
                "love this excellent perfect item"
                if good
                else "bad terrible broken poor quality"
            )
            rows.append(
                {
                    "user_id": f"u{u}",
                    "parent_asin": f"I{int(rng.integers(0, 80))}",
                    "main_category": str(cats[int(rng.integers(0, len(cats)))]),
                    "price": float(rng.uniform(8, 90)),
                    "product_title": f"sku-{k}",
                    "review_title": "ok" if good else "no",
                    "review_text": text,
                    "sort_timestamp": _ts(d),
                    "event_date": d,
                    "label_good": good,
                }
            )
    return rows


def run_all(raw: list[dict]) -> dict:
    by_user: dict[str, list] = {}
    for r in raw:
        by_user.setdefault(r["user_id"], []).append(r)
    emitted = []
    leaky_hits = 0
    for rows in by_user.values():
        emitted.extend(emit_causal_rows(rows))
        ordered = sorted(rows, key=lambda x: x["sort_timestamp"])
        for i, r in enumerate(ordered):
            feat = leaky_user_feat(ordered[: i + 1])
            leaky_hits += int((feat.get("review_cnt") or 0) != i)
    causal = leak_report(emitted)
    gap = shortcut_gap(raw)
    splits = causal["splits"]
    example = next(r for r in emitted if r["user_feat"]["review_cnt"] >= 2)
    return {
        "causal": causal,
        "leaky_current_rate": float(leaky_hits / max(len(emitted), 1)),
        "shortcut": gap,
        "example_user_text": build_user_text(example["user_feat"]),
        "example_item_text": build_item_text(example),
        "example_agod_blob_chars": len(agod_concat_blob(example)),
        "n_users": len(by_user),
        "split_cut": {
            "train_end": "2023-06-30",
            "valid_end": "2023-09-30",
            "n_train": splits["train"],
            "n_valid": splits["valid"],
            "n_test": splits["test"],
        },
    }


def plot_board(payload: dict, path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), facecolor="#f7f5f1")
    fig.suptitle("Recsys features are causal time-OOF, not a hashing bag", fontsize=13, fontweight="bold")

    ax = axes[0]
    ax.bar(
        ["causal\n(past only)", "leaky\n(+ current review)"],
        [payload["causal"]["leak_rate"], payload["leaky_current_rate"]],
        color=["#2B6CB0", "#C53030"],
    )
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("current review in user_feat")
    ax.set_title("Does y_t's own text enter X_t?")

    ax = axes[1]
    ax.bar(
        ["current review\ntext vs label", "causal\nreview_cnt vs label"],
        [
            payload["shortcut"]["corr_current_review_text"],
            abs(payload["shortcut"]["corr_causal_review_cnt"]),
        ],
        color=["#C53030", "#2B6CB0"],
    )
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("|corr| with label_good")
    ax.set_title("Label shortcut if you keep this review")

    fig.text(
        0.5,
        0.02,
        "JingxiangQU/mmoe-multimodal-rec  ·  emit then update  ·  SplitByDate  ·  AGOD concat throws this away",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.08, 1, 0.92])
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def write_docs(p: dict, path: Path):
    md = f"""# 推荐系统的特征做法（因果时间 OOF）

看的是 Amazon 多模态推荐那条流水线：[`JingxiangQU/mmoe-multimodal-rec`](https://github.com/JingxiangQU/mmoe-multimodal-rec) 的 `data4moe_beam.py` / `data4model.py`，也就是 AGOD 吃的 WebDataset [`jingxiang11111/amazon_reviews_for_rec`](https://huggingface.co/datasets/jingxiang11111/amazon_reviews_for_rec)。

特征合同是 **时间**，不是折号。和上一张 Super Learner OOF 是同一件事：t 时刻的 X 不能由 y_t 算出来。模型头上的 MMoE 才是 `π(x)`——特征和门不要混。

---

## 他们实际写了什么

`CausalPosNegByUser`（按 user 排序）是整条协议：

```
for each user, sort reviews by timestamp:
    user_feat ← 只看 *过去* 的评论
    emit (user_feat, item_meta, image, label_good/best)
    *然后* 才把本条评论折进 running state
```

| 字段 | 从哪来 | 能不能看当前条 |
|---|---|---|
| `cat_hist` / `review_cnt` / `price_mean,std` | 过去评论的 Welford / 计数 | 否（先 emit 再 update） |
| `history`（最多 3 条 title+text） | 过去评论 | 否 |
| item 文本 | 类目 / 标题 / 价格 / features / description（目录元数据） | 目录，不是评分目标编码 |
| 图像 | 主图 → 14×14 patch | 目录 |
| `label_good` (≥4★) / `label_best` (5★) | **当前** 评分 | 标签，不当特征 |
| 当前评论 title/text | Enrich 里有，写 WebDataset 时丢掉 | **不准进 X** |

负样本：同一份因果 `user_feat`，随机抽没见过的 `parent_asin`，标签置 0。切分是 `SplitByDate`：train ≤ **2023-06-30** < valid ≤ **2023-09-30**，不是随机行切。

`build_user_text` / `build_item_text` 只是把结构化画像 **展成一段话**，给句级交叉注意力用。这是序列化，不是第二种泄漏。

---

## 和 OOF stacking 的同一张表

| | Super Learner / 上一张 | 这条推荐流水线 |
|---|---|---|
| 禁止 | `Z_{{i,m}}=f_m(x_i)` in-sample | 当前评论文本、未来的 user 统计 |
| 允许 | `f_m^{{(-i)}}(x_i)` | t 之前的评论 + 当时的目录 |
| 切分 | 折 / probe vs holdout | **事件时间** |
| 典型泄漏 | 噪声专家背 y | 用本条 review 预测本条星级（几乎是情感分类） |
| 不是这件事 | — | MMoE 门 `softmax(W·mean(experts))` 看的是 x |

本例玩具流：因果路径 current-review 泄漏率 **{p['causal']['leak_rate']:.2f}**；若把当前条折进 `user_feat`，泄漏率 **{p['leaky_current_rate']:.2f}**。当前评论文本 vs `label_good` 相关 **{p['shortcut']['corr_current_review_text']:.2f}**，因果 `review_cnt` 只有 **{p['shortcut']['corr_causal_review_cnt']:.2f}**。留下本条评论，任务就塌成读这句话。

---

## 他们还留的洞（特征侧，不是 MMoE）

- **负样本 item 池**是全局抽 10k `parent_asin`，不按 t 过滤上架时间。
- **item meta 是快照价格**，不是事件时刻的价。
- **5★ 降采样在 GroupByUser 之前**，`cat_hist` 不是真实历史，是抽过的历史。
- 过去评论的 *文本* 带情感（「did not work」），这是故意的用户画像，只要它是 **过去** 的就不是 y_t 泄漏。
- 没有 user_id / item_id embedding：冷启动友好，协同信号丢掉。

---

## AGOD 现在把合同扔了

`scripts/run_agod_amazon_modality_lr.py` 做的是 `text = item + "\\n" + user`，再 `HashingVectorizer`。user 画像和 item 目录变成 **一个 ngram 袋**。因果切分还在样本行上，但模型看不到「这是过去的用户、那是当前的商品」。PO/VIMP 还在当前窗上 in-sample 拟合——相对这条推荐协议，那是 leaky stacking。

MMoE 头：`π=softmax(W mean(expert_vecs))`，专家是 user 文本 / item 文本 / 图 / 交叉。**那是 MoE。** 特征流水线不是。上一张说的 stacking `π(t)` 仍然不是这个门。

---

## 可复现

```bash
PYTHONPATH=. python3 -m tests.test_agod_recsys_features
PYTHONPATH=. python3 scripts/run_agod_recsys_features.py
```

源码对照：`CausalPosNegByUser` 先 yield 再 `hist.append`；`SplitByDate`；`build_user_text` 只用 `user_feat`。
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(md)


def write_latex(p: dict, path: Path):
    tex = (
        "% Causal recsys features vs leaky current-review\n"
        "\\begin{table}[t]\\centering\n"
        "\\caption{Amazon recsys feature contract is causal time-OOF. "
        "Keeping the current review makes label prediction a sentiment task.}\n"
        "\\label{tab:agod-recsys-features}\n"
        "\\begin{tabular}{lcc}\\toprule\n"
        "protocol & current-review leak rate & corr(text, label) \\\\\n"
        "\\midrule\n"
        f"causal (emit then update) & {p['causal']['leak_rate']:.2f} & "
        f"{p['shortcut']['corr_causal_review_cnt']:.2f} \\\\\n"
        f"leaky (+ current review) & {p['leaky_current_rate']:.2f} & "
        f"{p['shortcut']['corr_current_review_text']:.2f} \\\\\n"
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    raw = toy_users()
    payload = run_all(raw)
    (OUT / "agod_recsys_features.json").write_text(json.dumps(payload, indent=2))
    plot_board(payload, OUT / "AGOD_Recsys_Features_Board.png")
    write_docs(payload, OUT / "README.md")
    write_latex(payload, OUT / "AGOD_recsys_features_tables_only.tex")
    write_docs(payload, DOCS / "AGOD_recsys_features.md")
    write_latex(payload, DOCS / "AGOD_recsys_features_tables_only.tex")
    shutil.copy2(OUT / "AGOD_Recsys_Features_Board.png", ART / "AGOD_Recsys_Features_Board.png")
    shutil.copy2(OUT / "README.md", ART / "README.md")
    shutil.copy2(
        OUT / "AGOD_Recsys_Features_Board.png",
        Path("/opt/cursor/artifacts/agod_recsys_features_board.png"),
    )
    print("causal leak_rate", payload["causal"]["leak_rate"], flush=True)
    print("leaky", payload["leaky_current_rate"], "shortcut", payload["shortcut"], flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
