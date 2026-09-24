#!/usr/bin/env python3
"""NDCG@k frozen probe on OnlineRFPerm. IR, not quantile regression."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from online_rfperm_rank import make_rank_df, onlinePermOOB_rank  # noqa: E402


def _rate(key, runs):
    n = max(len(runs), 1)
    return float(sum(int(r.get(key, -1) not in (-1, None)) for r in runs)) / n


def study(n_reps=6, seed0=4):
    out = {}
    for kind in ("quiet", "rank_flip", "covariate"):
        runs = []
        for r in range(int(n_reps)):
            df, g = make_rank_df(
                n_queries=50, slate=8, kind=kind, onset_q=25, seed=int(seed0) + r
            )
            rec = onlinePermOOB_rank(
                df, g, k=4, ref_batch_size=160, batch_size=40, seed=int(seed0) + r
            )
            runs.append(rec)
        out[kind] = {
            "n": int(n_reps),
            "far_ndcg": _rate("ndcg_1", runs),
            "far_mse": _rate("mse_1", runs),
            "far_x": _rate("x_1", runs),
        }
    return out


def tex_table(eff) -> str:
    lines = [
        r"\begin{tabular}{@{}lccc@{}}",
        r"\toprule",
        r"scene & NDCG hop & MSE hop & MMD hop \\",
        r"\midrule",
    ]
    for k in ("quiet", "rank_flip", "covariate"):
        r = eff[k]
        lines.append(
            f"{k} & {r['far_ndcg']:.2f} & {r['far_mse']:.2f} & {r['far_x']:.2f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def write_docs(eff):
    table = tex_table(eff)
    body = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs}",
            r"\begin{document}",
            r"\paragraph{IR probe, not quantile regression.}",
            r"OnlineRFPerm's MSE is a pointwise rent check. Ranking systems pay rent in NDCG.",
            r"Freeze a scorer on $D_{\mathrm{ref}}$. Each batch is slates. Loss $=1-\mathrm{NDCG}@k$",
            r"vs the hold-out slice of $D_{\mathrm{ref}}$. $T$ is the batch label. $Y$ is graded",
            r"relevance, last column, never a feature. Slate ids are not $X$.",
            r"\paragraph{Read.}",
            r"A rank-rule flip should ring NDCG. An unused-column walk should ring MMD and",
            r"leave ranking quiet. Quiet FAR on NDCG is the IR analogue of the hop twitch;",
            r"ADDIS stays the primary mark if you put NDCG-loss on the $p$-stream.",
            r"\vspace{0.8em}",
            table,
            r"\end{document}",
        ]
    )
    (ROOT / "docs").mkdir(parents=True, exist_ok=True)
    (ROOT / "docs" / "online_rfperm_rank.tex").write_text(body, encoding="utf-8")
    res = ROOT / "results" / "online_rfperm_rank"
    res.mkdir(parents=True, exist_ok=True)
    (res / "board.tex").write_text(table, encoding="utf-8")
    (res / "effect.json").write_text(json.dumps(eff, indent=2), encoding="utf-8")


if __name__ == "__main__":
    print("[rank] 6 reps × quiet / rank_flip / covariate ...")
    eff = study()
    for k, v in eff.items():
        print(f"  {k:10s}  ndcg={v['far_ndcg']:.2f}  mse={v['far_mse']:.2f}  x={v['far_x']:.2f}")
    write_docs(eff)
    print("wrote docs/online_rfperm_rank.tex")
