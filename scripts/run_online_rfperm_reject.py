#!/usr/bin/env python3
"""Reject-inference readout and a quantile probe on OnlineRFPerm."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from online_rfperm_reject import (  # noqa: E402
    make_reject_df,
    make_tail_df,
    onlinePermOOB_quantile,
    onlinePermOOB_reject,
)


def _rate(flags_key, runs):
    n = max(len(runs), 1)
    return float(sum(int(r.get(flags_key, -1) not in (-1, None)) for r in runs)) / n


def reject_study(n_reps=6, seed0=7):
    kinds = ("quiet", "select", "concept", "covariate")
    out = {}
    for kind in kinds:
        runs = []
        for r in range(int(n_reps)):
            df = make_reject_df(
                n=400, p=6, kind=kind, onset=200, seed=int(seed0) + r
            )
            rec = onlinePermOOB_reject(
                df, ref_batch_size=140, batch_size=40, seed=int(seed0) + r
            )
            runs.append(rec)
        out[kind] = {
            "n": int(n_reps),
            "far_cc": _rate("cc_1", runs),
            "far_ips": _rate("ips_1", runs),
            "far_sel": _rate("sel_1", runs),
            "far_x": _rate("x_1", runs),
        }
    return out


def quantile_study(n_reps=6, seed0=11):
    quiet, tail = [], []
    for r in range(int(n_reps)):
        dq = make_tail_df(n=400, p=6, onset=10_000, seed=int(seed0) + r)
        dt = make_tail_df(n=400, p=6, onset=200, seed=int(seed0) + 50 + r)
        quiet.append(
            onlinePermOOB_quantile(dq, tau=0.9, ref_batch_size=140, batch_size=40, seed=r)
        )
        tail.append(
            onlinePermOOB_quantile(dt, tau=0.9, ref_batch_size=140, batch_size=40, seed=r)
        )
    return {
        "quiet": {"n": n_reps, "far_mse": _rate("mse_1", quiet), "far_pinball": _rate("pinball_1", quiet)},
        "tail": {"n": n_reps, "far_mse": _rate("mse_1", tail), "far_pinball": _rate("pinball_1", tail)},
    }


def tex_reject(rej, qtl) -> str:
    lines = [
        r"\begin{tabular}{@{}lcccc@{}}",
        r"\toprule",
        r"scene & complete-case hop & IPS hop & $P(S\mid X)$ hop & MMD hop \\",
        r"\midrule",
    ]
    for k in ("quiet", "select", "concept", "covariate"):
        r = rej[k]
        lines.append(
            f"{k} & {r['far_cc']:.2f} & {r['far_ips']:.2f} & {r['far_sel']:.2f} & {r['far_x']:.2f} \\\\"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\vspace{1.1em}",
        r"\begin{tabular}{@{}lcc@{}}",
        r"\toprule",
        r"scene & MSE hop & pinball $0.9$ hop \\",
        r"\midrule",
        f"quiet tail & {qtl['quiet']['far_mse']:.2f} & {qtl['quiet']['far_pinball']:.2f} \\\\",
        f"upper-tail hop & {qtl['tail']['far_mse']:.2f} & {qtl['tail']['far_pinball']:.2f} \\\\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    return "\n".join(lines)


def write_docs(rej, qtl):
    table = tex_reject(rej, qtl)
    body = "\n".join(
        [
            r"\documentclass[11pt]{article}",
            r"\usepackage[margin=1in]{geometry}",
            r"\usepackage{booktabs,amsmath}",
            r"\begin{document}",
            r"\paragraph{Not rejection inference.}",
            r"OnlineRFPerm asks whether a frozen probe still pays rent on the new batch.",
            r"$T$ is the batch label. Complete-case MSE on shown rows is not reject inference:",
            r"we do not impute $Y$ on $S=0$. Non-finite last-column $Y$ means unshown.",
            r"A policy hop $P(S=1\mid X)$ can move complete-case MSE while $P(Y\mid X)$ holds.",
            r"IPS-weighted MSE on $S=1$ is a reweight under selection on observables,",
            r"not a unique CATE, not Shapley.",
            r"\paragraph{Read.}",
            r"Select rings $P(S\mid X)$ (labeled-slice MMD $+$ shown rate) while full-X MMD stays at the quiet floor.",
            r"Concept rings the $Y$ probes with that same quiet floor on selection and MMD.",
            r"Covariate walks $X_0$, which also feeds $S$, so both MMD and selection hop --- that is the DGP.",
            r"Complete-case MSE FAR on quiet is still twitchy; ADDIS stays the primary mark on the $Y$ stream.",
            r"\paragraph{A second probe.}",
            r"MSE is $L^2$. Pinball at $\tau=0.9$ is a frozen upper-tail rent check.",
            r"A tail hop can ring pinball first.",
            r"\vspace{0.8em}",
            table,
            r"\end{document}",
        ]
    )
    docs = ROOT / "docs"
    docs.mkdir(parents=True, exist_ok=True)
    (docs / "online_rfperm_reject.tex").write_text(body, encoding="utf-8")
    res = ROOT / "results" / "online_rfperm_reject"
    res.mkdir(parents=True, exist_ok=True)
    (res / "board.tex").write_text(table, encoding="utf-8")
    (res / "effect.json").write_text(
        json.dumps({"reject": rej, "quantile": qtl}, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    print("[reject] 6 reps × quiet/select/concept/covariate ...")
    rej = reject_study()
    for k, v in rej.items():
        print(
            f"  {k:10s}  cc={v['far_cc']:.2f}  ips={v['far_ips']:.2f}  "
            f"sel={v['far_sel']:.2f}  x={v['far_x']:.2f}"
        )
    print("[quantile] 6 reps × quiet / upper-tail ...")
    qtl = quantile_study()
    for k, v in qtl.items():
        print(f"  {k:10s}  mse={v['far_mse']:.2f}  pinball={v['far_pinball']:.2f}")
    write_docs(rej, qtl)
    print("wrote docs/online_rfperm_reject.tex")
