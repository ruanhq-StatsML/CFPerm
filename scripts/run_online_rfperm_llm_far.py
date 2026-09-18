#!/usr/bin/env python3
"""Run onlinePermOOB_with_LLM on the tables we have. Print a FAR latex table."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT))

from online_rfperm_with_llm import (  # noqa: E402
    DETECTORS,
    compute_vimp,
    far_table_latex,
    make_rf_adapter,
    make_tabpfn_adapter,
    onlinePermOOB_with_LLM,
    run_far_study,
    split_by_vimp,
)

OUT = ROOT / "results" / "online_rfperm_llm"


def _frame_to_xy(df, y_raw):
    import pandas as pd

    Y = np.asarray(y_raw).ravel()
    if Y.dtype == object or str(Y.dtype).startswith("str") or str(Y.dtype) == "category":
        ys = pd.Series(Y.astype(str))
        pos = {"1", "UP", "True", "true", "Y", "yes"}
        if set(ys.unique()) <= {"0", "1"} or any(v in pos for v in ys.unique()):
            Y = ys.isin(sorted(pos)).astype(int).to_numpy()
        else:
            top = ys.value_counts().index[0]
            Y = (ys == top).astype(int).to_numpy()
    else:
        u = np.unique(Y)
        Y = (Y == (2 if 2 in u else u[np.argmax([np.sum(Y == v) for v in u])])).astype(int)
    blocks = []
    for c in df.columns:
        s = df[c]
        if str(s.dtype) in {"object", "category"} or s.dtype == object:
            codes, _ = pd.factorize(s.astype(str), sort=False)
            blocks.append(codes.astype(float))
        else:
            blocks.append(np.asarray(s, dtype=float))
    return np.column_stack(blocks), Y.astype(float)


def _cap(X, Y, n=4000):
    n = min(int(n), len(Y))
    return np.column_stack([X[:n], Y[:n]])


def load_tables():
    """Whatever is already on the board / fetchable. Skip what is missing."""
    out = {}
    try:
        from sklearn.datasets import fetch_openml

        bunch = fetch_openml("eeg-eye-state", version=1, as_frame=True, parser="auto")
        X, Y = _frame_to_xy(bunch.data.copy(), bunch.target)
        out["eeg"] = _cap(X, Y, 4000)
    except Exception as exc:
        print("skip eeg:", exc)
    try:
        from sklearn.datasets import fetch_openml

        bunch = fetch_openml("electricity", version=1, as_frame=True, parser="auto")
        df = bunch.data.copy()
        drop = [c for c in df.columns if str(c).lower() == "date"]
        X, Y = _frame_to_xy(df.drop(columns=drop), bunch.target)
        out["electricity"] = _cap(X, Y, 4000)
    except Exception as exc:
        print("skip electricity:", exc)
    try:
        from sklearn.datasets import fetch_openml

        bunch = fetch_openml("bank-marketing", version=1, as_frame=True, parser="auto")
        X, Y = _frame_to_xy(bunch.data.copy(), bunch.target)
        out["bankmarketing"] = _cap(X, Y, 4000)
    except Exception as exc:
        print("skip bankmarketing:", exc)
    try:
        from sklearn.datasets import fetch_covtype

        bunch = fetch_covtype()
        X = np.asarray(bunch.data, dtype=float)
        Y = (np.asarray(bunch.target) == 2).astype(float)
        out["covertype"] = _cap(X, Y, 4000)
    except Exception as exc:
        print("skip covertype:", exc)
    return out


def run_one(name, df, ref_batch_size=800, batch_size=50, seed=2026):
    llm = make_tabpfn_adapter(seed=seed)
    dl = make_rf_adapter()
    X_ref = df[:ref_batch_size, :-1]
    Y_ref = df[:ref_batch_size, -1]
    n_high = min(4, X_ref.shape[1] // 2)
    high_idx, low_idx = split_by_vimp(compute_vimp(X_ref, Y_ref, seed=seed), n_high=n_high)
    rec = onlinePermOOB_with_LLM(
        df, llm_adapter=llm, dl_adapter=dl,
        llm_cols=low_idx, dl_cols=high_idx,
        w_llm=0.5, w_dl=0.5,
        ref_batch_size=ref_batch_size, batch_size=batch_size,
        burnin=5, seed=seed,
    )
    rec["name"] = name
    rec["n"] = int(len(df))
    rec["p"] = int(df.shape[1] - 1)
    rec["ref_batch_size"] = int(ref_batch_size)
    rec["batch_size"] = int(batch_size)
    return rec


def dataset_table_latex(rows):
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{onlinePermOOB\_with\_LLM on the tables we have. Last column is $Y$. "
        r"high-VIMP RF + low-VIMP TabPFN-style $k$NN. Entries are first rejection batch ($-1$ = none).}",
        r"\label{tab:data}",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"table & $n$ & $p$ & ADDIS & SAFFRON & PH & ADWIN \\",
        r"\midrule",
    ]
    for rec in rows:
        lines.append(
            f"{rec['name']} & {rec['n']} & {rec['p']} & "
            f"{rec['addis_1']} & {rec['saffron_1']} & {rec['PH_1']} & {rec['ADWIN_1']} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    return "\n".join(lines)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    study = run_far_study(n_reps=12, n=2000, p=10, ref_batch_size=400, batch_size=40, seed0=2026)
    far_tex = far_table_latex(study)
    (OUT / "far.json").write_text(json.dumps(study, indent=2), encoding="utf-8")
    (OUT / "far_table.tex").write_text(far_tex, encoding="utf-8")

    tables = load_tables()
    rows = []
    for name, df in tables.items():
        ref = min(800, max(200, len(df) // 3))
        batch = 50 if len(df) - ref >= 150 else max(20, (len(df) - ref) // 6)
        print("running", name, df.shape, "ref", ref, "batch", batch)
        rec = run_one(name, df, ref_batch_size=ref, batch_size=batch)
        slim = {k: rec[k] for k in (
            "name", "n", "p", "ref_batch_size", "batch_size", "mse_ref",
            "addis_1", "saffron_1", "fix_1", "hop_1", "PH_1", "EWMA_1",
            "CUSUM_1", "DDM_1", "ADWIN_1", "martingale_1",
        )}
        rows.append(slim)
        print(" ", slim)

    data_tex = dataset_table_latex(rows)
    (OUT / "data_table.tex").write_text(data_tex, encoding="utf-8")
    (OUT / "data.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")

    note = (
        r"""% Compile: pdflatex docs/online_rfperm_llm_far.tex
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{microtype}
\title{onlinePermOOB with LLM routing}
\author{}
\date{}
\begin{document}
\maketitle
\thispagestyle{empty}

\paragraph{Setup.}
\texttt{onlinePermOOB\_with\_LLM(df, \ldots)}. Last column of \texttt{df} is $Y$, never a feature.
high-VIMP columns go to a frozen RF, low-VIMP columns go to a TabPFN-style $k$NN.
Fit once on the reference window. Trail batches score $T_t=\mathrm{MSE}_t-E_{\mathrm{ref}}$.
ADDIS is the primary mark.

"""
        + far_tex
        + "\n"
        + data_tex
        + r"""
\end{document}
"""
    )
    (ROOT / "docs" / "online_rfperm_llm_far.tex").write_text(note, encoding="utf-8")

    md = ["# onlinePermOOB_with_LLM", "", "## FAR", ""]
    md.append("| method | stationary | random-noise |")
    md.append("|---|---|---|")
    for name in DETECTORS:
        s = study["stationary"]["far"][name]
        n = study["random_noise"]["far"][name]
        md.append(f"| {name} | {100*s:.1f}% | {100*n:.1f}% |")
    md += ["", "## tables we have", "", "| table | n | p | ADDIS | SAFFRON | PH | ADWIN |", "|---|---|---|---|---|---|---|"]
    for rec in rows:
        md.append(
            f"| {rec['name']} | {rec['n']} | {rec['p']} | {rec['addis_1']} | "
            f"{rec['saffron_1']} | {rec['PH_1']} | {rec['ADWIN_1']} |"
        )
    (OUT / "FAR.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(far_tex)
    print(data_tex)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
