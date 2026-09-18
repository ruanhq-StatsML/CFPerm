#!/usr/bin/env python3
"""FAR of OnlineRFPerm-with-LLM routing on stationary and random-noise streams."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(ROOT))

from online_rfperm_with_llm import DETECTOR_ORDER, far_table_latex, run_far_study  # noqa: E402

OUT = ROOT / "results" / "online_rfperm_llm"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    study = run_far_study(
        n_reps=25,
        n_ref=400,
        n_new=40,
        n_batches=20,
        p=10,
        n_high=4,
        seed0=2026,
    )
    (OUT / "far.json").write_text(json.dumps(study, indent=2), encoding="utf-8")
    tex = far_table_latex(study)
    (OUT / "far_table.tex").write_text(tex, encoding="utf-8")
    note = (
        r"""% Compile: pdflatex docs/online_rfperm_llm_far.tex
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{microtype}
\title{OnlineRFPerm with LLM routing --- false-alarm rate}
\author{}
\date{}
\begin{document}
\maketitle
\thispagestyle{empty}

\paragraph{Setup.}
One row is one observation. $Y$ is the numeric outcome, never a feature.
$X$ is routed: high-VIMP columns go to a frozen RF, low-VIMP columns go to a
TabPFN-style context $k$NN. The blend is fit once on $D_{\mathrm{ref}}$.
Each new batch is $T_t=\mathrm{MSE}_t-E_{\mathrm{ref}}$. ADDIS is the primary
OnlineRFPerm mark. The remaining rows are the same MSE stream, other detectors.

Stationary: $Y=X^\top w+\varepsilon$, no onset. Random-noise: $X$ and $Y$ independent $N(0,1)$, no onset.
FAR is the share of replications that reject at least once. There is no labeled change, so any ring is a false alarm.

"""
        + tex
        + r"""

\end{document}
"""
    )
    (ROOT / "docs" / "online_rfperm_llm_far.tex").write_text(note, encoding="utf-8")
    lines = [
        "# OnlineRFPerm + LLM routing --- FAR",
        "",
        f"n_reps={study['stationary']['n_reps']}, n_batches={study['stationary']['n_batches']}.",
        "Routing: high-VIMP RF + low-VIMP TabPFN-style kNN.",
        "",
        "| method | stationary FAR | random-noise FAR |",
        "|---|---|---|",
    ]
    for name in DETECTOR_ORDER:
        s = study["stationary"]["far"][name]
        n = study["random_noise"]["far"][name]
        lines.append(f"| {name} | {100*s:.1f}% | {100*n:.1f}% |")
    (OUT / "FAR.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(tex)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
