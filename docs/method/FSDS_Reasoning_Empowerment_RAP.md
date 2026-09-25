# FSDS Reasoning Empowerment — RAP LaTeX prototype

**Main file:** [`FSDS_Reasoning_Empowerment_RAP.tex`](FSDS_Reasoning_Empowerment_RAP.tex)

Standalone writeup: how FSDS **marks expandable nodes** to guide the
next step inside Reasoning-via-Planning (RAP) / ToT.

Includes:

1. RAP bottleneck = next-step policy (expand / prune / reorder / backtrack)
2. Formal $X_i$, $Y_i=\mathbf{1}\{i\in\pi^\star\}$, $\mathrm{imp}(i)=\widehat{P}(Y_i=1\mid X_i)$
3. Overlap router + four FSDS components
4. Guide algorithm (expandable-node marker)
5. NetValue / ROI ledger
6. Justification (why FSDS, why $Y$, why RAP interface)
7. Scorecard numbers (success $13.3\%\to 93.3\%$, $\Delta\mathrm{NetValue}\approx +4.45$, ROI $\approx 29.7$)
8. Plug-in recipe for a real RAP loop

```bash
pdflatex docs/method/FSDS_Reasoning_Empowerment_RAP.tex
```

Code: `agod/fsds_reason_guide.py` · Summary: `docs/summaries/FSDS_Reasoning_Empowerment_NextStep.md` · PR #87
