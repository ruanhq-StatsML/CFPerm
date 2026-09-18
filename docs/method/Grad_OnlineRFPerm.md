# Grad-OnlineRFPerm — unified LaTeX writeup

**Main file:** [`Grad_OnlineRFPerm.tex`](Grad_OnlineRFPerm.tex)

Standalone formula + experiment note (not AGOD). Includes:

1. Single-stream口径: $g_t=\|\nabla_{\theta_U} L\|_2$
2. OnlineRFPerm $T/p$/FDR
3. Diagnostic layer shares
4. MVP lead table (5×5)
5. Extras: null/grace, α sweep, freeze closed-loop, extra packs

```bash
pdflatex docs/method/Grad_OnlineRFPerm.tex
```

PR: https://github.com/ruanhq-StatsML/CFPerm/pull/64
