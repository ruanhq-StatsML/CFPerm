# LLM / FSDS Landing Reports

Bilingual LaTeX briefings that map **RFPerm / OnlineRFPerm / PO-learner / FSDS / multi-layer graph methods** to ready LLM landing scenarios (inference monitoring, alignment gates, hallucination regimes, feature-bank quality, multimodal scheduling).

| File | Language |
|------|----------|
| [`LLM_FSDS_Landing_Report_CN.tex`](LLM_FSDS_Landing_Report_CN.tex) | Chinese |
| [`LLM_FSDS_Landing_Report_EN.tex`](LLM_FSDS_Landing_Report_EN.tex) | English |

## Build

```bash
cd docs/reports
xelatex LLM_FSDS_Landing_Report_CN.tex   # needs xeCJK + WenQuanYi Micro Hei
pdflatex LLM_FSDS_Landing_Report_EN.tex
```

## Sources

- [RFPerm](https://github.com/ruanhq-StatsML/RFPerm)
- [OnlineRFPerm](https://github.com/ruanhq-StatsML/OnlineRFPerm)
- [PO-learner / Causal Objective Permutation Test](https://github.com/ruanhq-StatsML/Causal_Objective_Permutation_Test)
- [FSDS](https://github.com/ruanhq-StatsML/FSDS-FeatureSelection_for_DistributionShift) (manuscript + Sep14 multi-layer metrics PDF)
- Local: `docs/po_risk/PO_risk_reweight_method.tex`, `docs/agod/AGOD_online_rfperm.md`
