# Labeled Hugging Face π-guided boards

Human feature-blocks → Stage-1 consensus π → Stage-2 (π-weighted RF, opinion pool, X+Z).
Y is the dataset label; W is the batch used for FSDS.

Compile the dashboard with:

```bash
cd docs/method && pdflatex -interaction=nonstopmode BlockAware_HF_dashboard.tex
```

## Datasets

### Adult census · income label · sex as batch

- Hub: `scikit-learn/adult-census-income`
- Label Y: income >50K
- Batch W: W = sex (Female vs Male); relationship dropped (Husband/Wife leak)
- n0=800, n1=800
- π: demography 0.480, work 0.290, hours 0.181, capital 0.049
- Observational AUC raw 0.855 vs π-RF 0.849

| setting | raw | π-RF | pool | X+Z | Δ π-RF | π-RF on GT |
|---|---:|---:|---:|---:|---:|---:|
| inject hours (α=0.40) | 0.950 | 0.987 | 0.988 | 0.987 | +0.037 | 0.947 |
| demography⊥ + inject hours | 0.945 | 0.987 | 0.986 | 0.983 | +0.042 | 0.984 |

### Yelp vs Amazon · sentiment label · marketplace as batch

- Hub: `fancyzhx/yelp_polarity vs fancyzhx/amazon_polarity`
- Label Y: binary sentiment
- Batch W: W = marketplace (Yelp vs Amazon); reviews truncated to 50 tokens
- n0=800, n1=800
- π: text 0.575, style 0.314, polarity 0.111
- Observational AUC raw 0.728 vs π-RF 0.725

| setting | raw | π-RF | pool | X+Z | Δ π-RF | π-RF on GT |
|---|---:|---:|---:|---:|---:|---:|
| inject polarity (α=0.70) | 0.821 | 0.817 | 0.827 | 0.834 | -0.004 | 0.602 |
| text⊥ + inject polarity | 0.794 | 0.810 | 0.823 | 0.816 | +0.016 | 0.691 |

### MultiNLI · entailment label · fiction vs telephone

- Hub: `nyu-mll/multi_nli validation_matched`
- Label Y: entailment vs rest
- Batch W: W = genre (fiction vs telephone)
- n0=800, n1=800
- π: premise 0.549, hypothesis 0.119, overlap 0.332
- Observational AUC raw 0.836 vs π-RF 0.818

| setting | raw | π-RF | pool | X+Z | Δ π-RF | π-RF on GT |
|---|---:|---:|---:|---:|---:|---:|
| inject overlap (α=0.80) | 0.924 | 0.969 | 0.971 | 0.977 | +0.045 | 0.836 |
| premise⊥ + inject overlap | 0.902 | 0.977 | 0.977 | 0.967 | +0.075 | 0.926 |

