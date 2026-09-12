# Amazon Reviews 2023: 9-batch TSS vs LR clocks

Same next-epoch TSS map as the multimodal probe, on a scalar text head. Nine McAuley product categories are consecutive batches (gift cards → digital music → beauty → software → subscription boxes → instruments → magazines → handmade → fashion). `X` is TF-IDF frozen on batch 0; `Y` is the 1–5 star rating. Linear MSE, predict-then-update.

Do not read this as retrieval or ranking. Category hops are a domain stream, not a timestamp stream.

## Estimators

- \(c\): \(|d|\) of the TF-IDF coordinate-mean under consecutive categories.
- \(\delta\): two-fold excess ridge MSE after mean-aligning \(X\).
- \(\eta^\star=\eta_0(1+\beta\hat\delta/s)/(1+\lambda\hat c\, s)\), \(s=\sqrt{n_{\mathrm{iter}}}\). Quiet holds last \(\eta\). Train batch \(t\) with \(\eta^{(t)}\), then set \(\eta^{(t+1)}\).

On this extract (6 seeds, 240 reviews/category): mean \(\hat c=0.309\), mean \(\hat\delta=0.387\). Both channels are live. TSS net-shrinks \(\bar\eta\) to \(0.039\) (inv-\(c\) to \(0.034\)) because \(\lambda c s\) dominates \(\beta\delta/s\).

## What to compare

Online MSE, BWT (MSE on batch 0 after the stream), and regret vs the hindsight-best constant \(\eta\) inside the TSS box \([0.04, 2.5\eta_0]\). Unconstrained larger constants still cut online MSE: the linear probe is underfitting, which is not a scheduler ranking.

| method | online MSE | regret | BWT MSE | \(\bar\eta\) |
| --- | --- | --- | --- | --- |
| plateau | 1.597 (0.055) | 0.212 (0.146) | 0.879 (0.121) | 0.058 |
| cosine | 1.604 (0.056) | 0.275 (0.126) | 0.955 (0.125) | 0.044 |
| constant | 1.607 (0.059) | 0.296 (0.108) | 0.903 (0.125) | 0.100 |
| Polyak | 1.632 (0.061) | 0.492 (0.108) | 0.841 (0.138) | 0.109 |
| TSS | 1.634 (0.059) | 0.515 (0.173) | 0.818 (0.103) | 0.039 |
| inv-\(c\) | 1.644 (0.060) | 0.594 (0.161) | 0.784 (0.108) | 0.034 |

Hindsight-best constant in the box is \(\eta=0.25\) on every seed. Plateau has the lowest online MSE and regret (it only halves after two stalls). TSS tracks inv-\(c\): shrinking under \(\hat c\) helps BWT (0.818 / 0.784 vs constant 0.903) and costs a bit of online MSE. Polyak raises \(\eta\) when the loss spikes, the wrong sign for a covariate-heavy category hop.

Paired Wilcoxon (6 seeds, all signs agree): TSS is worse than plateau / cosine / constant on online MSE and regret (\(p=0.03\)) and better on BWT (\(p=0.03\)). TSS vs inv-\(c\): lower regret, higher BWT MSE (\(p=0.03\)). Do not oversell TSS vs plateau on this stream.

## Run

```
python3 scripts/run_amazon_continuous_batches.py
python3 scripts/run_amazon_continuous_batches.py --synthetic
python3 scripts/run_amazon_continuous_batches.py --quick
```

Heads of the McAuley JSONL files are range-fetched at runtime (`data/amazon/`, gitignored). Method writeup: `docs/method/Amazon_continuous_batches_note.tex`. Table: `Amazon_continuous_batches.tex`.
