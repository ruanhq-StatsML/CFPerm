# Amazon Reviews 2023: 9-batch TSS vs LR clocks

Same next-epoch TSS map as the multimodal probe, on a scalar text head. Nine McAuley product categories are consecutive batches. `X` is TF-IDF frozen on batch 0; `Y` is the 1–5 star rating. Linear MSE, predict-then-update.

The **batch-relationship heatmap** is the same \(\mathrm{cosine}(B_i,B_j)\) object as the MSR-VTT modality boards (one text head; batches = categories). Sparse TF-IDF pairwise cosine is near 0, so the board is cosine of category-mean vectors. TSS \(\hat c\) is the hop the trainer walks: \(1-\cos(\bar x_{t-1},\bar x_t)\). Gift cards is the outlier (cosine \(\approx 0.50\)–\(0.61\) with every other category); music/magazines and beauty/instruments sit close (\(\approx 0.91\)–\(0.93\)). Mean heatmap \(\hat c=0.263\), mean \(\hat\delta=0.387\).

Do not read this as retrieval or ranking. Category hops are a domain stream, not a timestamp stream.

## What to compare

Online MSE, BWT (MSE on batch 0 after the stream), and regret vs the hindsight-best constant \(\eta\) inside the TSS box \([0.04, 2.5\eta_0]\). Unconstrained larger constants still cut online MSE: the linear probe is underfitting, which is not a scheduler ranking.

| method | online MSE | regret | BWT MSE | \(\bar\eta\) |
| --- | --- | --- | --- | --- |
| plateau | 1.597 (0.055) | 0.212 (0.146) | 0.879 (0.121) | 0.058 |
| cosine | 1.604 (0.056) | 0.275 (0.126) | 0.955 (0.125) | 0.044 |
| constant | 1.607 (0.059) | 0.296 (0.108) | 0.903 (0.125) | 0.100 |
| TSS | 1.631 (0.059) | 0.485 (0.150) | 0.883 (0.099) | 0.046 |
| Polyak | 1.632 (0.061) | 0.492 (0.108) | 0.841 (0.138) | 0.109 |
| inv-\(c\) | 1.633 (0.060) | 0.506 (0.142) | 0.849 (0.112) | 0.039 |

Hindsight-best constant in the box is \(\eta=0.25\) on every seed. Plateau still has the lowest online MSE and regret. Heatmap \(\hat c\) makes TSS track the Gift→Music / Software→SubBox hops (shrink) and hold when Music–Beauty–Software sit close. Do not oversell TSS vs plateau on this stream.

## Run

```
python3 scripts/run_amazon_continuous_batches.py
python3 scripts/run_amazon_continuous_batches.py --synthetic
```

Lead figure: `amazon_batch_relationship_heatmap.png`. Comparison: `amazon_vs_schedulers.png`. Method: `docs/method/Amazon_continuous_batches_note.tex`.
