# Post-hoc localization on the audit tables

Directly on the csv: pull subset indices, **look at the mean**, then call `MMD()` and `po_risk()`.
HH chosen is not Y.

```bash
PYTHONPATH=. python3 scripts/run_posthoc_on_dataset.py
```

| Dataset | subsets from | n | sig pairs |
|---|---|---:|---:|
| HH two-stream | column T | 2400 | 1 |
| HH multi-step hops | step clipped to 0/1/2+ | 1200 | 3 |
| HH helpful policy hop | batch >= 4 | 1200 | 1 |
| BeaverTails vs ToxicChat | column T | 2400 | 1 |

Look at the mean first (already computed). Then pairwise MMD / PO-risk.

## Conditional mean (already computed — look at this)

| Stream | subset | n | mean Y | mean top x |
|---|---|---:|---:|---:|
| HH two-stream | T0 | 1200 | 0.491 | 1.67 |
| HH two-stream | T1 | 1200 | 0.491 | 1.42 |
| HH two-stream | Q0 | 623 | 0.520 | 0.339 |
| HH two-stream | Q1 | 577 | 0.390 | 0.816 |
| HH two-stream | Q2 | 600 | 0.498 | 1.49 |
| HH two-stream | Q3 | 600 | 0.550 | 3.54 |
| HH multi-step hops | T0 | 380 | 0.395 | 0.158 |
| HH multi-step hops | T1 | 379 | 0.417 | 0.175 |
| HH multi-step hops | T2 | 441 | 0.370 | 0.152 |
| HH multi-step hops | Q0 | 306 | 0.000 | 0.0345 |
| HH multi-step hops | Q1 | 296 | 0.000 | 0.0833 |
| HH multi-step hops | Q2 | 299 | 0.579 | 0.158 |
| HH multi-step hops | Q3 | 299 | 0.997 | 0.372 |
| HH helpful policy hop | T0 | 320 | 0.519 | 1.65 |
| HH helpful policy hop | T1 | 880 | 0.515 | 1.68 |
| HH helpful policy hop | Q0 | 303 | 0.624 | 0.354 |
| HH helpful policy hop | Q1 | 309 | 0.589 | 0.897 |
| HH helpful policy hop | Q2 | 296 | 0.436 | 1.68 |
| HH helpful policy hop | Q3 | 292 | 0.408 | 3.85 |
| BeaverTails vs ToxicChat | T0 | 1200 | 0.427 | 0.299 |
| BeaverTails vs ToxicChat | T1 | 1200 | 0.780 | 0.662 |
| BeaverTails vs ToxicChat | Q0 | 612 | 0.632 | 0.0938 |
| BeaverTails vs ToxicChat | Q1 | 593 | 0.553 | 0.268 |
| BeaverTails vs ToxicChat | Q2 | 599 | 0.487 | 0.455 |
| BeaverTails vs ToxicChat | Q3 | 596 | 0.740 | 1.11 |

Read the mean first. Pairwise MMD and PO-risk are the significance next to it.

## Pairwise subset MMD / PO-risk

| Dataset | pair | n | MMD | MMD p | PO-risk | PO p | mean Y |
|---|---|---|---:|---:|---:|---:|---|
| HH two-stream | T0 vs T1 | 1200/1200 | 0.00492 | 0.0385 | 0.013 | 0.0385 | 0.491/0.491 |
| HH multi-step hops | T0 vs T1 | 380/379 | 0.0563 | 0.0385 | 2.57e-09 | 1 | 0.395/0.417 |
| HH multi-step hops | T0 vs T2 | 380/441 | 0.271 | 0.0385 | 0.000169 | 1 | 0.395/0.370 |
| HH multi-step hops | T1 vs T2 | 379/441 | 0.121 | 0.0385 | 0.000169 | 1 | 0.417/0.370 |
| HH helpful policy hop | T0 vs T1 | 320/880 | 0 | 1 | 0.00649 | 0.0385 | 0.519/0.515 |
| BeaverTails vs ToxicChat | T0 vs T1 | 1200/1200 | 0.279 | 0.0385 | 0.00546 | 0.0385 | 0.427/0.780 |

Quartile splits of `x_n_toks` are in the mean table (Q0–Q3). Nothing else.

