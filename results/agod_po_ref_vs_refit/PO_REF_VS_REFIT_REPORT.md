# Reference PO vs re-fit PO-learner (OnlineRFPerm reject)

Blunt eval: gate opens → score three PO variants on the OOD batch →
compare (1) ranking quality vs oracle hardness, (2) next-MSE on sig-only.

## Downstream next-MSE (sig-only)

| dataset | n_sig | unif | ref_po | probe_po | **refit_po** | dre | best |
|---|---:|---:|---:|---:|---:|---:|---|
| `metro_interstate` | 6 | 8.602e+05 | 7.925e+05 | 7.839e+05 | **7.57e+05** | 1.642e+06 | `refit_po` |
| `beijing_pm25` | 5 | 2364 | 2261 | 2639 | **2678** | 3577 | `ref_po` |
| `stocks_AAPL` | 6 | 0.0007784 | 0.0008489 | 0.0008779 | **0.0009023** | 0.0008401 | `uniform` |
| `waymo_proxy` | 22 | 0.009992 | 0.01039 | 0.01112 | **0.01078** | 0.01154 | `uniform` |
| `stocks_MSFT` | 8 | 0.0006164 | 0.0006494 | 0.0007019 | **0.0007034** | 0.0005892 | `dre` |
| `stocks_IWM` | 7 | 0.0005025 | 0.0005518 | 0.0005562 | **0.0005362** | 0.0005117 | `uniform` |

**Wins (sig-only):** `uniform`=3, `ref_po`=1, `probe_po`=0, `refit_po`=1, `dre`=1
**refit_po < ref_po:** `2/6`
**refit_po < probe_po:** `3/6`

## PO ranking quality on reject batches

| dataset | spearman ref | probe | **refit** | topk ref | probe | **refit** |
|---|---:|---:|---:|---:|---:|---:|
| `metro_interstate` | 0.302 | 0.316 | **0.320** | 0.286 | 0.364 | **0.350** |
| `beijing_pm25` | 0.084 | 0.393 | **0.399** | 0.350 | 0.520 | **0.470** |
| `stocks_AAPL` | 0.287 | 0.670 | **0.665** | 0.521 | 0.707 | **0.693** |
| `waymo_proxy` | 0.115 | 0.549 | **0.529** | 0.259 | 0.595 | **0.566** |
| `stocks_MSFT` | 0.378 | 0.753 | **0.721** | 0.494 | 0.700 | **0.633** |
| `stocks_IWM` | 0.608 | 0.725 | **0.689** | 0.650 | 0.714 | **0.714** |

### How to read this

- **PO itself** is just absolute residual risk — intentionally blunt.
- **Spearman / topk**: does the PO score pick the same hard rows as an
  oracle residual? Higher ⇒ better OOD instance ranking.
- **sig-only MSE**: after IPTW with √PO, does next-batch error drop on
  the batches where the gate actually fired?
- Expectation: `refit_po` ≥ `probe_po` ≥ `ref_po` on ranking; MSE lift is
  softer and may still lose to uniform on some packs.
