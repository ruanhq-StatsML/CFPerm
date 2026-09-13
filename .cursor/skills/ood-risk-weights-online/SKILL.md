---
name: ood-risk-weights-online
description: >-
  Quantify and compare online learning methods that use OOD risk / proximity /
  gradient-alignment scores as sample or domain weights (hop IW, DGA, Mahalanobis,
  attr adapter). Use when the user asks to evaluate OOD-as-weights, online MSE
  and regret on streaming boards (Amazon, Diabetes, interstate/stock surrogates),
  or to extend the OOD-weight online eval harness.
---

# OOD-risk as weights — online quantification

## When to use

- User wants **OOD degree / proximity / gradient alignment as training weights**
- Compare **hop IW vs DGA vs Mahalanobis vs uniform** on streaming boards
- Report **online MSE, cumulative regret vs uniform, BWT, oracle gap** + plots
- Extend boards (new dataset) without inventing a new loss

## Hard rules

1. **Do not invent a new loss.** Only change sample/domain weights into existing Ridge / adapter.
2. **Do not claim Waymo / Interstate / raw stock** unless those feeds are actually in the repo. Use `interstate` / `stock` **surrogates** or say they are missing.
3. Amazon champion stays **`hop`** unless a board clearly beats it; report numbers honestly.
4. Regret comparator for closed-form Ridge = **uniform `ridge_past`** on the same seed (not SGD η-hindsight unless using the Amazon continuous-batches scheduler board).

## Mental model

```
target hop risk ↓
    │
    ├─ feature proximity to last batch?  → hop cosine IW
    ├─ Mahalanobis / energy OOD score?   → mahal IW (far → downweight)
    ├─ gradient aligns with D_spe?       → DGA (Fan–Grangier–Ablin)
    └─ typed (ĉ, δ̂) gate bank⊕hop?      → attr adapter
```

All of the above write **weights**; the learner stays weighted Ridge (or the existing adapter).

## Run (default)

```bash
python3 scripts/run_ood_weight_online_eval.py
python3 scripts/run_ood_weight_online_eval.py --quick
python3 scripts/run_ood_weight_online_eval.py --boards amazon,diabetes,interstate,stock
```

Artifacts land in `results/ood_weight_online_eval/`:

| file | content |
| --- | --- |
| `ood_weight_online_eval.json` | tables + per-seed rows + paths |
| `ood_weight_online_eval.png` | MSE + regret bars per board |
| `path_<board>.png` | online MSE trajectories |
| `README.md` | markdown scoreboard |

## Code map

| path | role |
| --- | --- |
| `Python/src/ood_weight_online_eval.py` | boards, weight runners, metrics, plots |
| `Python/src/attribution_adapter.py` | `run_hop_ridge`, `run_dga_ridge`, `run_attr_adapter`, `run_ridge_past` |
| `scripts/run_ood_weight_online_eval.py` | CLI |
| `datasets/datasets.zip` | Diabetes source/target CSVs |

### Stream contract

Anything with `.X`, `.y`, `.batch` (int hop ids) works — same as `AmazonStream`.

### Metrics

| metric | definition |
| --- | --- |
| `online_mse` | mean per-hop MSE predicting batch \(t\) from past |
| `cum_mse` | sum of per-hop MSE |
| `regret` | `cum_mse(method) − cum_mse(uniform)` on same seed |
| `bwt` | MSE on batch 0 after fitting on all past with that weight rule |
| `oracle_gap` | `cum_mse − best method on that seed` |

## Adding a board

1. Write `load_<name>_stream(...) -> AmazonStream` in `ood_weight_online_eval.py`.
2. Register in `load_board` + `BOARDS`.
3. Prefer **real** public tables when available; label synthetic clearly with `meta["surrogate"]=True`.
4. Re-run CLI; update `results/ood_weight_online_eval/README.md` numbers from the JSON (do not hand-wave).

## Adding a weight channel

1. Implement `run_<name>_ridge(stream) -> summary` with keys `online_mse`, `cum_mse`, `online_path`.
2. Register in `run_weight_method` + `WEIGHT_METHODS`.
3. Keep the polarity explicit in the docstring: **OOD-far → downweight** (adaptation) vs **contaminant reject**.

## Interpretation cheat-sheet

| pattern | read |
| --- | --- |
| hop regret ≪ 0 on Amazon | heatmap proximity IW is the right channel |
| dga ≤ uniform, hop fails | density-ratio surrogate weak; use gradient alignment |
| mahal ≈ hop | geometry of last-batch mean drives the gain |
| all regrets ≥ 0 on a board | OOD weighting does not help — report that |

## Related docs

- `docs/method/DGA_hop_note.tex` — DGA mapping
- `docs/method/Attribution_adapter_note.tex` — hop / attr adapter
- `results/attribution_adapter/README.md` — Amazon hop vs DGA live table
