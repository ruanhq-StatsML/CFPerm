# Recent-control PO re-fit vs OOD batch

## Roles (keep them separate)

When OnlineRFPerm rejects at stream index `t`:

| role | batches | job |
|---|---|---|
| **recent control** `R` | `batch_{t-k}, …, batch_{t-1}` | re-fit μ0 under the *recent* regime |
| **OOD batch** `O` | `batch_t` (the significant one) | score PO, build IPTW weights, fit downstream RF |

`k` = `--n-control` / `n_recent` (default 1 → previous batch only).

Non-reject steps: everyone stays **uniform**; those steps are **not** in the sig-only MSE.

## Pipeline

```text
OnlineRFPerm(f_ref, batch_t) → p_t, FDR
if not reject:
    w = 1                          # uniform
else:
    μ0 ← fit(R)                    # recent control
    PO_i = |Y_i − μ0(X_i)|         # on O only (+ optional μ-gap)
    w_i ∝ √PO_i  (or PO^{1/3})
    fit RF on O with sample_weight=w
evaluate next-batch MSE; average only over reject steps
```

## Why this cut

- Gate already decided **current** is the OOD batch — that is the only place weights should fire.
- μ0 must be refreshed from **recent** data, not a stale probe or a frozen burn-in ref.
- Putting `prev` into the treated pool (`prev∪cur`) mixes recent ID with OOD; prefer scoring PO on `O` alone.

## Window modes

- `recent_ood` (**default**): `R = last k`, `O = current`.
- `prev_cur` (legacy): `T0 = before the pair`, `T1 = prev∪cur`, weights still taken from the current slice.

## Code

- `agod/po_refit.py` — `build_recent_ood_windows`, `RefitWindows`, weight helpers
- `scripts/run_agod_rfperm_po_refit.py` — `--n-control`, `--window-mode`
- `scripts/run_agod_po_cbrt_compare.py` — same flags for √ vs ∛

```bash
PYTHONPATH=. python3 scripts/run_agod_rfperm_po_refit.py \
  --window-mode recent_ood --n-control 1
```
