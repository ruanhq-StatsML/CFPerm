# Mag×Var eval — 2×2 + metrics + landing

## Scope

Resolve **benign vs malignant shift** when both raise mean OOD/PO.  
Eval object is a **2×2 contingency** (nothing else as the primary scorecard):

|  | pred benign | pred malignant |
| --- | --- | --- |
| **true benign** | TN | FP |
| **true malignant** | FN | TP |

Quiet/REF is calibration only — not a third class in this table.

## Metrics list

| metric | meaning |
| --- | --- |
| `contingency_2x2` | counts above |
| `accuracy` | (TN+TP)/N |
| `balanced_accuracy` | (mal_recall + ben_recall)/2 |
| `malignant_recall` | TP/(TP+FN) |
| `malignant_precision` | TP/(TP+FP) |
| `benign_recall` | TN/(TN+FP) |
| `J_gap` | $J_{\mathrm{mal}}-J_{\mathrm{ben}}$ merge separation |
| `d_mu_benign` / `d_mu_malignant` | Magnitude components |
| `d_s_benign` / `d_s_malignant` | Variance components |
| `magnitude_only_acc` | weak baseline (Δμ midpoint split) |

Gate = Magnitude × Variance (+ ΔError when Y exists). Peak is not the headline.

## Boards (dataloaders)

`Python/src/as_dta_shift_data.py` — one API:

```python
from as_dta_shift_data import load_board, load_boards, BOARD_NAMES
board = load_board("hc3_finance")
```

| board | source |
| --- | --- |
| `hc3_all` / `hc3_finance` / `hc3_medicine` | HC3 AIGC |
| `news20_sci` / `news20_talk` | 20 Newsgroups topic pairs |

## Landing (why OOD Mag×Var still has room)

1. **AIGC moderation** — conformal vs freeze without “OOD↑ ⇒ retrain”
2. **Trust & safety text** — synonym drift vs adversarial paraphrase
3. **Industrial verticals** — finance/medicine transfer gates
4. **Uplift / RL monitors** — Mag×Var stratum before actuator spend
5. **API / channel gates** — per-route resolution, CPU bootstrap

## Run

```bash
python3 scripts/run_as_dta_mag_var_boards.py
python3 scripts/run_as_dta_mag_var_boards.py --quick
```
