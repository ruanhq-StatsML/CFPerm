# AS-DTA — AIGC detection diagnosis board

As observed, benign evasion inflates the PO-risk uniformly, while malignant
adversarial injection induces extreme heteroscedasticity (peaks) in the
PO-risk distribution.

## REF detector (HC3, sklearn TF-IDF)

| split | acc | auc | n |
| --- | --- | --- | --- |
| train | 1.000 | 1.000 | 2080 |
| val | 0.960 | 0.994 | 520 |
| LIVE ref | 0.966 | 0.994 | 1400 |

Backend: word 1–2 grams (syntactic proxy) + char_wb 3–5 (lexical / surface proxy).

## LIVE diagnoses (K=28 bootstrap)

| batch | mean PO | peak gap | CVaR | kurtosis | acc | action |
| --- | --- | --- | --- | --- | --- | --- |
| ref | 0.1457 | 0.0186 | 0.1592 | 1.25 | 0.966 | `hold` |
| benign_backtranslate | 0.1483 | 0.0186 | 0.1615 | 1.11 | 0.965 | `hold` |
| benign_synonym | 0.1549 | 0.0185 | 0.1682 | 1.09 | 0.963 | `conformal_calibrate` |
| malignant_attack | 0.2317 | 0.0290 | 0.2587 | 0.13 | 0.874 | `freeze_lexical_adapt_syntactic` |

## After targeted adaptation

| batch | mean PO | peak gap | acc | note |
| --- | --- | --- | --- | --- |
| malignant_attack | 0.2317 → 0.2968 | **0.0290 → 0.0168** | **0.874 → 0.920** | freeze lexical + adv-aug syntactic head; peaks & accuracy improve |
| benign_synonym | — | — | 0.963 → 0.965 | conformal thr≈0.44, FPR≈5.4%, recall≈0.97 |

## Quadrant read

- **Bottom-left** REF / back-translate → hold
- **Bottom-right** synonym (mean↑, peak flat) → conformal
- **Top-right** malignant (mean↑, peak↑, CVaR↑) → freeze lexical / adapt syntactic

## Run

```bash
python3 scripts/run_as_dta_aigc.py --quick
python3 scripts/run_as_dta_aigc.py --max-pairs 2000
```

Code: `Python/src/as_dta_aigc.py`  
Artifacts: `results/as_dta_aigc/as_dta_diagnosis_dashboard.png`,
`as_dta_diagnosis_board.html`
