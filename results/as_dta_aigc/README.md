# AS-DTA — AIGC detection diagnosis board

As observed, benign evasion inflates the PO-risk uniformly, while malignant adversarial injection induces extreme heteroscedasticity (peaks) in the PO-risk distribution.

## REF detector

```
{
  "accuracy": 0.9596153846153846,
  "f1": 0.9596928982725528,
  "auc": 0.9943195266272189,
  "n": 520,
  "backend": "sklearn_tfidf",
  "train_accuracy": 0.9995192307692308,
  "train_auc": 1.0,
  "n_train": 2080,
  "n_live": 1400,
  "n_pairs": 2000
}
```

## LIVE diagnoses

| batch | mean PO | peak gap | CVaR | kurtosis | acc | action |
| --- | --- | --- | --- | --- | --- | --- |
| ref | 0.1457 | 0.0186 | 0.1592 | 1.25 | 0.966 | `hold` |
| benign_backtranslate | 0.1483 | 0.0186 | 0.1615 | 1.11 | 0.965 | `hold` |
| benign_synonym | 0.1549 | 0.0185 | 0.1682 | 1.09 | 0.963 | `conformal_calibrate` |
| malignant_attack | 0.2317 | 0.0290 | 0.2587 | 0.13 | 0.874 | `freeze_lexical_adapt_syntactic` |

## After adaptation

| batch | mean PO | peak gap | acc |
| --- | --- | --- | --- |
| malignant_attack | 0.2968 | 0.0168 | 0.920 |
| benign_synonym | 0.1549 | 0.0162 | 0.965 |

## Run

```bash
python3 scripts/run_as_dta_aigc.py --quick
python3 scripts/run_as_dta_aigc.py
```

Open `as_dta_diagnosis_board.html` for the interactive-ish board,
or `as_dta_diagnosis_dashboard.png` for the scatter + peaks figure.
