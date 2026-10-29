# AS-DTA business metrics — industrial HC3 + precision-recall

Train on **general HC3**, LIVE on **finance / medicine** (industrial integrity
verticals). Drift shape (Mean×Peak) routes the actuator; value is read on
**precision / recall / FPR / expected cost**, not only accuracy.

## Cost model

- FP (human → flagged AI): **$50** — complaint / trust / legal
- FN (AI missed): **$20** — integrity / compliance miss
- Expected cost = $c_{FP}\cdot FP + c_{FN}\cdot FN$ on the LIVE batch

## Drift-shape detection

Accuracy **100%** (quiet / benign / malignant from Mean×Peak route).

```
{
  "quiet": {
    "quiet": 4,
    "benign": 0,
    "malignant": 0
  },
  "benign": {
    "quiet": 0,
    "benign": 2,
    "malignant": 0
  },
  "malignant": {
    "quiet": 0,
    "benign": 0,
    "malignant": 2
  }
}
```

## Precision-Recall before → after

| vertical | LIVE | recall | precision | FPR | cost/1k | actuator |
| --- | --- | --- | --- | --- | --- | --- |
| finance | benign_synonym | 0.981→**0.960** | 0.767→0.862 | 0.297→0.153 | $7629→$4233 | `conformal_domain_calibrated` |
| finance | malignant_attack | 0.887→**0.960** | 0.753→0.867 | 0.291→0.147 | $8404→$4067 | `freeze_lexical_adv_aug` |
| medicine | benign_synonym | 0.987→**0.987** | 0.874→0.874 | 0.142→0.142 | $3675→$3696 | `conformal_domain_calibrated` |
| medicine | malignant_attack | 0.951→**0.973** | 0.872→0.921 | 0.140→0.083 | $3992→$2358 | `freeze_lexical_adv_aug` |

## ROI rollup

- Cost saved (oracle actuators): **$22430**
- Cost saved (routed / detected shape): **$22430**
- Malignant mean recall lift: **+0.047**

Oracle actuators show medicine value (right drug); routed actuators show end-to-end AS-DTA value (detect shape → pick drug). ROI = avoided FP/FN dollar cost on industrial LIVE.

## Run

```bash
python3 scripts/run_as_dta_business_pr.py --quick
python3 scripts/run_as_dta_business_pr.py
```
