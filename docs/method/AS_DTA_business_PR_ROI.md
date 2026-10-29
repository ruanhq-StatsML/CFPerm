# AS-DTA business metrics — industrial HC3 + precision-recall + ROI

Train on **general HC3**, LIVE on **finance / medicine** (industrial integrity
verticals). Drift shape (Mean×Peak) routes the actuator; value is read on
**precision / recall / FPR / expected cost**, not only accuracy.

## Cost model (business)

| error | meaning | default USD |
| --- | --- | --- |
| FP | human flagged as AI (complaint / trust / legal) | **$50** |
| FN | AI missed (integrity / compliance) | **$20** |

Expected cost on a LIVE batch $= c_{\mathrm{FP}}\cdot\mathrm{FP} + c_{\mathrm{FN}}\cdot\mathrm{FN}$.

## Can drift shape be detected?

**Yes — 100%** on this suite (8 LIVE cells: 2 verticals × 4 kinds).

| true \ pred | quiet | benign | malignant |
| --- | --- | --- | --- |
| quiet | 4 | 0 | 0 |
| benign | 0 | 2 | 0 |
| malignant | 0 | 0 | 2 |

Quiet = REF / back-translate · Benign = synonym · Malignant = TextFooler-lite collision.
Routed actuators matched oracle actuators → detection is not the bottleneck; medicine is.

## Precision–Recall before → after (oracle = routed)

| vertical | LIVE | recall | precision | FPR | cost/1k USD | actuator |
| --- | --- | --- | --- | --- | --- | --- |
| finance | synonym | 0.981→0.960 | 0.767→**0.862** | ↓ | $7629→**$4233** | domain conformal |
| finance | malignant | 0.887→**0.960** | 0.753→**0.867** | ↓ | $8404→**$4067** | freeze+aug |
| medicine | synonym | 0.987→0.987 | 0.874→0.874 | ≈ | $3675→$3696 | domain conformal |
| medicine | malignant | 0.951→**0.973** | 0.872→**0.921** | ↓ | $3992→**$2358** | freeze+aug |

### Recall specifically (your ask)

- **Malignant adaptation** lifts recall in both verticals: finance **+7.3pp**, medicine **+2.2pp** (mean **+4.7pp**), *and* precision rises — not a recall–precision tradeoff.
- **Benign conformal** is a *threshold* drug: finance trades a little recall (−2.1pp) for a large precision lift (+9.5pp) and ~**45%** cost/1k cut (fewer expensive FPs). Medicine already near-calibrated → flat.

## ROI justification

On the eval LIVE batches alone:

| lens | number |
| --- | --- |
| Cost saved (oracle actuators) | **$22,430** |
| Cost saved (routed / detected shape) | **$22,430** (same — detection 100%) |
| Malignant mean recall lift | **+0.047** |

**How to pitch**

1. **Detection value** — Mean×Peak separates quiet / benign / malignant without labels on the attack type; wrong drug is expensive (conformal cannot fix feature-collision; freeze+aug is overkill for synonym).
2. **Medicine value** — Right actuator moves the PR curve and dollar cost: malignant freeze+aug raises recall *and* precision; benign conformal cuts FP-heavy cost under domain shift.
3. **Industrial transfer** — Detector trained on general chat; LIVE is finance/medicine. Domain conformal (match vertical REF FPR) beats blindly targeting 5% on general val.
4. **Unit economics** — At $50 FP / $20 FN, one malignant finance batch (~1.2k items) drops from ~$10.1k → ~$4.9k expected loss; annualize by LIVE volume.

## Run

```bash
python3 scripts/run_as_dta_business_pr.py --quick
python3 scripts/run_as_dta_business_pr.py
python3 scripts/run_as_dta_business_pr.py --cost-fp 80 --cost-fn 25
```

Artifacts: `results/as_dta_business_pr/`
