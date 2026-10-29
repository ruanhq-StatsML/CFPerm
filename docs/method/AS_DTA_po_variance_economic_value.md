# Economic value of PO-risk variance routing (justify)

Cost model (industrial AIGC moderation):  
$C = 50\cdot\mathrm{FP} + 20\cdot\mathrm{FN}$ USD per item batch, reported as **cost / 1k**.

PO-variance routing (Magnitude → Std/Peak/CVaR/Error) does not invent a new detector —
it **picks the right cheap drug** and **avoids the wrong expensive one**.

## 1. Measured dollar lift (HC3 finance / medicine LIVE)

| LIVE cell | shape via PO var | actuator | cost/1k before → after | $ saved / 1k |
| --- | --- | --- | --- | --- |
| finance · synonym | benign | conformal | 7629 → **4233** | **≈ $3.4k (45%)** |
| finance · malignant | malignant | freeze+aug | 8404 → **4067** | **≈ $4.3k (52%)** |
| medicine · malignant | malignant | freeze+aug | 3992 → **2358** | **≈ $1.6k (41%)** |
| medicine · synonym | benign | conformal | ~flat | (already calibrated) |

**Suite total (oracle = routed, shape detection 100%):** **$22,430** avoided on eval LIVE.  
Malignant recall lift mean **+4.7pp** *and* precision up — not a recall tax.

## 2. Where the money comes from (three ledgers)

1. **Benign ledger (FP-heavy)** — uniform OOD raises false AI flags.  
   Variance quiet → conformal only (move $\tau$). Finance synonym: FPR-driven cost cut **~45%**.  
   Full retrain / freeze here = wasted GPU + risk of overfitting attack patterns that aren't there.

2. **Malignant ledger (FN-heavy)** — concentrated evasion punches FN.  
   Variance / Peak / ΔError spike → freeze-lexical + adv-aug. Finance malignant: cost/1k **−52%**, recall **0.887→0.960**.  
   Conformal alone cannot fix feature-collision; threshold moves don't restore lexical invariance.

3. **Routing ledger (wrong-drug tax)** — treating all OOD the same destroys value:
   - Malignant → conformal: FN stays → keep paying ~$20×FN (integrity).
   - Benign → freeze/retrain: pay adapt cost + possible FP churn for no Peak disease.
   - “OOD high ⇒ retrain always”: mixes ledger 1 and 2 → overspend on benign, under-fix malignant.

PO variance is the **classifier of which ledger you are in**. Detection accuracy on this suite: **8/8**.

## 3. Contrast vs anti-patterns (economic)

| Policy | Economic failure |
| --- | --- |
| Hard-filter top-PO | Deletes Peak/CVaR signal → silent malignant; selection-biased uplift; FN cost hidden |
| Pseudo-label high-PO | Fake metric lift; retrain on noise; dollars look better on paper, production FN/FP don't |
| Always retrain on OOD↑ | Benign cells burn adapt budget; malignant still needs the *right* head — variance skipped |
| Hold forever | Finance synonym alone leaves ~$3.4k/1k on the table |

## 4. Unit economics (pitch)

- At **$50 FP / $20 FN**, one finance malignant batch (~1.2k) drops expected loss **~$10.1k → ~$4.9k**.
- Shape gate is **CPU bootstrap of PO** (seconds), not a GPU retrain — routing cost ≪ medicine cost.
- Annualize: $\text{saved/1k} \times \text{LIVE volume/1k}$. Even sparse malignant weeks dominate ROI if unrouted FN accrues.

## 5. One-line justify

> PO-risk variance is an **allocation rule for moderation dollars**: quiet variance → spend on threshold (cut FP); spiked variance → spend on lexical freeze (cut FN). On industrial HC3 LIVE that allocation recovered **~$22k** on the eval set with **100%** shape hit-rate — the economic value is right-drug selection, not a bigger model.
