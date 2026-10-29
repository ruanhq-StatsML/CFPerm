# AS-DTA on AIGC detection — opportunity & product board

**AS-DTA** = Adversarial-Shift Diagnosis & Targeted Adaptation.

Not “another AI detector” — a **LIVE-batch diagnosis layer** on top of an
AIGC detector that routes *benign paraphrase* vs *malignant adversarial*
shifts to different actuators.

## Pipeline

1. **REF fine-tune** — HC3 Human vs ChatGPT; sklearn TF-IDF (CPU) or HF encoder.
2. **LIVE batches** — `ref` · `benign_backtranslate` · `benign_synonym` · `malignant_attack`.
3. **PO-risk proxy + Peaks** — bootstrap sub-envs → mean, CVaR95, Kurtosis, Peak Gap.
4. **Route** — hold / conformal / freeze-lexical+adapt-syntactic.
5. **Dashboard** — Mean × Peak Gap scatter + histograms + HTML board.

PO-risk proxy (fast): $\mathrm{PO}_i = 1 - \hat p(Y_i\mid X_i)$.

Feature blocks: **lexical** = char_wb n-grams (surface / embedding proxy);
**syntactic** = word 1–2 grams (longer-range / CLS proxy).

## Expected quadrant story

| region | batches | action |
| --- | --- | --- |
| low mean, low peak | REF, back-translate | hold |
| high mean, low peak | synonym (uniform inflate) | conformal calibrate FPR |
| high mean, high peak | TextFooler-lite attack | freeze lexical, adapt syntactic |

Paper line: *benign evasion inflates PO-risk uniformly; malignant injection
induces extreme heteroscedasticity (peaks).*

## Commercial justification (short)

- Buyers: LLM platforms, education/hiring integrity, UGC moderation.
- Differentiator: **shape of drift → different medicine** (not another AUROC).
- Same product grammar as FSDS→RL: typed signal → actuator (here Peaks → freeze / conformal).

## Measured (this checkout, 2000 HC3 pairs)

| batch | mean PO | peak | acc | routed action |
| --- | --- | --- | --- | --- |
| ref | 0.146 | 0.019 | 0.966 | hold |
| back-translate | 0.148 | 0.019 | 0.965 | hold |
| synonym | 0.155 | 0.019 | 0.963 | conformal_calibrate |
| malignant | 0.232 | 0.029 | 0.874 | freeze_lexical_adapt_syntactic |

After freeze+aug on malignant: Peak **0.029→0.017**, acc **0.874→0.920**.
Conformal on synonym: FPR≈5.4%, recall≈0.97.

## Run

```bash
python3 scripts/run_as_dta_aigc.py --quick
python3 scripts/run_as_dta_aigc.py --max-pairs 2000
```

Code: `Python/src/as_dta_aigc.py` · Results: `results/as_dta_aigc/`
· Board: `as_dta_diagnosis_dashboard.png` / `as_dta_diagnosis_board.html`
