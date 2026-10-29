# Watermark forms × uplift monitoring — negative effects & justify

Sparse watermark evidence is useful as a **gate**, not as a **label**. When the
same wash / missing / gray pattern is used to monitor an uplift model
(`W=batch`, `Y=R` / PO-risk), each watermark *kind* can inject selection bias,
label noise, or policy harm. The critical split is **intent of the wash**:

| Intent of wash / absence | Observational twin | Correct posture | Wrong posture (negative effect) |
| --- | --- | --- | --- |
| **Evade / break rules** (destroy WM to fool detector, spoof human) | soft/hard wash below diagonal; strip api_log | treat as *attack-shaped* → PO+FN / freeze-lexical / audit | treat as human → FN, poisoned uplift Y |
| **Privacy / autonomy** (open model, redact logs, strip WM, local rewrite) | same wash / missing / no api_log | stratify; do **not** punish; never force Y=AI | treat as attack → FP, trust/legal harm, biased τ̂ |
| **Benign OOD** (synonym / domain paraphrase) | soft wash + uniform OOD Δμ | conformal / calibrate | full retrain or attack actuator |

## 1. Why this matters for uplift monitors

Uplift / meta-learner monitors estimate treatment effects on reward
(`R_PO = E[τ̂²]`, ATE-proxy). Watermark fields enter as:

1. **Covariate / stratum** — fine, if missingness is modeled.
2. **Hard filter** (drop `strong_pos`, only score `missing`) — **selection bias**.
3. **Pseudo-label for Y** (washed ⇒ human) — **label noise**, breaks PO+FN.

Negative effects are not “watermark is bad”; they are **mis-using a sparse
provenance channel as a causal/label channel**.

## 2. Per-kind negative effects (elaborate)

### `api_log` (provenance)
- **Upside:** strongest short-circuit; almost no wash.
- **Negative:** covers closed API only. Privacy users on open/local models are
  systematically `missing` → monitor over-weights vendor traffic; τ̂ on
  “AI-heavy” batches is conditioned on *logged* AI, not all AI.
- **Justify:** use as positive evidence when present; **never** interpret
  absence as human or as attack.

### `hard_wm` / `soft_wm` (crypto / statistical)
- **Upside:** p-value + robust probe separates clear positives.
- **Negative (evasion):** adversary paraphrases / rewrites to wash bits → looks
  like gray/neg. If policy says “no WM ⇒ allow”, rule is broken.
- **Negative (privacy):** same rewrite to remove tracking WM → identical
  geometry. Punishing wash = punishing privacy-preserving transform.
- **Justify:** wash is **ambiguous intent**. Route gray → human / secondary
  PO+FN with strong task Y; do not auto-escalate wash → sanction.

### `fingerprint` (passive)
- **Upside:** works somewhat on open models (where api_log fails).
- **Negative:** high FP under OOD / privacy style shift (non-native writers,
  domain jargon). Feeding FP into uplift Y invents fake AI mass → false
  uplift alarms and bad RL rewards.
- **Justify:** never `strong_pos` alone; at most gray prior.

### `user_declare` (policy)
- **Upside:** rare honest strong_pos.
- **Negative:** strategic non-declare (evasion) vs privacy non-declare look
  identical; incentive to lie if declare triggers sanctions.
- **Justify:** treat as sparse policy bit; missing declare ≠ negative test.

### `missing`
- **Negative if misread:** “detector off / too short / wrong family” ≠
  “watermark failed attack” ≠ “privacy redaction”.
- **Justify:** structural absence → hand to PO+FN with **task labels**, not
  watermark-implied Y.

## 3. Same observational pattern, opposite policy

```
wash / no-api / gray
        │
        ├─ prior: attack surface high, spoof history  → intentional_injection path
        ├─ prior: privacy mode, open-weight, redact   → hold / stratify / no sanction
        └─ prior: domain OOD only                     → conformal_calibrate
```

**Same evidence, three intents.** The watermark kind alone cannot justify the
actuator; you need (a) task Y / PO+FN shape, (b) product context (privacy mode
flag, channel), (c) never collapse wash → single label.

## 4. Justify: what watermark *should* do under uplift monitoring

1. **Gate, don’t label.** Strength routes *how you look* (short-circuit /
   vague / PO+FN); task Y decides *whether adaptation worked*.
2. **Stratify τ̂ by WM stratum** (`strong_pos` / `gray` / `missing`) to detect
   coverage bias — don’t pool and call it global uplift.
3. **Separate cost accounts:**
   - Evade-missed (FN on malicious wash) → security cost.
   - Privacy-punished (FP sanction on redaction) → trust/legal cost.
   Mixing them is the same fallacy as “OOD high ⇒ retrain”.
4. **Sparse evidence is a feature:** low coverage forces honesty that most
   LIVE rows have no watermark Y — uplift must stand on outcome R and strong
   task labels, not on WM.

## 5. One-line product claim

> Watermark forms reduce **provenance uncertainty** when present; they increase
> **policy risk** when wash/missing is treated as intent. For uplift monitors,
> WM is a stratum and a short-circuit gate — never a substitute for Y, and never
> a single-bit “attack vs privacy” classifier.
