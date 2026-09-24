# Adjacent board: DiffDB / Tencent ops·content / 真驱动 / 上线

Not “AUC high = good.” This note elaborates the **logic** behind the panels
and what still sits outside the board.

## 1. DiffusionDB logic

| piece | here |
|---|---|
| \(X\) | prompt TF-IDF tokens |
| \(Y\) | continuous `image_nsfw` → early-quantile binary |
| grain | every 1000 / 2000 prompts in time order |

**Expected reading: weak transfer (~0.60 HGB / LogReg).**

Why that is *correct*, not a failure:

- Style tokens (`mucha`, `alphonse`, …) dominate SelectKBest on chunk \(t\), but
  NSFW is only loosely tied to those tokens; the association is shallow.
- `top_cmean` often shows stopwords / resolution tokens (`8k`, `4k`, `and`) —
  **mean-shift ≠ predictive rank**. Empirically `FSDS∩cmean ≈ 0` on this smoke —
  the two rankings disagree completely. That is the DiffDB tell.
- HGB ≈ LogReg (~0.60 / ~0.62): the weakness is the **association**, not an
  HGB quirk.
- So the board says: “do **not** treat prompt-token FSDS as a stable NSFW
  filter across adjacent batches.” That is the product signal.

If you wanted a stronger DiffDB story later: concat CFG/step/sampler, or CLIP
embeddings — still the same transfer question, different \(X\).

## 2. Tencent OPS panel logic

Ops = convert-leak dropped, **volume kept** (`e_n_exp`, `e_log1p_exp`, …).

**Expected reading: AUC ≈ 0.99 HGB, tops = exposure intensity, high Jaccard.**

Logic:

- Convert is rare; exposure counts are strong, stable correlates.
- Next-chunk transfer of “more exp → more convert” is almost tautological for
  an ops dashboard — useful as a **baseline of what intensity alone buys**.
- HGB ≫ LogReg on ops@2000 (≈0.999 vs ≈0.66) can happen when volume interacts
  nonlinearly; content panel (below) is the check that we are not worshipping
  the nonlinear probe.
- It is **not** a content tip and **not** a claim that we should ship HGB on
  volume features (they may be policy/allocation endogenous).
## 3. Tencent content panel logic

Content = ops − volume. Credit / share / covisit ranks remain.

**Expected reading: early-pair AUC drops (e.g. 0.99 → ~0.65), tops flip to
`i_credit_*` / `i_share_*`; mean AUC may recover later.**

Smoke snapshot: ops−content gap @1000 ≈ **0.16** (“volume explains most
transfer”); @2000 ≈ **0.10**. Tops flip
`e_log1p_exp…` → `i_credit_last / i_share_last`. On content, HGB ≈ LogReg
(~0.84) — residual transfer is not an HGB artifact.
Logic of the **gap** `ops_auc − content_auc`:

| gap | reading |
|---|---|
| large (≥0.1) | most transferable signal was intensity |
| small | non-volume structure still carries next-chunk predictivity |

Content panel is the **hypothesis generator** for attribution-ish features
(credit/share). It still does **not** prove those are true drivers — only that
they transfer better than chance once volume is removed.

## 4. Are there others? (yes)

| pack | board logic |
|---|---|
| **Waymo proxy** | Planted gradual drift. High AUC + high Jaccard = sanity that the probe recovers persistent kinematics. If this failed, distrust the board. |
| **Metro @1000** | High AUC + **low** Jaccard → weather/calendar still predicts traffic, but *which* feats rotate (regime). Grain \(N\) matters: @2000 Jaccard rises. |
| **Beijing PM2.5** | Smooth pollution regimes → high transfer + mid Jaccard. Lagged \(Y\) in \(X\) is intentional meteo practice; treat as known persistence, not a discovery. |

These are the contrasts that make DiffDB/ops/content readable: weak / intensity /
content-residual / planted / shifting / smooth.

## 5. Beyond AUC (this iteration)

| signal | role |
|---|---|
| \(\Delta\bar Y\) | level shift across the cut |
| HGB AUC | nonlinear transfer probe |
| **LogReg AUC** | same selected \(X\); if both move together, not an HGB quirk |
| **Brier** | crude calibration of the probe (still not a ship metric) |
| Jaccard(top-5 across chunks) | driver-set stability |
| **FSDS∩cmean Jaccard** | selection method agreement (rank vs mean-shift) |
| **ops−content gap** | how much of Tencent transfer is volume |

## 6. 真驱动 (true driver) — what the board cannot do

Transfer ≠ causation.

To claim a **true driver** you need at least one of:

1. **Ablation / hold-out feature family** under a fixed policy (does \(\hat Y\)
   and downstream KPI move when the family is removed?).
2. **PO / IPTW / risk post-train** style estimands already in the AGOD docs —
   with Acc constraints, not unconstrained AUC chase.
3. **Intervention or quasi-experiment** (traffic, ranking, creative).

Board contribution only: **shortlist** families (ops intensity vs credit/share
vs prompt tokens) and say which ones *travel*. Shortlist ≠ driver certificate.

Rule of thumb on this board:

- DiffDB tokens: weak travelers → poor driver candidates for NSFW.
- Tencent volume: strong travelers → strong *correlate*, weak *driver claim*
  without policy exogeneity.
- Tencent credit/share: mid travelers after volume drop → **candidates** for
  localize / PO follow-up, not finished science.

## 7. HGB 该不该上线 — ship gate

Hard rule encoded in `summary.json → ship_gate`:

```text
promote_HGB_to_production = false
```

Why, even when AUC ≈ 1.0:

| missing for ship | why board does not cover it |
|---|---|
| serve skew / label delay | adjacent chunk ≠ online arrival process |
| Acc / latency / cost budget | probe is max_depth=3 smoke, not the serve graph |
| calibration under decision threshold | Brier here is diagnostic only |
| abstain / gray / rollback | coverage docs are separate |
| Acc constraint vs PO risk | AGOD post-train path, not this runner |

What the board **can** greenlight:

- investigate a feature family
- drop a weak-transfer family from a tip story
- pick business grain \(N=1000\) vs \(2000\) for dashboards

What needs other tools:

- **true driver** → ablation / PO / experiment  
- **ship model** → holdout Acc + calibration + cost + gray rollback  

## 8. Practical reading card

1. Look at DiffDB: weak? good — do not overclaim tokens.  
2. Look at ops−content gap: large? intensity baseline owned.  
3. Look at content tops: credit/share? queue for localize/PO, not HGB ship.  
4. Look at HGB vs LogReg: agree? association is real for the probe; else dig.  
5. Refuse “AUC high → 上线.”

Artifacts: `results/sample_chunk_adjacent_board/summary.json`
(`ops_content_gap`, `ship_gate`), `SAMPLE_CHUNK_BOARD.md`.

Reason codes (auto): `results/sample_chunk_adjacent_board/reason_codes/` —
see [`Board_Reason_Codes.md`](./Board_Reason_Codes.md).
