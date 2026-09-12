# Online modality selection & LR adjustment — how it works, how well

## Decision rule (implemented)

```
each online window Dt vs reference D0:
  for m in {text, image}:
      AUC_m, VIMP_m, PO_m  ← RF domain / residual probes on X_m
      g_m ← Normalize(AUC_m·VIMP_m + γ·PO_m)
  α ← EMA( Softmax(g / τ) )          # soft selection
  LR_m ← lr0 · (β + (1-β)·α_m·|M|)   # continuous adjustment
  if α_m < θ: drop ∂L/∂θ_m           # hard selection (save BWD)
```

| Piece | Role |
|---|---|
| MSG `g_m` | ranks which modality is drifting |
| Softmax `α_m` | turns ranks into a budget over modalities |
| `LR_m` | spends more gradient steps on high-α modalities |
| Gate `α_m<θ` | refuses to update near-zero modalities → FLOPs↓ |

B1 disables selection (α=½). B2 selects on AUC only. B3 selects on full MSG.

## What “effect” means here

1. **Selection fidelity (upstream):** does α track the drifting channel?
   On Amazon category shifts, α_image typically rises (image pack changes more
   across categories than English review text templates).
2. **Compute effect:** hard gate skips one proj BWD → rel adapt FLOPs ~0.6–0.8
   (encoder FWD still paid). Wall-clock shrinks less than FLOPs on CPU.
3. **Task effect (downstream):** holdout ΔAcc after the update.
   Smoke often shows B3 ≈ B1 on ΔAcc (sometimes slightly worse) while FLOPs↓
   ⇒ **efficiency win, not yet an accuracy win**.
4. **Cost-utility:** holdΔAcc / relFLOPs. If ΔAcc holds and FLOPs fall, utility↑.

## Reading the smoke (typical)

- B3 **selects image** more often (`frac_image_selected` high, text gated).
- B3 **FLOPs_ratio vs B1 < 1**.
- B3 **holdΔacc_gap vs B1 ≈ 0 or slightly negative** on short smoke.
- Therefore: online selection/adjustment is doing the *intended control*;
  claiming SOTA accuracy from this alone would be overclaim. Fair claim:
  **iso-loss / iso-ΔAcc adaptation with lower update cost**.

## Failure modes to watch

- Over-gate: θ too high / τ too sharp → skip the modality that still carries label signal.
- AUC-only (B2): can chase covariate shift that is irrelevant to Y.
- Encoder-dominated cost: gating proj BWD helps, but big savings need also
  skipping unused modality encoder FWD (stricter hard path).
