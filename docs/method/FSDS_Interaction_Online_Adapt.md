# FSDS × Interaction Value Online Adaptation

**Primary eval = multimodal `interaction(i,j)` visualization** in the
shared embedding — not MSE. MSE is a footnote for the clever-covariate
actuator only.

```
Z = shared_embed(text, video, audio)
Φ(i,j) = f(Z) − f(Z_{-i}) − f(Z_{-j}) + f(Z_{-i,-j})
```

PO-risk still uses RF outcome crossfit for the online gate; the board you
read is the Φ heatmaps (`interaction_TI_early_late.png`, diag profiles).

## Primary: Φ(i,j) Interaction Shift

| scenario | TI diag drop | TI ‖Δ‖F | TA diag drop | IA diag drop | top Δ TI (i,j) |
| --- | ---: | ---: | ---: | ---: | --- |
| healthy | 0.0028 | 0.050 | -0.0014 | -0.0025 | (0,0)=-0.0396 |
| interaction_break | 0.0048 | 0.087 | -0.0002 | -0.0035 | (0,0)=-0.0677 |
| alignment_attack | 0.0028 | 0.049 | 0.0009 | -0.0016 | (0,0)=-0.0393 |
| compound | 0.0047 | 0.086 | 0.0007 | -0.0020 | (0,0)=-0.0675 |

Figures: `interaction_TI_matrix.png` (early | late | Δ),
`interaction_diag_profiles.png` (aligned Φ(k,k)).

## Hop interpretations (`fsds_clever_cov`)

After each hop: RF-PO modality π, top Φ tokens, triangle, typed (c,δ), read.

## Run

```
python3 scripts/run_fsds_interaction_online_adapt.py
python3 scripts/run_fsds_interaction_online_adapt.py --quick
python3 scripts/run_fsds_interaction_online_adapt.py --viz-only
```

Code: `Python/src/fsds_interaction_online_adapt.py`

Eval contract: **look at Φ(i,j)**; use RF-PO only to gate / interpret.
