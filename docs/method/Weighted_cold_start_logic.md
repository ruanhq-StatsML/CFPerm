# Cold-start logic (locked): fused sample_w × subset_w

## One sentence
Pure **membership** diagnoses *where* the shift lives; cold-start **adapts** with **fused weights** = soft `sample_w` × hard observation `subset_w`, after a boundary gate.

## Pipeline
1. **Detect** — LIVE/REF error rise (OnlineRFPerm-style gate).
2. **Attribute** — FSDS RF-domain VIMP → top-K subspace → membership (diagnostic).
3. **Boundary** — Jaccard(label supports) low → `incomparable` → skip directed adapt.
4. **Fuse** — membership ∪ private-live labels (discovered / pseudo ids).
5. **Weights (default cold-start channel)**
   - `sample_w = 1 + λ_eff · membership_fused` (soft, mean-1)
   - `subset_w` = boost_eff on hard observation subset `{i : membership ≥ q}`
   - **`fused_w = sample_w × subset_w`** ← what `partial_fit` uses
   - **temper**: if membership is diffuse, shrink fused_w toward uniform and clip ratio
6. **Adapt** — `partial_fit(..., sample_weight=fused_w)` + mild amplify of attributed features.

## Why this split
| channel | role |
| --- | --- |
| membership only | diagnosis / ablation — where the footprint is |
| sample_w | soft emphasis on high-membership rows |
| subset_w | hard observation subset — concentrate budget on the localized set |
| fused | cold-start default — soft × hard |
| incomparable | honesty — do not directed-adapt when label supports barely overlap |

Digits/amazon show weights are **not** unconditionally helpful: when uniform is already strong or the subset is wrong, fused can over-emphasize sparse novel and hurt known retention. Boundary + sparsity still gate the update.

## Code
- `weighted_cold_start.adapt_batch(..., weight_mode="fused")` — default
- `fsds_cfperm_weight_benchmark` — amazon / digits / news20 / synth_cgcd ablation

```bash
python3 scripts/run_weighted_cold_start.py --quick
python3 scripts/run_fsds_cfperm_weight_benchmark.py --quick
```
