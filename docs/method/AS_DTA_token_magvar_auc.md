# Token perturb → Mag×Var AUC — benign vs malignant shift

## Scope

**No multimodal.** Separate **benign shift** vs **malignant shift** with
PO-risk **Magnitude × Variance**, plus token-wise **perturb → re-forward**
attribution as the lexical mechanism check. Scorecard = **AUC** (+ 2×2).

## Your logic (followed)

1. Fit detector on REF.  
2. LIVE: benign = uniform synonym; malignant = concentrated subset attack.  
3. Unit $\mathrm{PO}_i=1-\hat p(Y_i\mid x_i)$.  
4. Magnitude $\mu=\mathbb{E}[\mathrm{PO}]$; Variance $s=\mathrm{Std}(r_b)$ on bootstrap sub-envs.  
5. **Token-wise:** drop token $j$ → re-forward → $a_j=|\mathrm{PO}^{\mathrm{drop}\,j}-\mathrm{PO}|$.  
6. Merge $J=\Delta\mu\cdot(1+\kappa\Delta s_{+})$; evaluate with **AUC** across replicates × datasets.

## Why this justifies Mag×Var

- Synonym and attack can both raise $\mu$; **$s$ / $J$** track whether PO mass is
  uniform or concentrated (mixture attack).  
- Token perturb→re-forward asks the same question on the lexical axis:
  malignant triggers concentrate $|a_j|$; benign edits spread more evenly.  
- **AUC** is the clean public metric: can a Mag×Var score rank malignant above benign?

## Datasets

`hc3_all`, `hc3_finance`, `hc3_medicine`, `news20_sci`, `news20_talk`

## Results (quick suite)

Pooled AUC (malignant=1):

| score | AUC |
| --- | --- |
| $J$ (Mag×Var merge) | **0.984** |
| $\Delta s$ (variance) | **0.926** |
| $\Delta\mu$ (magnitude) | 0.993 |
| token top-share | ~0.52 (mechanistic; board-dependent) |

Per-board AUC($J$) ≈ 0.97–1.00. Magnitude alone can look strong when attacks
are heavy; **variance remains the resolution axis** when means sit in the same
band — $J$ merges both for ops.

## Run

```bash
python3 scripts/run_as_dta_token_magvar_auc.py
python3 scripts/run_as_dta_token_magvar_auc.py --quick
```
