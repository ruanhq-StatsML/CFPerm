# Cold-start logic (v3): representation-first

## One sentence
**Attribute to find subspace \(S\) → separate representations on \(S\) before touching the head → expand the head only when novel mass / footprint is ready → optional TTA restricted to \(S\).**

Weights (soft / hard / fused) still matter — they fire at **step 3** (expand). Step 2 is geometry, not importance weighting.

## Four steps

### 1. Detect + attribute (RFPerm → FSDS / CFPerm)
- Gate: LIVE/REF error rise.
- Attribute: RF-domain VIMP on batch \(W\in\{0,1\}\) → top-\(K\) subspace \(\mathcal{S}\) → membership \(m_i\).
- Role: **focus**, not classify. On CIFAR-like boards this is a color/texture block; on Office-Home-like boards a domain-style block. Novel already separates in \(\mathcal{S}\); the model just has not used it.
- Boundary: low label Jaccard → `incomparable` → skip directed adapt.

### 2. Representation separation (**key / "instant learn"**)
- Sparse novel: **do not** rush the classification head — a few rows overfit and drag old boundaries.
- On \(\mathcal{S}\), compute known-class prototypes from REF.
- For high-membership LIVE rows whose distance to the nearest old prototype is \(<\) margin, apply a margin push:
  \[
  x_i^{\mathcal{S}} \leftarrow x_i^{\mathcal{S}} + \eta\, m_i\,(\mathrm{margin}-d_i)\,\widehat{(x_i^{\mathcal{S}}-\mu_{\mathrm{nn}})}.
  \]
- Persist a mean \(\mathcal{S}\)-offset for inference.
- This step only needs a few samples: it does **not** teach “what the new class is”, only “stop collapsing into old centers”.

### 3. Dynamic head expand ("learn" completion)
- Buffer novel / high-mem rows across batches.
- Expand when `buffer ≥ N*` **or** attributed footprint is stable (VIMP-top mass + membership concentration).
- Then joint-fit new nodes with **novel buffer + known REF replay**, using the usual weight channel (`auto` → fused if sharp, else soft).
- This is the completion of learning — it happens **after** separation.

### 4. Optional TTA (cautious)
- TENT / SHOT-style entropy min is allowed, but **only on parameters tied to \(\mathcal{S}\)** (here: `coef_[:, S]` + light intercept). Freeze the rest.
- Revert if entropy rises.

## Why this order
| step | answers | if you skip it |
| --- | --- | --- |
| 1 attribute | *where* the shift lives | weights / sep aim at the wrong coords |
| 2 sep | *stop collapsing* with few rows | head overfits or pulls old boundary |
| 3 expand | *name* the new class with enough mass | novel_acc stays open-set / rejected |
| 4 TTA | *fast* BN/head nudge on \(S\) | optional; can drift if global |

## Relation to sample-level weights
- Feature-level FSDS = subspace picker (\(W\): LIVE vs REF).
- Sample-level \(w_i\) = SGD share at **expand**.
- Repr-sep = geometric preconditioner so expand does not fight collapsed features.
- Downsample vs \(s(w)\) audit still applies at expand time.

## Code
- `Python/src/repr_sep_cold_start.py` — `adapt_batch_repr_sep`, `ColdStartState`
- `category_shift_boards.run_cgcd(..., mode="repr_sep")`
- Legacy weight-only path: `weighted_cold_start.adapt_batch` (still default ablation)

```bash
python3 scripts/run_weighted_cold_start.py --quick
python3 scripts/run_repr_sep_cold_start.py --quick
```
