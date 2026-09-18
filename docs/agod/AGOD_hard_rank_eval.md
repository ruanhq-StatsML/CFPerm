# Hard-sample ranking quality for PO-risk

PO-risk is blunt on purpose: `score_i = |Y_i − μ(X_i)|`.  
When OnlineRFPerm rejects, the useful claim is **not** “MSE always drops”, but:

> **Does this score put the hard OOD rows near the top?**

That is a **ranking** question. Downstream √PO-IPTW is only as good as this ranking.

---

## 1. What counts as a “hard sample”

On reject batch `O = batch_t`, with recent control `R`:

```text
μ_oracle ← fit(R ∪ O) with uniform weights   # diagnostic only
truth_i  = |Y_i − μ_oracle(X_i)|             # higher = harder
```

`μ_oracle` **sees current labels** — so `truth` is a proxy for residual difficulty,  
not something you can deploy. It answers: *if we knew a good fit on this batch,  
which rows would still be hard?*

Candidate scores (no access to `μ_oracle`):

| score | μ used |
|---|---|
| `ref_po` | frozen OnlineRFPerm `f_ref` |
| `probe_po` | rolling prev-batch probe |
| `refit_po` | μ0 **re-trained** on recent control `R`, scored on `O` |

---

## 2. How ranking quality is measured

Let `n = |O|`, default `k = ⌈0.2 n⌉` (top 20%).

### Full-order concordance

- **Spearman ρ** — `corr(rank(score), rank(truth))`  
  Whole-batch ordering. Insensitive to absolute scale (good for residual risk).

- **Pearson** — raw `corr(score, truth)`  
  Secondary; more scale-sensitive.

### Top-hard recovery (the main “难样本” lens)

- **Precision@k = Recall@k** — when both tops have size `k`:

  ```text
  |Top_k(score) ∩ Top_k(truth)| / k
  ```

  “Of the rows PO calls hardest, how many are truly hardest?”

- **Lift@k** — `Precision@k / (k/n)`  
  1.0 = random; 2.0 = twice random recovery of the hard set.

- **NDCG@k** — graded relevance = `truth_i`  
  Rewards putting *very* hard rows higher than merely hard ones.

- **AUROC (top-q% positive)** — label `truth` top-q% as class 1;  
  AUROC of `score` ranking that class. Threshold-free top-hard detection.

We also report Precision@10% / @20% / @30% as a small grid.

---

## 3. Protocol (when we average)

1. Stream batches; OnlineRFPerm burn-in then monitor.  
2. **Only on reject** (`gate=1`): compute `truth`, `ref/probe/refit` scores, metrics.  
3. Non-reject ≡ uniform → **excluded** (would dilute the hard-sample story).  
4. Report **mean over reject batches** per dataset.

Primary table for this claim = ranking metrics.  
Sig-only next-MSE is secondary (IPTW may not convert ranking into MSE).

---

## 4. How to read the numbers (from the current run)

| pack | Spearman ref → probe/refit | story |
|---|---|---|
| stocks (AAPL/MSFT/IWM) | ~0.3–0.6 → **~0.67–0.75** | re-fit/probe clearly better at ordering hard rows |
| waymo_proxy | 0.12 → **~0.53** | frozen `f_ref` almost useless; recent μ helps a lot |
| beijing_pm25 | 0.08 → **~0.40** | same pattern, milder |
| metro | ~0.30 → ~0.32 | already OK; little headroom |

**Takeaway:** emphasize **Precision@k / Spearman / AUROC**, not overall MSE.  
Re-training the PO-learner after RFPerm is justified when it lifts hard-row ranking  
even if √PO-IPTW does not always win next-MSE vs uniform.

---

## 5. Code

- `agod/hard_rank_metrics.py` — `hard_rank_metrics`, `aggregate_hard_rank`
- `agod/po_ref_vs_refit.py` — scores + thin wrapper
- `scripts/run_agod_po_ref_vs_refit.py` — end-to-end reject-only eval
