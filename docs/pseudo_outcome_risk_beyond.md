# Blockwise Weights for Distribution-Shift Adaptation

Not online learning, and not a larger claim. PO-risk / VIMP of $\hat\tau$ gives **blockwise weights**. Those weights are what you use to adapt for distribution shift on the next batch.

The rest of CFPerm is unchanged: $\varphi$ is a discrepancy, not a CATE; unique CS / CD decomposition is impossible; subset localization is all you need. No prototype here.

---

## The claim

Partition $x$ into blocks (e.g. video $\oplus$ audio $\oplus$ text, or any coordinate groups). Fit the usual DR / R-learner residual $\varphi$ and $\hat\tau(x)$ between $D_0$ and $D_t$. Blockwise variable importance (or LOCO) of $\hat\tau$ is a weight vector

\[
w = (w_1,\ldots,w_B),\qquad w_b \ge 0,\ \sum_b w_b = 1.
\]

Adaptation for the shift is **reweighting / correcting by $w$**: more weight on blocks that carry the discrepancy, less on the rest. Then predict on the next batch with that block-weighted correction. The generator (or the rest of the model) is not the object being updated.

That is the whole method.

---

## What this is not

- Not online learning of $\theta$.
- Not a new identification result for CS vs CD.
- Not “the LLM infers the shift.”
- Not a claim that $w$ is a unique root cause — only a weight for adaptation.

---

## What we are willing to claim

Blockwise $w$ from PO-risk VIMP is a practical input to distribution-shift adaptation. Use it as weights. Stop there.
