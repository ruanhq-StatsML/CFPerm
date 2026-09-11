# Blockwise Weights on a Few Batches (or Streaming Data)

Not online learning. That phrase does not apply. There are three or four batches — or, if they arrive one after another, call it streaming data. Same object either way.

PO-risk / VIMP of $\hat\tau$ gives **blockwise weights**. Use those weights to adapt for distribution shift. That is the claim. Nothing else.

$\varphi$ remains a discrepancy, not a CATE. Unique CS / CD decomposition remains impossible. No prototype here.

---

Partition $x$ into blocks. On the batches you have (three or four, or a stream of the same), fit the usual DR / R-learner residual and $\hat\tau(x)$. Blockwise importance of $\hat\tau$ is

\[
w = (w_1,\ldots,w_B),\qquad w_b \ge 0,\ \sum_b w_b = 1.
\]

Adapt by reweighting / correcting with $w$. If another batch shows up, recompute $w$. That is still the same batch procedure, not online learning.

Stop there.
