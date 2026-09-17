# PO-risk board

PO+MSE both break → freeze that layer's training.
RF PO-risk should not collapse first; serving MSE is likelier to break.
MSE broken, PO quiet → not concept drift; read MMD²(X_new, X_ref).
MMD is vs the reference batch, not pairwise history, not layer reps.
PO broken, MSE holds → watch. No online-bootstrap.

- [electricity freeze board](electricity/board.html)
- [electricity PO×MSE trend](electricity/po_mse_trend.png)
- [electricity MA gate](electricity/batch_size_gate.png)
- [covertype freeze board](covertype/board.html)
- [covertype PO×MSE trend](covertype/po_mse_trend.png)
- [covertype MA gate](covertype/batch_size_gate.png)
- [airlines MA gate](airlines/batch_size_gate.png)
