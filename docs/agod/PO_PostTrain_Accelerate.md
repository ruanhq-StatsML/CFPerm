# Accelerating Post-Training with PO-Risk

> **LaTeX:** [`PO_PostTrain_Accelerate.tex`](PO_PostTrain_Accelerate.tex)

Two algorithms (PO only):

- **A — Modality:** freeze low chronic \(L\), dump steps on spike \(S\)
- **B — Sample:** if batch rejected, \(w_i\propto\sqrt{\mathrm{PO}_i}\) on OOD-large rows; else \(w=1\)

Goal: lower \(\mathrm{FLOPs}\) / \(T(\mathrm{Acc}^\star)\); Acc is a constraint.

**Smoke:** Affec `po_gated` FLOPs 0.70; synthetic `po_fuse` 0.78 & \(T^\star=5\); gated √PO wins 6/6 packs.
