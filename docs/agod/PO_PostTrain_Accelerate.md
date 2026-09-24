# Accelerating Post-Training with PO-Risk

> **LaTeX:** [`PO_PostTrain_Accelerate.tex`](PO_PostTrain_Accelerate.tex)

Two algorithms (PO only):

- **A — Modality:** freeze low chronic \(L\), dump steps on spike \(S\)
- **B — Sample:** if batch rejected, \(w_i\propto\sqrt{\mathrm{PO}_i}\) on OOD-large rows; else \(w=1\)

Goal: lower \(\mathrm{FLOPs}\) / \(T(\mathrm{Acc}^\star)\); Acc is a constraint.

**Optional by design:** freeze / gated √PO only bite when residual mass is asymmetric. No bite → FLOPs=1 (noop). Bite + Acc held → acceleration. Noop packs are not counterexamples.

**Smoke (multi-pack, Acc constraint \(\Delta\mathrm{Acc}\ge-0.005\)):**
- **Held FLOPs↓ (gate bites):** Affec `po_gated` 0.70 / `po_proto` 0.83; synthetic `po_fuse` 0.78 & \(T^\star=5\); Amazon modality B3 0.50; COCO outdoor/indoor B5g 0.81.
- **Noop (\(\mathrm{FLOPs}=1\)):** Food-101 / Fashion-IQ / COCO in \(M{=}2\) PO-risk compare.
- **FLOPs↓ but Acc fail (negative controls):** MSR-VTT packed B5g; Fashion/MM/Indiana/Micro B5g; Amazon online-select B3.
- **Alg B:** gated √PO wins 6/6 packs vs always-√ (next-MSE / \(T^\star\), not FLOPs).
