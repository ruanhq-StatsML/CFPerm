# PO-Risk for Post-Training (consolidated)

> **LaTeX:** [`AGOD_PO_PostTrain_Consolidated.tex`](AGOD_PO_PostTrain_Consolidated.tex)  
> MVP card: [`PO_Boost_MVP.md`](PO_Boost_MVP.md) · entry: `python3 scripts/po_boost_mvp.py`

## One-line

PO → (modality freeze/dump) ⊥ (reject ⇒ \(w_i\propto\sqrt{\mathrm{PO}_i}\) on OOD-large rows) → faster post-train under Acc constraint.

## Data snapshot

**Affec \(M{=}5\):** `po_gated` FLOPs **0.70**, ΔAcc −0.39% (ship); `po_proto` FLOPs 0.83, ΔAcc +0.20%.

**Food/Fashion/COCO \(M{=}2\):** FLOPs = 1 — honest noop.

**Synthetic MVP:** `po_fuse` FLOPs **0.778**, \(T(\mathrm{Acc}^\star)=5\), ship=Y; reject = 1×hop_oos + 11×proxy.

**Sample weights:** gated √PO beats always-√ on 6/6 prior packs.
