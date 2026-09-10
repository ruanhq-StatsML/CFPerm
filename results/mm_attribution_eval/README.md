# Multimodal attribution evaluation

## How to evaluate the logic

1. **Define GT modality from the batch mechanism**, not from the scores.
   - view shift (Frontal/Lateral) → GT = image
   - caption-keyword split → GT = text (constructed control)
   - time-order only → no strong GT (baseline / null)
   - synthetic inject on one block → GT = that block

2. **Score methods** with RF Domain (CS), coord-MMD (CS), PO-risk (CD).

3. **Judge recovery** via
   - `mass_on_gt` (modality VIMP share on GT)
   - `selection_AUC` (coords in GT block as positives)
   - `topk_precision_k20` / full-block precision
   - RF `domain_auc` (was the shift actually detectable?)

4. **Pass criteria (practical)**
   - strong CS controls: `mass_on_gt ≥ 0.6` and `selection_AUC ≥ 0.6`
   - inject recoveries: same, with `domain_auc` clearly > 0.5
   - time-order baseline: near-balanced mass OK if AUC≈0.5–0.6

## Results

| Scenario | GT | Method | domain AUC / PO | mass_on_gt | sel AUC | P@20 |
|----------|----|--------|-----------------|------------|---------|------|
| indiana_view_control | image | rf_domain | 1.0 | 0.9999 | 0.8076 | 1.0 |
| indiana_view_control | image | coord_mmd | — | 0.9848 | 0.9522 | 1.0 |
| indiana_view_control | image | po_risk | 0.1331 | 0.7587 | 0.6101 | 0.8 |
| coco_keyword_control | text | rf_domain | 0.9999 | 0.6765 | 0.6223 | 0.65 |
| coco_keyword_control | text | coord_mmd | — | 0.579 | 0.5683 | 0.65 |
| coco_keyword_control | text | po_risk | 3.3828 | 0.6178 | 0.5517 | 0.7 |
| coco_time_order_baseline | — | rf_domain | 0.5682 | — | — | — |
| coco_time_order_baseline | — | coord_mmd | — | — | — | — |
| coco_time_order_baseline | — | po_risk | 16.9315 | — | — | — |
| coco_time_inject_image | image | rf_domain | 1.0 | 0.9999 | 0.7195 | 1.0 |
| coco_time_inject_image | image | coord_mmd | — | 0.9839 | 0.8998 | 1.0 |
| coco_time_inject_image | image | po_risk | 0.1975 | 0.522 | 0.4991 | 0.45 |
| coco_time_inject_text | text | rf_domain | 1.0 | 0.9997 | 0.7453 | 1.0 |
| coco_time_inject_text | text | coord_mmd | — | 0.9831 | 0.8991 | 1.0 |
| coco_time_inject_text | text | po_risk | 0.31 | 0.6885 | 0.5471 | 0.8 |
