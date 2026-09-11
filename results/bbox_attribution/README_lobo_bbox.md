# COCO post-hoc LOBO bounding-box attribution

- form: `data/img_txt/coco_time_order/coco_bbx_text_form/samples_bbx_text.jsonl`
- PO-risk (fit) = 0.0018
- f* = ['nn_center_dist_mean', 'largest_area_frac', 'area_gini', 'cx_std', 'cy_mean']

## Named selection (PO VIMP top)

| rank | feature | vimp |
|---:|---|---:|
| 1 | `nn_center_dist_mean` | 0.06811 |
| 2 | `cx_std` | 0.05569 |
| 3 | `cy_mean` | 0.05416 |
| 4 | `cx_mean` | 0.05326 |
| 5 | `aspect_std` | 0.05113 |
| 6 | `aspect_mean` | 0.04726 |
| 7 | `cy_std` | 0.04579 |
| 8 | `largest_area_frac` | 0.04425 |
| 9 | `h_mean` | 0.04362 |
| 10 | `nn_center_dist_min` | 0.04349 |

## Top-box category counter (among LOBO images)

| category | n_top | mean |Δ PO| |
|---|---:|---:|
| person | 4 | 0.9696 |
| car | 3 | 0.9631 |
| traffic light | 3 | 0.9635 |
| bench | 2 | 0.6768 |
| bed | 2 | 0.9933 |
| toilet | 2 | 0.9903 |
| surfboard | 2 | 0.9963 |
| chair | 2 | 0.9708 |
| dining table | 2 | 0.9857 |
| bowl | 2 | 0.9508 |

## Examples

| image | batch | bbx | category | Δ PO rel | text |
|---:|---:|---|---|---:|---|
| 500859 | 0 | bbx_1 | bench | 0.9694 | bench |
| 330531 | 0 | bbx_2 | pizza | 0.9069 | pizza: A pizza with various toppings is pictured uncooked. |
| 230976 | 0 | bbx_2 | car | 0.9334 | car |
| 270918 | 0 | bbx_2 | person | 0.9207 | person |
| 235788 | 0 | bbx_2 | bus | 0.9036 | bus: A vintage VW bus that is red and white. |
| 391735 | 0 | bbx_2 | person | 0.9588 | person |
| 234147 | 0 | bbx_1 | carrot | 1.0000 | carrot |
| 109869 | 0 | bbx_1 | bench | -0.3841 | bench: A cat standing on top of a wooden bench. |
| 558577 | 0 | bbx_1 | bed | 0.9912 | bed: A man lying down on a bed showing his clean sneaker bottoms. |
| 520437 | 0 | bbx_1 | suitcase | 0.9158 | suitcase |
| 191691 | 0 | bbx_2 | toothbrush | 0.4569 | toothbrush |
| 169226 | 0 | bbx_1 | car | 0.9827 | car |
