# COCO CLIP early vs late image_id (order proxy)

Not a keyword split. Sorted by `image_path` / numeric id from
`cat-state/clip-embeddings` COCO open-clip ViT-B/32:

- early window: smallest 4000 ids → sample 3000
- late window: largest 4000 ids → sample 3000

This is an **order / pseudo-temporal** batch (COCO has no real timestamps).

Compare with `coco_outdoor_indoor/` (keyword-defined CS, text-leaning by construction).
