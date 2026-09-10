# COCO CLIP outdoor vs indoor subset

Built from `cat-state/clip-embeddings` open-clip ViT-B/32 COCO embeds.

- `img_feats.npy` / `txt_feats.npy`: (6000, 512)
- Batch W: outdoor caption keywords vs indoor caption keywords (disjoint)
- `labels.npy`: factorized first-token of caption (for PO-risk path)

```bash
python3 scripts/run_multimodal_fsds_align.py
```
