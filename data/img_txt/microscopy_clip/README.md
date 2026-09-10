# Microscopy CLIP image + text

From Hugging Face `kvriza8/clip_microscopy_image_text_embeddings`.

- `img_feats.npy` / `txt_feats.npy`: (20936, 512)
- Batch W: short vs long caption (median split)
- `labels.npy`: factorized `caption_summary`

```bash
python3 scripts/run_multimodal_fsds_align.py
```
