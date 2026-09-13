#!/usr/bin/env python3
"""Download Food-101 (HF) validation and pack CLIP image/text features.

  hf download ethz/food101 --repo-type dataset \\
      --include 'data/validation-*.parquet' --local-dir data/raw/food101
  PYTHONPATH=. python3 scripts/pack_food101_clip.py

Writes data/img_txt/food101/{img_feats,txt_feats,labels}.npy
Text = CLIP("a photo of {class}, a type of food").
"""
from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from PIL import Image
from transformers import CLIPModel, CLIPProcessor

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "food101" / "data"
OUT = ROOT / "data" / "img_txt" / "food101"
N_MAX = 5000
SEED = 2026
BATCH = 32
CLIP_ID = "openai/clip-vit-base-patch32"


def _feat(out):
    if torch.is_tensor(out):
        return out
    if getattr(out, "pooler_output", None) is not None:
        return out.pooler_output
    raise TypeError(type(out))


def _decode(cell):
    if isinstance(cell, dict) and "bytes" in cell:
        return Image.open(io.BytesIO(cell["bytes"])).convert("RGB")
    if isinstance(cell, Image.Image):
        return cell.convert("RGB")
    raise TypeError(type(cell))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = sorted(RAW.glob("validation-*.parquet"))
    if not parts:
        raise FileNotFoundError(
            f"missing {RAW}/validation-*.parquet — run hf download first"
        )
    df = pd.concat(
        [pd.read_parquet(p, columns=["image", "label"]) for p in parts],
        ignore_index=True,
    )
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(df), size=min(N_MAX, len(df)), replace=False)
    df = df.iloc[idx].reset_index(drop=True)

    try:
        from datasets import load_dataset_builder

        names = load_dataset_builder("ethz/food101").info.features["label"].names
    except Exception:
        names = [f"class_{i}" for i in range(101)]

    proc = CLIPProcessor.from_pretrained(CLIP_ID)
    model = CLIPModel.from_pretrained(CLIP_ID).eval()
    labels = df["label"].to_numpy().astype(np.int64)
    texts = [f"a photo of {names[int(y)]}, a type of food" for y in labels]

    img_feats, txt_feats = [], []
    with torch.no_grad():
        for s in range(0, len(df), BATCH):
            sl = slice(s, s + BATCH)
            imgs = [_decode(x) for x in df["image"].iloc[sl]]
            tinp = {
                k: v
                for k, v in proc(
                    images=imgs, return_tensors="pt", padding=True
                ).items()
            }
            fi = _feat(model.get_image_features(**tinp))
            fi = fi / fi.norm(dim=-1, keepdim=True)
            img_feats.append(fi.cpu().numpy().astype(np.float32))

            txt_in = {
                k: v
                for k, v in proc(
                    text=texts[s : s + BATCH],
                    return_tensors="pt",
                    padding=True,
                    truncation=True,
                ).items()
            }
            ft = _feat(model.get_text_features(**txt_in))
            ft = ft / ft.norm(dim=-1, keepdim=True)
            txt_feats.append(ft.cpu().numpy().astype(np.float32))
            if (s // BATCH) % 25 == 0:
                print(f"  encoded {min(s + BATCH, len(df))}/{len(df)}", flush=True)

    img = np.concatenate(img_feats, 0)
    txt = np.concatenate(txt_feats, 0)
    np.save(OUT / "img_feats.npy", img)
    np.save(OUT / "txt_feats.npy", txt)
    np.save(OUT / "labels.npy", labels)
    pd.DataFrame(
        {"label": labels, "label_name": [names[int(y)] for y in labels]}
    ).to_csv(OUT / "df_metadata.csv", index=False)
    (OUT / "pack_info.json").write_text(
        json.dumps(
            {
                "source": "ethz/food101",
                "split": "validation",
                "n": int(len(labels)),
                "clip": CLIP_ID,
                "text_template": "a photo of {label}, a type of food",
                "n_class_raw": int(len(np.unique(labels))),
            },
            indent=2,
        )
    )
    print(f"wrote {OUT} img={img.shape} txt={txt.shape}")


if __name__ == "__main__":
    main()
