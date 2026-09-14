"""Frozen ViT embeddings for image-OOD PO-risk.

Freeze a pretrained ViT, cache pooled embeddings, then fit a
pseudo-outcome learner on ID embeddings (see ``agod.image_ood``).
"""
from __future__ import annotations

import io
import json
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np


def load_vit(model_name: str = "vit_tiny_patch16_224", pretrained: bool = True):
    import timm
    import torch
    from torchvision import transforms

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = timm.create_model(model_name, pretrained=pretrained, num_classes=0)
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    model.to(device)

    cfg = timm.data.resolve_model_data_config(model)
    size = int(cfg.get("input_size", (3, 224, 224))[-1])
    tfm = transforms.Compose(
        [
            transforms.Resize(size + 32),
            transforms.CenterCrop(size),
            transforms.ToTensor(),
            transforms.Normalize(
                mean=cfg.get("mean", (0.485, 0.456, 0.406)),
                std=cfg.get("std", (0.229, 0.224, 0.225)),
            ),
        ]
    )
    return model, tfm, device


def _pil_from_bytes(b: bytes):
    from PIL import Image

    return Image.open(io.BytesIO(b)).convert("RGB")


def embed_pils(model, tfm, device, images: Sequence, *, batch_size: int = 32) -> np.ndarray:
    import torch

    outs: List[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for i in range(0, len(images), batch_size):
            batch = images[i : i + batch_size]
            x = torch.stack([tfm(im) for im in batch]).to(device)
            z = model(x)
            if isinstance(z, (tuple, list)):
                z = z[0]
            outs.append(z.detach().cpu().numpy().astype(np.float32))
            if (i // batch_size) % 20 == 0:
                print(
                    f"  embed {min(i + batch_size, len(images))}/{len(images)}",
                    flush=True,
                )
    return np.vstack(outs) if outs else np.zeros((0, 1), np.float32)


def embed_food101_parquets(
    parquet_dir: Path,
    out_dir: Path,
    *,
    model_name: str = "vit_tiny_patch16_224",
    max_n: Optional[int] = 8000,
    batch_size: int = 32,
    seed: int = 0,
) -> dict:
    """Extract ViT embeddings from Food101 parquet shards."""
    import pandas as pd
    import pyarrow.parquet as pq

    parquet_dir = Path(parquet_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    feat_path = out_dir / "vit_feats.npy"
    if feat_path.is_file() and (out_dir / "labels.npy").is_file():
        X = np.load(feat_path)
        y = np.load(out_dir / "labels.npy")
        return {
            "out_dir": str(out_dir),
            "n": int(len(y)),
            "d": int(X.shape[1]),
            "model": model_name,
            "cached": True,
        }

    paths = sorted(parquet_dir.glob("*.parquet"))
    if not paths:
        raise FileNotFoundError(f"no parquet under {parquet_dir}")

    full = pd.concat([pq.read_table(p).to_pandas() for p in paths], ignore_index=True)
    rng = np.random.default_rng(seed)
    if max_n is not None and len(full) > max_n:
        idx = rng.choice(len(full), size=max_n, replace=False)
        full = full.iloc[np.sort(idx)].reset_index(drop=True)

    print(f"decoding {len(full)} Food101 images …", flush=True)
    images, labels = [], []
    for i, row in full.iterrows():
        im = row["image"]
        if not (isinstance(im, dict) and "bytes" in im):
            raise TypeError(type(im))
        images.append(_pil_from_bytes(im["bytes"]))
        labels.append(int(row["label"]))
        if (i + 1) % 1000 == 0:
            print(f"  decoded {i + 1}/{len(full)}", flush=True)

    model, tfm, device = load_vit(model_name)
    print(f"embedding with {model_name} on {device} …", flush=True)
    X = embed_pils(model, tfm, device, images, batch_size=batch_size)
    y = np.asarray(labels, np.int64)
    np.save(feat_path, X)
    np.save(out_dir / "labels.npy", y)
    pd.DataFrame({"label": y, "idx": np.arange(len(y))}).to_csv(
        out_dir / "df_metadata.csv", index=False
    )
    info = {
        "out_dir": str(out_dir),
        "n": int(len(y)),
        "d": int(X.shape[1]),
        "model": model_name,
        "n_classes": int(len(np.unique(y))),
        "cached": False,
    }
    (out_dir / "embed_meta.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    return info


def make_class_holdout_domains(
    labels: np.ndarray,
    *,
    id_frac_classes: float = 0.7,
    seed: int = 0,
) -> np.ndarray:
    y = np.asarray(labels).ravel().astype(int)
    classes = np.sort(np.unique(y)).copy()
    rng = np.random.default_rng(seed)
    rng.shuffle(classes)
    n_id = max(1, int(round(len(classes) * id_frac_classes)))
    id_set = set(classes[:n_id].tolist())
    return np.asarray(["id" if int(yi) in id_set else "ood" for yi in y], dtype=object)
