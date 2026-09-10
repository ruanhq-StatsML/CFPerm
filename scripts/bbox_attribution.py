"""Bounding-box block helpers for hierarchical attribution boards.

Builds a fixed-length bbox feature vector per image from COCO-style
annotations, and exposes subgroup slices for L2b mass reporting:

  geo | category_hist | top_boxes

This is a board module (evidence). Business interpretation stays outside.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any

import numpy as np


def build_bbox_features(
    image_ids: np.ndarray,
    id2hw: dict[int, tuple[int, int]],
    anns_by_image: dict[int, list[dict[str, Any]]],
    category_ids: list[int],
    *,
    top_m: int = 5,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return (n, d_bbox) float32 matrix + schema."""
    cat2i = {c: i for i, c in enumerate(category_ids)}
    d_cat = len(category_ids)
    d_geo = 8
    d_top = top_m * 6
    d_bbox = d_geo + d_cat + d_top
    X = np.zeros((len(image_ids), d_bbox), dtype=np.float32)

    for i, iid in enumerate(map(int, image_ids)):
        W, H = id2hw.get(iid, (1, 1))
        anns = anns_by_image.get(iid, [])
        if not anns:
            continue
        areas, cxs, cys, rows = [], [], [], []
        hist = np.zeros(d_cat, dtype=np.float32)
        tot = 0.0
        for a in anns:
            x, y, w, h = a["bbox"]
            area = (w * h) / max(W * H, 1)
            cx = (x + w / 2) / max(W, 1)
            cy = (y + h / 2) / max(H, 1)
            areas.append(area)
            cxs.append(cx)
            cys.append(cy)
            tot += area
            hist[cat2i[a["category_id"]]] += 1
            rows.append(
                (
                    area,
                    cx,
                    cy,
                    w / max(W, 1),
                    h / max(H, 1),
                    cat2i[a["category_id"]] / max(d_cat - 1, 1),
                )
            )
        areas = np.asarray(areas)
        cxs = np.asarray(cxs)
        cys = np.asarray(cys)
        X[i, 0] = np.log1p(len(anns))
        X[i, 1] = min(tot, 5.0)
        X[i, 2] = float(areas.mean())
        X[i, 3] = float(areas.max())
        X[i, 4] = float(cxs.mean())
        X[i, 5] = float(cys.mean())
        X[i, 6] = float(cxs.std()) if len(cxs) > 1 else 0.0
        X[i, 7] = float(cys.std()) if len(cys) > 1 else 0.0
        if hist.sum() > 0:
            hist = hist / hist.sum()
        X[i, d_geo : d_geo + d_cat] = hist
        rows.sort(reverse=True)
        for j, r in enumerate(rows[:top_m]):
            base = d_geo + d_cat + j * 6
            X[i, base : base + 6] = r

    schema = {
        "d_bbox": d_bbox,
        "geo_slice": [0, d_geo],
        "category_hist_slice": [d_geo, d_geo + d_cat],
        "top_boxes_slice": [d_geo + d_cat, d_bbox],
        "top_M": top_m,
        "n_categories": d_cat,
        "category_ids": list(map(int, category_ids)),
        "geo_names": [
            "n_boxes_log1p",
            "total_area_frac",
            "mean_area_frac",
            "max_area_frac",
            "mean_cx",
            "mean_cy",
            "std_cx",
            "std_cy",
        ],
        "top_box_layout": ["cx", "cy", "w", "h", "area_frac", "cat_norm"] * top_m,
    }
    return X, schema


def group_annotations(instances: dict) -> tuple[dict, dict, list[int]]:
    id2hw = {im["id"]: (im["width"], im["height"]) for im in instances["images"]}
    by = defaultdict(list)
    for a in instances["annotations"]:
        by[a["image_id"]].append(a)
    cats = sorted({a["category_id"] for a in instances["annotations"]})
    return id2hw, by, cats
