"""Handcrafted interpretable bbox features for attribution boards.

Business rule: rank feature NAMES (n_people, area_frac_vehicle, ...),
never opaque embedding coordinates (index 10/768).
"""
from __future__ import annotations

# Re-export builder entry via the runner module path for notebooks.
from pathlib import Path

__all__ = ["KEY_CLASSES", "SUPER_GROUPS"]

KEY_CLASSES = [
    "person",
    "car",
    "truck",
    "bus",
    "bicycle",
    "motorcycle",
    "dog",
    "cat",
    "chair",
    "couch",
    "dining table",
    "tv",
    "laptop",
    "cell phone",
    "bottle",
    "cup",
    "bowl",
    "book",
]

SUPER_GROUPS = [
    "person",
    "vehicle",
    "animal",
    "outdoor",
    "accessory",
    "sports",
    "kitchen",
    "food",
    "furniture",
    "electronic",
    "appliance",
    "indoor",
]

README = Path(__file__).with_name("run_bbox_named_attribution.py")
