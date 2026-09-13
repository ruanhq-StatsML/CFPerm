"""Default cfg for the TencentGR-10M disk subset (not the full 10M dump)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


DEFAULT_CFG: Dict[str, Any] = {
    "root": "/workspace/data/tencent_subset",
    "cache_dir": "/workspace/data/tencent_subset_cache",
    "seq_shards": ["seq/part-00000-*.parquet"],
    "dataset_name": "TAAC2025/TencentGR-10M",
    "mm_emb_configs": ["mm_emb_82_1024", "mm_emb_84_32"],
    "mm_emb_dims": {"mm_emb_82_1024": 1024, "mm_emb_84_32": 32},
    "max_seq_len": 50,
    "emb_dim": 1056,
    "split_ratio": 0.8,
    "batch_size": 128,
    "lr": 1e-3,
    "epochs": 5,
    "device": "cpu",
    "max_users": 8000,
    "user_feat_dim": 20,
    "tower_dim": 64,
    "hidden_dim": 128,
    "pad_id": 0,
}


def mm_emb_dirname(cfg_name: str) -> str:
    """Map mm_emb_82_1024 -> emb_82_1024_parquet."""
    name = cfg_name.strip()
    if name.startswith("mm_emb_"):
        name = "emb_" + name[len("mm_emb_") :]
    if not name.endswith("_parquet"):
        name = name + "_parquet"
    return name


def resolve_root(cfg: Dict[str, Any]) -> Path:
    return Path(cfg.get("root", DEFAULT_CFG["root"])).expanduser().resolve()
