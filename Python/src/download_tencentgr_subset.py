#!/usr/bin/env python3
"""Download a TencentGR-10M slice that stays under ~40GB.

Full Hub dump is ~413GB. This pulls:
  seq/part-00000 only (~1.0G)
  item_feat/* (~348M)
  user_feat/* (~88M)
  mm_emb/emb_82_1024_parquet/part-000[0-9][0-9]-* (~18G, shards 00-99)
  mm_emb/emb_84_32_parquet/* (~3.7G)
  indexer.pkl (~503M)

Does not train anything.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--local-dir", default="/workspace/data/tencent_subset")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    Path(args.local_dir).mkdir(parents=True, exist_ok=True)
    cmd = [
        "hf",
        "download",
        "TAAC2025/TencentGR-10M",
        "--repo-type",
        "dataset",
        "--local-dir",
        args.local_dir,
        "--include",
        "seq/part-00000-*.parquet",
        "--include",
        "item_feat/*",
        "--include",
        "user_feat/*",
        "--include",
        "mm_emb/emb_82_1024_parquet/part-000[0-9][0-9]-*",
        "--include",
        "mm_emb/emb_84_32_parquet/*",
        "--include",
        "indexer.pkl",
    ]
    if args.dry_run:
        cmd.append("--dry-run")
    print(" ".join(cmd), flush=True)
    subprocess.check_call(cmd)


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError:
        print("hf CLI not found; install with: pip install huggingface_hub", file=sys.stderr)
        raise
