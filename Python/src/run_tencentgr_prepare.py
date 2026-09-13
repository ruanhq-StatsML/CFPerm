#!/usr/bin/env python3
"""Load the on-disk TencentGR subset and dump seq/item/user/indexer/mm/samples.

Does not train the three-tower.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from tencentgr.config import DEFAULT_CFG
from tencentgr.dataset import TencentGRDataset


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=DEFAULT_CFG["root"])
    p.add_argument("--cache-dir", default=DEFAULT_CFG["cache_dir"])
    p.add_argument("--max-users", type=int, default=int(DEFAULT_CFG["max_users"]))
    p.add_argument("--max-seq-len", type=int, default=int(DEFAULT_CFG["max_seq_len"]))
    p.add_argument("--split", default="all", choices=["train", "test", "all"])
    p.add_argument("--split-ratio", type=float, default=float(DEFAULT_CFG["split_ratio"]))
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = dict(DEFAULT_CFG)
    cfg.update(
        {
            "root": args.root,
            "cache_dir": args.cache_dir,
            "max_users": args.max_users,
            "max_seq_len": args.max_seq_len,
            "split_ratio": args.split_ratio,
        }
    )
    ds = TencentGRDataset.from_raw(cfg, split=args.split)
    cache = ds.dump_cache(args.cache_dir)
    meta = json.loads((cache / "meta.json").read_text())
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
