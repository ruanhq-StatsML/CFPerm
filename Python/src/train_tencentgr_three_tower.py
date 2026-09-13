#!/usr/bin/env python3
"""Train the TencentGR three-tower on a dumped cache.

This script is the training entrypoint. It does not download data and it
does not run unless you invoke it. Default device is CPU.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import torch
from torch.utils.data import DataLoader

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from tencentgr.config import DEFAULT_CFG
from tencentgr.dataset import TencentGRDataset
from tencentgr.three_tower import ThreeTowerModel


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default=DEFAULT_CFG["cache_dir"])
    p.add_argument("--out-dir", default="experiments/tencentgr/three_tower")
    p.add_argument("--batch-size", type=int, default=int(DEFAULT_CFG["batch_size"]))
    p.add_argument("--lr", type=float, default=float(DEFAULT_CFG["lr"]))
    p.add_argument("--epochs", type=int, default=int(DEFAULT_CFG["epochs"]))
    p.add_argument("--device", default=str(DEFAULT_CFG["device"]))
    p.add_argument("--tower-dim", type=int, default=int(DEFAULT_CFG["tower_dim"]))
    p.add_argument("--hidden-dim", type=int, default=int(DEFAULT_CFG["hidden_dim"]))
    p.add_argument("--num-workers", type=int, default=0)
    p.add_argument("--max-steps", type=int, default=0, help="0 = full epoch; >0 caps steps (debug).")
    p.add_argument("--dry-run", action="store_true", help="One forward pass, then exit. No training.")
    return p.parse_args()


def move_batch(batch, device: torch.device):
    return {k: v.to(device) if torch.is_tensor(v) else v for k, v in batch.items()}


def main() -> None:
    args = parse_args()
    device = torch.device(args.device)
    train_ds = TencentGRDataset.from_cache(args.cache_dir, split="train")
    test_ds = TencentGRDataset.from_cache(args.cache_dir, split="test")
    model = ThreeTowerModel(
        user_dim=train_ds.user_feat_dim,
        emb_dim=train_ds.emb_dim,
        hidden_dim=args.hidden_dim,
        tower_dim=args.tower_dim,
    ).to(device)
    loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=args.num_workers)
    if args.dry_run:
        batch = move_batch(next(iter(loader)), device)
        out = model(batch)
        skip = {"logits", "click_logit", "conversion_logit"}
        print(json.dumps(
            {k: float(v.detach().cpu()) for k, v in out.items() if k not in skip and torch.is_tensor(v) and v.ndim == 0},
            indent=2,
        ))
        print("dry-run ok; not training.")
        return
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history = []
    for epoch in range(args.epochs):
        model.train()
        running = 0.0
        n_step = 0
        for step, batch in enumerate(loader, start=1):
            batch = move_batch(batch, device)
            opt.zero_grad(set_to_none=True)
            out = model(batch)
            out["loss"].backward()
            opt.step()
            running += float(out["loss"].detach().cpu())
            n_step += 1
            if args.max_steps and step >= args.max_steps:
                break
        row = {
            "epoch": epoch + 1,
            "loss": running / max(n_step, 1),
            "steps": n_step,
            "n_train": len(train_ds),
            "n_test": len(test_ds),
            "heads": "retrieval+click+conversion",
            "last_action": "dropped",
        }
        history.append(row)
        print(json.dumps(row))
        ckpt = {
            "epoch": epoch + 1,
            "model": model.state_dict(),
            "cfg": train_ds.cfg,
            "user_dim": train_ds.user_feat_dim,
            "emb_dim": train_ds.emb_dim,
        }
        torch.save(ckpt, out_dir / "last.pt")
    (out_dir / "history.json").write_text(json.dumps(history, indent=2))


if __name__ == "__main__":
    main()
