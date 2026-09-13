#!/usr/bin/env python3
"""Fit a DR pseudo-outcome learner on the dumped TencentGR cache.

Y = 1{any action in the window is click/conversion} by default.
W = 1{last_timestamp >= median} (late vs early batch).
X = user_vec ⊕ mean(history_emb) ⊕ target_emb, optionally compacted.

Does not train the three-tower. Not a CATE estimate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from tencentgr.config import DEFAULT_CFG
from tencentgr.dataset import TencentGRDataset, design_matrix_from_dataset
from tencentgr.dr_po_learner import compact_blocks, fit_dr_pseudo_outcome


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default=DEFAULT_CFG["cache_dir"])
    p.add_argument("--out-dir", default="experiments/tencentgr")
    p.add_argument("--max-n", type=int, default=4000)
    p.add_argument("--n-splits", type=int, default=3)
    p.add_argument("--compact", action="store_true", default=True)
    p.add_argument("--full-x", action="store_true", help="Use full 1056-d mm blocks (slow, high-p).")
    p.add_argument("--keep", type=int, default=32)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--ridge-alpha", type=float, default=1.0)
    p.add_argument("--clip-e", type=float, default=0.01)
    p.add_argument("--y-mode", default="any_click", choices=["any_click", "last_click", "action_rate"])
    return p.parse_args()


def _jsonable(d: dict) -> dict:
    skip = {"phi", "tau_hat", "e_hat", "mu0_hat", "mu1_hat", "scaler_mean", "scaler_scale", "tau_coef"}
    out = {}
    for k, v in d.items():
        if k in skip:
            continue
        if isinstance(v, (np.floating, np.integer)):
            out[k] = v.item()
        else:
            out[k] = v
    return out


def main() -> None:
    args = parse_args()
    ds = TencentGRDataset.from_cache(args.cache_dir, split="all")
    X, Y, W, user_ids = design_matrix_from_dataset(ds, max_n=args.max_n, y_mode=args.y_mode)
    use_compact = args.compact and not args.full_x
    if use_compact:
        X = compact_blocks(X, ds.user_feat_dim, ds.emb_dim, keep=args.keep)
    fit = fit_dr_pseudo_outcome(
        X,
        Y,
        W,
        n_splits=args.n_splits,
        clip_e=args.clip_e,
        ridge_alpha=args.ridge_alpha,
        seed=args.seed,
    )
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_dir / "dr_po_arrays.npz",
        phi=fit["phi"],
        tau_hat=fit["tau_hat"],
        e_hat=fit["e_hat"],
        Y=Y,
        W=W,
        user_ids=user_ids[: len(Y)],
        tau_coef=fit["tau_coef"],
    )
    summary = _jsonable(fit)
    summary.update(
        {
            "compact_x": bool(use_compact),
            "x_dim": int(X.shape[1]),
            "user_feat_dim": int(ds.user_feat_dim),
            "emb_dim": int(ds.emb_dim),
            "n_cache_samples": int(len(ds.samples)),
            "y_definition": {
                "any_click": "1{any action in history+target >= 1}",
                "last_click": "1{target_action >= 1}",
                "action_rate": "mean(action_id over history+target)",
            }[args.y_mode],
            "y_mode": args.y_mode,
            "w_definition": "1{last_timestamp >= median}",
            "last_click_note": (
                "On this slice the last seq event is almost always exposure (action=0), "
                "so last_click is near-degenerate. Default Y is any_click."
            ),
        }
    )
    path = out_dir / "dr_po_subset.json"
    path.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
