"""Ground-truth evaluator for localized item subsets.

Plug in a list of GT orders / items → hit / precision / recall @k.
CSV columns accepted (any one is enough):
  - item_id
  - oid / order_id / sku_id / goods_id  (aliases → item_id)
Optional: order_id for order-level coverage.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Union

import numpy as np
import pandas as pd

PathLike = Union[str, Path]

_ITEM_ALIASES = ("item_id", "oid", "sku_id", "goods_id", "product_id", "iid")
_ORDER_ALIASES = ("order_id", "ord_id", "txn_id", "purchase_id")


def _pick_col(df: pd.DataFrame, aliases: Sequence[str]) -> Optional[str]:
    lower = {c.lower(): c for c in df.columns}
    for a in aliases:
        if a.lower() in lower:
            return lower[a.lower()]
    return None


def load_gt_items(path: PathLike) -> Set[int]:
    """Load ground-truth item ids from CSV / parquet / txt (one id per line)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    if p.suffix.lower() in {".parquet", ".pq"}:
        df = pd.read_parquet(p)
    elif p.suffix.lower() in {".txt", ".list"}:
        vals = [ln.strip() for ln in p.read_text().splitlines() if ln.strip()]
        return {int(float(v)) for v in vals}
    else:
        df = pd.read_csv(p)
    col = _pick_col(df, _ITEM_ALIASES)
    if col is None:
        # single-column file
        if df.shape[1] == 1:
            col = df.columns[0]
        else:
            raise ValueError(
                f"GT file {p} needs an item column among {_ITEM_ALIASES}; got {list(df.columns)}"
            )
    return {int(x) for x in pd.to_numeric(df[col], errors="coerce").dropna().astype(int)}


def load_gt_orders(path: PathLike) -> pd.DataFrame:
    """Load GT orders with at least item_id; order_id optional."""
    p = Path(path)
    if p.suffix.lower() in {".parquet", ".pq"}:
        df = pd.read_parquet(p)
    else:
        df = pd.read_csv(p)
    icol = _pick_col(df, _ITEM_ALIASES)
    if icol is None:
        raise ValueError(f"GT orders need item col among {_ITEM_ALIASES}; got {list(df.columns)}")
    out = pd.DataFrame({"item_id": pd.to_numeric(df[icol], errors="coerce")})
    ocol = _pick_col(df, _ORDER_ALIASES)
    if ocol is not None:
        out["order_id"] = df[ocol]
    out = out.dropna(subset=["item_id"])
    out["item_id"] = out["item_id"].astype(int)
    return out.reset_index(drop=True)


def evaluate_item_subset(
    predicted: Sequence[int],
    gt_items: Iterable[int],
    *,
    ks: Sequence[int] = (50, 100, 200),
) -> Dict[str, float]:
    """Hit / precision / recall of localized items vs GT item set."""
    pred = [int(x) for x in predicted]
    gt = {int(x) for x in gt_items}
    n_gt = len(gt)
    out: Dict[str, float] = {
        "n_pred": float(len(pred)),
        "n_gt": float(n_gt),
        "n_hit_all": float(len(set(pred) & gt)),
    }
    if n_gt == 0:
        out["recall_all"] = float("nan")
        out["precision_all"] = float("nan")
        for k in ks:
            out[f"hit@{k}"] = float("nan")
            out[f"precision@{k}"] = float("nan")
            out[f"recall@{k}"] = float("nan")
        return out

    hit_all = set(pred) & gt
    out["recall_all"] = float(len(hit_all) / n_gt)
    out["precision_all"] = float(len(hit_all) / max(len(pred), 1))

    for k in ks:
        kk = min(int(k), len(pred))
        top = set(pred[:kk]) if kk else set()
        inter = top & gt
        out[f"hit@{k}"] = float(len(inter))
        out[f"precision@{k}"] = float(len(inter) / max(kk, 1))
        out[f"recall@{k}"] = float(len(inter) / n_gt)
    return out


def evaluate_order_coverage(
    predicted: Sequence[int],
    gt_orders: pd.DataFrame,
    *,
    ks: Sequence[int] = (50, 100, 200),
) -> Dict[str, float]:
    """Fraction of GT orders that touch ≥1 localized item."""
    if "order_id" not in gt_orders.columns:
        return {"order_coverage": float("nan"), "n_orders": float("nan")}
    pred_set_full = set(int(x) for x in predicted)
    orders = gt_orders.dropna(subset=["order_id", "item_id"]).copy()
    orders["item_id"] = orders["item_id"].astype(int)
    n_orders = int(orders["order_id"].nunique())
    out: Dict[str, float] = {"n_orders": float(n_orders)}
    if n_orders == 0:
        out["order_coverage"] = float("nan")
        return out

    def _cov(pred_set: Set[int]) -> float:
        touched = (
            orders[orders["item_id"].isin(pred_set)]["order_id"].nunique() if pred_set else 0
        )
        return float(touched / n_orders)

    out["order_coverage"] = _cov(pred_set_full)
    for k in ks:
        kk = min(int(k), len(predicted))
        out[f"order_coverage@{k}"] = _cov(set(int(x) for x in predicted[:kk]))
    return out


def evaluate_gt(
    predicted_items: Sequence[int],
    *,
    gt_items_path: Optional[PathLike] = None,
    gt_orders_path: Optional[PathLike] = None,
    ks: Sequence[int] = (50, 100, 200),
) -> Dict:
    """One-shot evaluator: optional item list and/or order list."""
    blob: Dict = {"ks": list(ks), "available": False}
    pred = [int(x) for x in predicted_items]
    gt_items: Set[int] = set()
    if gt_items_path:
        gt_items |= load_gt_items(gt_items_path)
    orders_df = None
    if gt_orders_path:
        orders_df = load_gt_orders(gt_orders_path)
        gt_items |= set(orders_df["item_id"].astype(int).tolist())
        blob["orders"] = evaluate_order_coverage(pred, orders_df, ks=ks)
    if not gt_items:
        blob["note"] = "no GT items/orders provided"
        return blob
    blob["available"] = True
    blob["items"] = evaluate_item_subset(pred, gt_items, ks=ks)
    return blob
