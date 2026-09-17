"""Online layer-freeze CV: read PO-risk to see where to start updating.

Pretrain AnyMLP on D_ref (T=0, n_ref=10000). Clone k+1 models.
On each incoming batch (T=1): CosineAnnealingLR from the top i layers,
then read PO-risk on (ref ∪ new). i* = argmin_i PO-risk is the layer
to start updating from. Incoming n_new must be large; a small batch
would need online-bootstrap, which is too expensive for the board.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from dl_model_registry import (
    AnyMLP,
    DLModelRegistry,
    apply_train_top_i,
    construct_dataloader,
    spawn_layer_models,
)
from streaming_po_risk import REF_N, streaming_po_risk


def _brier(y, p) -> float:
    y = np.asarray(y, dtype=float).ravel()
    p = np.asarray(p, dtype=float).ravel()
    return float(np.mean((p - y) ** 2))


def _mu_fn(registry: DLModelRegistry, model, task="classification") -> Callable:
    def fn(X):
        pred = registry.predict_numpy(model, X, task=task)
        pred = np.asarray(pred, dtype=float)
        if pred.ndim > 1:
            if pred.shape[-1] == 1:
                pred = pred.reshape(-1)
            else:
                pred = pred[:, 1] if pred.shape[-1] == 2 else pred.max(axis=1)
        return pred

    return fn


def pretrain_mlp(
    X_ref,
    Y_ref,
    *,
    hidden_dims=(64, 32),
    registry: DLModelRegistry | None = None,
    batch_size: int = 128,
    val_frac: float = 0.15,
    seed: int = 2026,
) -> tuple[AnyMLP, DLModelRegistry]:
    rng = np.random.default_rng(seed)
    n = len(Y_ref)
    idx = rng.permutation(n)
    n_val = max(1, int(n * val_frac))
    va, tr = idx[:n_val], idx[n_val:]
    registry = registry or DLModelRegistry(epochs=6, patience=3, scheduler="cosine", verbose=False)
    adapter = registry.make_any_mlp(input_dim=X_ref.shape[1], output_dim=1, hidden_dims=hidden_dims)
    model = adapter.build().to(registry.device)
    train_loader = construct_dataloader(X_ref[tr], Y_ref[tr], batch_size=batch_size, shuffle=True)
    val_loader = construct_dataloader(X_ref[va], Y_ref[va], batch_size=batch_size, shuffle=False)
    model = adapter.fit(model, train_loader, val_loader)
    return model, registry


def run_layer_freeze_cv(
    X,
    Y,
    *,
    n_ref: int = REF_N,
    batch_size_stream: int = 2500,
    hidden_dims=(64, 32),
    online_epochs: int = 2,
    loader_batch: int = 128,
    seed: int = 2026,
    max_batches: int | None = 12,
    registry: DLModelRegistry | None = None,
    n_ref_eval: int | None = None,
):
    """Stream T=1 batches; return per-layer PO-risk rows plus i*(t)."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    if len(Y) < n_ref + 50:
        raise ValueError(f"need n_ref={n_ref} plus a stream, got n={len(Y)}")
    X_ref, Y_ref = X[:n_ref], Y[:n_ref]
    X_stream, Y_stream = X[n_ref:], Y[n_ref:]
    pretrained, registry = pretrain_mlp(
        X_ref, Y_ref, hidden_dims=hidden_dims, registry=registry, batch_size=loader_batch, seed=seed
    )
    models = spawn_layer_models(pretrained)
    k = pretrained.n_layer_groups
    take = int(n_ref if n_ref_eval is None else min(n_ref_eval, n_ref))
    rng = np.random.default_rng(seed + 7)
    ref_eval = np.arange(n_ref) if take == n_ref else rng.choice(n_ref, size=take, replace=False)

    rows = []
    n_stream = len(Y_stream)
    n_batches = n_stream // batch_size_stream
    if max_batches is not None:
        n_batches = min(n_batches, int(max_batches))
    cursor = 0
    for t in range(n_batches):
        sl = slice(cursor, cursor + batch_size_stream)
        cursor += batch_size_stream
        Xb, Yb = X_stream[sl], Y_stream[sl]
        po_data = streaming_po_risk(X_ref[ref_eval], Y_ref[ref_eval], Xb, Yb, mu_fn=None)
        star_i, star_po = 0, float("inf")
        for model in models:
            i = int(model._freeze_i)
            apply_train_top_i(model, i)
            if i > 0:
                loader = construct_dataloader(Xb, Yb, batch_size=loader_batch, shuffle=True)
                registry._train_loop(
                    model,
                    loader,
                    val_loader=None,
                    task="classification",
                    epochs=online_epochs,
                    restore_best=False,
                )
            mu_fn = _mu_fn(registry, model)
            po_fit = streaming_po_risk(X_ref[ref_eval], Y_ref[ref_eval], Xb, Yb, mu_fn=mu_fn)
            p_new = mu_fn(Xb)
            p_ref = mu_fn(X_ref[ref_eval])
            rec = {
                "t": t,
                "i": i,
                "name": f"model_{i}",
                "n_new": int(len(Yb)),
                "n_ref": int(take),
                "po_stream": float(po_data),
                "po_fit": float(po_fit),
                "brier_new": _brier(Yb, p_new),
                "brier_ref": _brier(Y_ref[ref_eval], p_ref),
                "mean_Y_new": float(np.mean(Yb)),
                "mean_Y_ref": float(np.mean(Y_ref[ref_eval])),
            }
            rows.append(rec)
            if po_fit < star_po:
                star_po, star_i = po_fit, i
        for rec in rows[-len(models) :]:
            rec["i_star"] = int(star_i)
            rec["po_star"] = float(star_po)
    return {
        "k": int(k),
        "n_ref": int(n_ref),
        "n_batches": int(n_batches),
        "hidden_dims": list(hidden_dims),
        "rows": rows,
        "layer_names": [f"model_{i}" for i in range(k + 1)],
        "recommend_i": int(np.round(np.median([r["i_star"] for r in rows[:: k + 1]]))) if rows else 0,
    }
