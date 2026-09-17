"""Layer freeze is only for a large PO-risk deviation that MSE also confirms.

PO × MSE contrast on the board:
  both quiet → keep training every layer
  PO broken, MSE holds → watch (do not freeze yet)
  both broken → freeze (this update strategy is not working)

Tabular PO-risk keeps its own outcome model μ(Y|X) and propensity e(T|X).
No online-bootstrap. Causal MA of both series is the stability readout.

When we do freeze, keep PO-risk and MSE together per freeze-depth:

    PO_Dict  = {layer0: array, layer1: array, …, layer_k: array}
    MSE_Dict = {layer0: array, layer1: array, …, layer_k: array}
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
from streaming_po_risk import (
    ACTION_FREEZE,
    REF_N,
    annotate_moving_average,
    annotate_po_mse_contrast,
    batch_mse,
    large_deviation,
    ma_window,
    moving_average,
    po_mse_action,
    ref_split_baseline,
    streaming_po_and_mse,
    streaming_po_risk,
)


def layer_key(i: int) -> str:
    return f"layer{int(i)}"


def stack_layer_metric_dicts(rows, k: int) -> dict:
    """PO-risk and MSE per freeze-depth, aligned to the stream index.

    layer i = model_i (train the top i groups). Missing hops are NaN —
    clones are not trained on a quiet MA.
    """
    n = len(rows)
    keys = [layer_key(i) for i in range(int(k) + 1)]
    po = {name: np.full(n, np.nan) for name in keys}
    mse = {name: np.full(n, np.nan) for name in keys}
    for t, rec in enumerate(rows):
        for x in rec.get("layers") or []:
            name = layer_key(x["i"])
            if name not in po:
                continue
            if x.get("po_fit") is not None:
                po[name][t] = float(x["po_fit"])
            if x.get("mse") is not None:
                mse[name][t] = float(x["mse"])
    return {"PO_Dict": po, "MSE_Dict": mse}


def attach_layer_dicts(result: dict) -> dict:
    """Replay-safe: rebuild PO_Dict / MSE_Dict and the PO×MSE contrast."""
    k = int(result["k"])
    result.update(stack_layer_metric_dicts(result.get("rows") or [], k))
    n_new = result.get("n_new")
    rows = result.get("rows") or []
    if n_new is None and rows:
        n_new = rows[0].get("n_new")
    if rows and "po_base" in result:
        result.update(
            annotate_po_mse_contrast(
                rows, result["po_base"], mse_base=result.get("mse_base"), n_new=n_new
            )
        )
    return result


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


def _refresh_clones(full: AnyMLP):
    return spawn_layer_models(full)


def run_deviation_gate(
    X,
    Y,
    *,
    n_ref: int = REF_N,
    batch_size_stream: int = 5000,
    max_batches: int | None = None,
    seed: int = 2026,
    n_ref_eval: int | None = None,
):
    """PO-risk gate only. Does not decide when to update — that is business logic."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    X_ref, Y_ref = X[:n_ref], Y[:n_ref]
    X_stream, Y_stream = X[n_ref:], Y[n_ref:]
    take = int(n_ref if n_ref_eval is None else min(n_ref_eval, n_ref))
    rng = np.random.default_rng(seed + 7)
    ref_eval = np.arange(n_ref) if take == n_ref else rng.choice(n_ref, size=take, replace=False)
    po_base = ref_split_baseline(X_ref[ref_eval], Y_ref[ref_eval], seed=seed)
    n_batches = len(Y_stream) // batch_size_stream
    if max_batches is not None:
        n_batches = min(n_batches, int(max_batches))
    rows = []
    cursor = 0
    for t in range(n_batches):
        sl = slice(cursor, cursor + batch_size_stream)
        cursor += batch_size_stream
        po_stream = streaming_po_risk(
            X_ref[ref_eval], Y_ref[ref_eval], X_stream[sl], Y_stream[sl], mu_fn=None, seed=seed + t
        )
        large = large_deviation(po_stream, po_base)
        rows.append(
            {
                "t": t,
                "n_new": int(batch_size_stream),
                "po_stream": float(po_stream),
                "po_base": float(po_base),
                "large_deviation": bool(large),
            }
        )
    ma_info = annotate_moving_average(rows, po_base, n_new=batch_size_stream)
    return {
        "n_ref": int(n_ref),
        "n_new": int(batch_size_stream),
        "n_batches": int(n_batches),
        "po_base": float(po_base),
        "frac_large": float(np.mean([r["large_deviation"] for r in rows])) if rows else 0.0,
        "po_mean": float(np.mean([r["po_stream"] for r in rows])) if rows else 0.0,
        "po_std": float(np.std([r["po_stream"] for r in rows])) if rows else 0.0,
        "rows": rows,
        **ma_info,
    }


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
    """Stream T=1 batches. Freeze only when PO-risk and MSE both break."""
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    if len(Y) < n_ref + 50:
        raise ValueError(f"need n_ref={n_ref} plus a stream, got n={len(Y)}")
    X_ref, Y_ref = X[:n_ref], Y[:n_ref]
    X_stream, Y_stream = X[n_ref:], Y[n_ref:]
    full, registry = pretrain_mlp(
        X_ref, Y_ref, hidden_dims=hidden_dims, registry=registry, batch_size=loader_batch, seed=seed
    )
    apply_train_top_i(full, full.n_layer_groups)
    k = full.n_layer_groups
    take = int(n_ref if n_ref_eval is None else min(n_ref_eval, n_ref))
    rng = np.random.default_rng(seed + 7)
    ref_eval = np.arange(n_ref) if take == n_ref else rng.choice(n_ref, size=take, replace=False)
    po_base = ref_split_baseline(X_ref[ref_eval], Y_ref[ref_eval], seed=seed)
    serve = _mu_fn(registry, full)
    mse_base = batch_mse(Y_ref[ref_eval], serve(X_ref[ref_eval]))
    w = ma_window(batch_size_stream)

    rows = []
    po_hist, mse_hist = [], []
    n_stream = len(Y_stream)
    n_batches = n_stream // batch_size_stream
    if max_batches is not None:
        n_batches = min(n_batches, int(max_batches))
    cursor = 0
    for t in range(n_batches):
        sl = slice(cursor, cursor + batch_size_stream)
        cursor += batch_size_stream
        Xb, Yb = X_stream[sl], Y_stream[sl]
        po_stream = streaming_po_risk(X_ref[ref_eval], Y_ref[ref_eval], Xb, Yb, mu_fn=None, seed=seed + t)
        mse_stream = batch_mse(Yb, serve(Xb))
        po_hist.append(float(po_stream))
        mse_hist.append(float(mse_stream))
        po_ma = float(moving_average(po_hist, w)[-1])
        mse_ma = float(moving_average(mse_hist, w)[-1])
        po_broken = large_deviation(po_ma, po_base)
        mse_broken = large_deviation(mse_ma, mse_base)
        action = po_mse_action(po_broken, mse_broken)
        rec = {
            "t": t,
            "n_new": int(len(Yb)),
            "n_ref": int(take),
            "po_stream": float(po_stream),
            "po_base": float(po_base),
            "po_ma": float(po_ma),
            "mse_stream": float(mse_stream),
            "mse_base": float(mse_base),
            "mse_ma": float(mse_ma),
            "large_deviation": bool(large_deviation(po_stream, po_base)),
            "po_broken": bool(po_broken),
            "mse_broken": bool(mse_broken),
            "action": action,
            "all_trainable": action != ACTION_FREEZE,
            "freeze_from": None,
            "i_star": int(k) if action != ACTION_FREEZE else None,
            "layers": [],
        }
        if action != ACTION_FREEZE:
            loader = construct_dataloader(Xb, Yb, batch_size=loader_batch, shuffle=True)
            apply_train_top_i(full, k)
            registry._train_loop(
                full, loader, val_loader=None, task="classification", epochs=online_epochs, restore_best=False
            )
            rec["i_star"] = int(k)
            rec["freeze_from"] = None
        else:
            clones = _refresh_clones(full)
            star_i, star_po = 0, float("inf")
            star_mse_i, star_mse = 0, float("inf")
            layer_rows = []
            for model in clones:
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
                po_fit, mse = streaming_po_and_mse(
                    X_ref[ref_eval], Y_ref[ref_eval], Xb, Yb, mu_fn=mu_fn, seed=seed + t + i
                )
                layer_rows.append(
                    {
                        "i": i,
                        "name": f"model_{i}",
                        "layer": layer_key(i),
                        "po_fit": float(po_fit),
                        "mse": float(mse),
                    }
                )
                if po_fit < star_po:
                    star_po, star_i = po_fit, i
                if mse < star_mse:
                    star_mse, star_mse_i = mse, i
            rec["layers"] = layer_rows
            rec["i_star"] = int(star_i)
            rec["po_star"] = float(star_po)
            rec["i_star_mse"] = int(star_mse_i)
            rec["mse_star"] = float(star_mse)
            rec["freeze_from"] = None if star_i == k else f"model_{star_i}"
            rec["all_trainable"] = star_i == k
            winner = next(m for m in clones if int(m._freeze_i) == star_i)
            full.load_state_dict(winner.state_dict())
        rows.append(rec)
        serve = _mu_fn(registry, full)
    n_large = sum(1 for r in rows if r["large_deviation"])
    freeze_rows = [r for r in rows if r["action"] == ACTION_FREEZE]
    out = {
        "k": int(k),
        "n_ref": int(n_ref),
        "n_new": int(batch_size_stream),
        "n_batches": int(n_batches),
        "hidden_dims": list(hidden_dims),
        "po_base": float(po_base),
        "mse_base": float(mse_base),
        "n_large": int(n_large),
        "rows": rows,
        "layer_names": [f"model_{i}" for i in range(k + 1)],
        "layer_keys": [layer_key(i) for i in range(k + 1)],
        "recommend_i": int(k) if not freeze_rows else int(
            np.round(np.median([r["i_star"] for r in freeze_rows]))
        ),
    }
    out.update(annotate_po_mse_contrast(rows, po_base, mse_base=mse_base, n_new=batch_size_stream))
    out.update(stack_layer_metric_dicts(rows, k))
    return out
