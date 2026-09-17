#!/usr/bin/env python3
"""Grad-OnlineRFPerm monitor: one unfrozen ``||∇L||_2`` stream vs MSE/MMD/PO.

Engineering口径 (no cross-param multiple testing):
  g_t = ||∇_{θ_U} L||_2   # ℓ₂ of *all* unfrozen grads as one vector
  one OnlineRFPerm on T_t = g_t - e_ref
  per-layer shares = diagnostics only (freeze-depth), not separate tests

Datasets (Sep17 MVP): synthetic, Covertype, bank-marketing, electricity,
eeg-eye-state.

  PYTHONPATH=. python3 scripts/run_agod_grad_rfperm_monitor.py \\
    --datasets synthetic covertype bank electricity eeg \\
    --seeds 0 1 2 3 4 --batch-size 128 --n-batches 48 --n-burn 8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.datasets import fetch_covtype, fetch_openml
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import LabelEncoder, StandardScaler

from agod.grad_rfperm import (
    first_reject_index,
    init_grad_rfperm,
    layer_grad_norms,
    lead_time,
    relative_grad_shares,
    top_share_layers,
    unfrozen_grad_l2,
    update_grad_rfperm,
)
from agod.online_rfperm import online_fdr_step, rank_pvalue, update_online_rfperm

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class StreamMLP(nn.Module):
    """Named Linear layers so ``layer_grad_norms`` sees fc1/fc2/fc3."""

    def __init__(self, d: int, hidden: Sequence[int] = (64, 32)):
        super().__init__()
        h1, h2 = hidden
        self.fc1 = nn.Linear(d, h1)
        self.fc2 = nn.Linear(h1, h2)
        self.fc3 = nn.Linear(h2, 1)
        self.act = nn.ReLU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.fc3(self.act(self.fc2(self.act(self.fc1(x))))).squeeze(-1)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def _to_xy_binary_or_reg(X, y, *, max_n: int, seed: int) -> Tuple[np.ndarray, np.ndarray, dict]:
    rng = np.random.default_rng(seed)
    X = np.asarray(X, dtype=np.float64)
    # encode categoricals if object/string
    if hasattr(X, "dtype") and X.dtype == object:
        raise TypeError("pass numeric X")
    y = np.asarray(y)
    if y.dtype.kind in "OSU":
        y = LabelEncoder().fit_transform(y).astype(np.float64)
    else:
        y = y.astype(np.float64)
        # covertype-style multi-class → keep as float label (MSE break still works)
        if np.unique(y).size > 20:
            # rare; binary-ize by median
            y = (y > np.median(y)).astype(np.float64)
    n = min(max_n, len(X))
    # preserve time order when present; else shuffle once then take prefix
    # OpenML packs here are not strictly time-ordered for bank/eeg; electricity is.
    idx = np.arange(len(X))
    if len(X) > n:
        # contiguous window for stream-like drift (prefer later half + earlier)
        start = int(rng.integers(0, max(1, len(X) - n)))
        idx = idx[start : start + n]
    X = X[idx].astype(np.float32)
    y = y[idx].astype(np.float64)
    # drop non-finite
    mask = np.isfinite(X).all(axis=1) & np.isfinite(y)
    X, y = X[mask], y[mask]
    X = StandardScaler().fit_transform(X).astype(np.float32)
    meta = {"n": int(len(X)), "d": int(X.shape[1]), "y_unique": int(np.unique(y).size)}
    return X, y, meta


def load_synthetic(
    max_n: int,
    seed: int,
    *,
    shift_batch: int = 20,
    batch_size: int = 128,
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Controlled concept+covariate drift starting at stream batch ``shift_batch``."""
    rng = np.random.default_rng(seed)
    d = 16
    n = max_n
    X = rng.normal(size=(n, d)).astype(np.float32)
    y = np.zeros(n, dtype=np.float64)
    t_cut = int(shift_batch * batch_size)
    t_cut = max(batch_size * 4, min(t_cut, n - batch_size * 4))
    w0 = rng.normal(size=d)
    w0 /= np.linalg.norm(w0) + 1e-9
    w1 = rng.normal(size=d)
    w1 /= np.linalg.norm(w1) + 1e-9
    # pre: stable linear; post: rotated + amplified + covariate shift
    y[:t_cut] = X[:t_cut] @ w0 + rng.normal(0, 0.2, t_cut)
    X[t_cut:] = X[t_cut:] * 1.6 + 0.75
    y[t_cut:] = X[t_cut:] @ w1 * 2.2 + rng.normal(0, 0.45, n - t_cut)
    # scale using pre-shift stats only (avoid leaking post regime into scaler)
    mu = X[:t_cut].mean(axis=0, keepdims=True)
    sd = X[:t_cut].std(axis=0, keepdims=True) + 1e-6
    X = ((X - mu) / sd).astype(np.float32)
    meta = {
        "n": n,
        "d": d,
        "y_unique": -1,
        "shift_batch": int(t_cut // batch_size),
        "shift_idx": t_cut,
        "name": "synthetic",
    }
    return X, y, meta


def load_covertype(max_n: int, seed: int) -> Tuple[np.ndarray, np.ndarray, dict]:
    bun = fetch_covtype()
    # binary: Spruce/Fir (1) vs rest — classic streaming OOD proxy
    y = (bun.target == 1).astype(np.float64)
    X, y, meta = _to_xy_binary_or_reg(bun.data, y, max_n=max_n, seed=seed)
    meta["name"] = "covertype"
    return X, y, meta


def load_openml_named(name: str, max_n: int, seed: int) -> Tuple[np.ndarray, np.ndarray, dict]:
    bun = fetch_openml(name, version=1, as_frame=True, parser="auto")
    df = bun.data.copy()
    # numericize
    for c in df.columns:
        if df[c].dtype.kind in "OSUb":
            df[c] = LabelEncoder().fit_transform(df[c].astype(str).fillna("NA"))
        else:
            df[c] = df[c].astype(float)
    X = df.to_numpy(dtype=np.float64)
    y = bun.target
    X, y, meta = _to_xy_binary_or_reg(X, y, max_n=max_n, seed=seed)
    meta["name"] = name.replace("-", "_")
    return X, y, meta


def _loaders(batch_size: int):
    return {
        "synthetic": lambda max_n, seed: load_synthetic(
            max_n, seed, shift_batch=20, batch_size=batch_size
        ),
        "covertype": load_covertype,
        "bank": lambda max_n, seed: load_openml_named("bank-marketing", max_n, seed),
        "electricity": lambda max_n, seed: load_openml_named("electricity", max_n, seed),
        "eeg": lambda max_n, seed: load_openml_named("eeg-eye-state", max_n, seed),
    }


def make_stream(X, y, bs: int, n_batches: int) -> List[Tuple[np.ndarray, np.ndarray]]:
    need = bs * n_batches
    X, y = X[:need], y[:need]
    if len(X) < need:
        raise ValueError(f"need {need} rows, got {len(X)}")
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, need, bs)]


# ---------------------------------------------------------------------------
# Metrics helpers
# ---------------------------------------------------------------------------


def rbf_mmd2(X0: np.ndarray, X1: np.ndarray, *, max_n: int = 256, rng: Optional[np.random.Generator] = None) -> float:
    rng = rng or np.random.default_rng(0)
    X0 = np.asarray(X0, np.float64)
    X1 = np.asarray(X1, np.float64)
    if len(X0) > max_n:
        X0 = X0[rng.choice(len(X0), max_n, replace=False)]
    if len(X1) > max_n:
        X1 = X1[rng.choice(len(X1), max_n, replace=False)]
    Z = np.vstack([X0, X1])
    # median heuristic on subsample
    if len(Z) > 400:
        Zs = Z[rng.choice(len(Z), 400, replace=False)]
    else:
        Zs = Z
    d2 = np.sum((Zs[:, None, :] - Zs[None, :, :]) ** 2, axis=-1)
    med = float(np.median(d2[d2 > 0])) if np.any(d2 > 0) else 1.0
    gamma = 1.0 / (2.0 * med + 1e-12)

    def k(A, B):
        d = np.sum((A[:, None, :] - B[None, :, :]) ** 2, axis=-1)
        return np.exp(-gamma * d)

    Kxx = k(X0, X0)
    Kyy = k(X1, X1)
    Kxy = k(X0, X1)
    n, m = len(X0), len(X1)
    # unbiased
    np.fill_diagonal(Kxx, 0.0)
    np.fill_diagonal(Kyy, 0.0)
    return float(
        Kxx.sum() / (n * (n - 1) + 1e-12)
        + Kyy.sum() / (m * (m - 1) + 1e-12)
        - 2.0 * Kxy.mean()
    )


def batch_po(X0, y0, X1, y1, seed: int) -> float:
    if len(X0) < 16 or len(X1) < 16:
        return 0.0
    rf = RandomForestRegressor(
        n_estimators=20, max_depth=4, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    rf.fit(X0, y0)
    e0 = float(np.mean(np.abs(y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(y1 - rf.predict(X1))))
    return max(e1 - e0, 0.0)


class ScalarRFPerm:
    """Lightweight OnlineRFPerm on a scalar series (MMD / PO)."""

    def __init__(self, e_ref: float = 0.0):
        self.e_ref = float(e_ref)
        self.T_hist: List[float] = []
        self.p_hist: List[float] = []
        self.reject_hist: List[int] = []
        self.wealth = 1.0
        self.n_burn = 0

    def step(self, value: float, *, burn_in: bool, alpha: float = 0.05) -> dict:
        T = float(value) - self.e_ref
        out = {"T": T, "p": 1.0, "reject": False, "burn_in": burn_in}
        if burn_in:
            self.T_hist.append(T)
            self.n_burn += 1
            self.p_hist.append(1.0)
            self.reject_hist.append(0)
            # refresh e_ref as burn mean of raw value
            n = len(self.T_hist)
            # recover raw ≈ T + e_ref_old; simpler: set e_ref to running mean of values
            self.e_ref = float(self.e_ref * (n - 1) / n + value / n) if n else float(value)
            # rewrite last T relative to updated ref for consistency next steps
            self.T_hist[-1] = float(value) - self.e_ref
            return out
        p = rank_pvalue(T, self.T_hist, ewma=True)
        rej = online_fdr_step(self, p, alpha=alpha, procedure="alpha_investing")  # type: ignore[arg-type]
        self.T_hist.append(T)
        self.p_hist.append(p)
        out.update({"p": p, "reject": bool(rej)})
        return out


# ---------------------------------------------------------------------------
# Stream run
# ---------------------------------------------------------------------------


def _pretrain_mlp(
    model: StreamMLP,
    opt: torch.optim.Optimizer,
    stream: List[Tuple[np.ndarray, np.ndarray]],
    n_burn: int,
    *,
    steps: int = 80,
) -> None:
    """Fit MLP on burn-in batches so grad norms aren't dominated by cold-start."""
    Xs = np.concatenate([stream[t][0] for t in range(n_burn)], axis=0)
    ys = np.concatenate([stream[t][1] for t in range(n_burn)], axis=0)
    xt = torch.from_numpy(np.asarray(Xs, np.float32))
    yt = torch.from_numpy(np.asarray(ys, np.float32))
    model.train()
    for _ in range(steps):
        loss = ((model(xt) - yt) ** 2).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()


@torch.no_grad()
def _clear_grads(model: nn.Module) -> None:
    for p in model.parameters():
        p.grad = None


def serving_grad_snapshot(model: nn.Module, X: np.ndarray, y: np.ndarray) -> dict:
    """One backward on frozen f_ref → global unfrozen ℓ₂ + layer shares (diag)."""
    model.train()
    _clear_grads(model)
    xt = torch.from_numpy(np.asarray(X, np.float32))
    yt = torch.from_numpy(np.asarray(y, np.float32))
    loss = ((model(xt) - yt) ** 2).mean()
    loss.backward()
    g = unfrozen_grad_l2(model)
    norms = layer_grad_norms(model, unfrozen_only=True)
    shares = relative_grad_shares(norms) if norms else {}
    _clear_grads(model)
    model.eval()
    return {"g": g, "layer_norms": norms, "shares": shares, "loss": float(loss.detach().item())}


def run_dataset(
    name: str,
    stream: List[Tuple[np.ndarray, np.ndarray]],
    *,
    n_burn: int,
    alpha: float,
    seed: int,
    train_steps: int,
    lr: float,
    use_share: bool,
) -> dict:
    del use_share  # shares are diagnostic only; kept in API for CLI compat
    torch.manual_seed(seed)
    d = stream[0][0].shape[1]
    # Frozen serving MLP = Grad OnlineRFPerm's f_ref (weights fixed after pretrain)
    f_ref = StreamMLP(d)
    opt_ref = torch.optim.Adam(f_ref.parameters(), lr=lr)
    _pretrain_mlp(f_ref, opt_ref, stream, max(n_burn, 1), steps=max(80, 10 * train_steps))
    for p in f_ref.parameters():
        p.requires_grad_(True)  # need grads for monitoring; never stepped after pretrain

    layer_names = ["fc1", "fc2", "fc3"]
    # ONE OnlineRFPerm stream on global unfrozen ||∇||_2 — no per-layer tests
    grad_state = init_grad_rfperm("unfrozen_l2")

    X0, y0 = stream[0]
    from agod.online_rfperm import fit_online_rfperm

    mse_state = fit_online_rfperm(X0, y0, seed=seed)
    mmd_mon = ScalarRFPerm(e_ref=0.0)
    po_mon = ScalarRFPerm(e_ref=0.0)

    mse_traj, mmd_traj, po_traj = [], [], []
    mse_p, mmd_p, po_p = [], [], []
    mse_rej, mmd_rej, po_rej = [], [], []
    grad_p: List[float] = []
    grad_g: List[float] = []
    grad_rej: List[int] = []
    share_traj: Dict[str, List[float]] = {n: [] for n in layer_names}
    layer_g: Dict[str, List[float]] = {n: [] for n in layer_names}
    train_mse: List[float] = []
    top_share_at_reject: Optional[List[Tuple[str, float]]] = None

    Xp, yp = X0, y0
    for t in range(len(stream)):
        Xc, yc = stream[t]
        burn = t < n_burn

        snap = serving_grad_snapshot(f_ref, Xc, yc)
        train_mse.append(float(snap["loss"]))
        ug = update_grad_rfperm(grad_state, float(snap["g"]), burn_in=burn, alpha=alpha)
        grad_p.append(float(ug["p"]))
        grad_g.append(float(snap["g"]))
        grad_rej.append(int(ug["reject"]))

        shares = snap["shares"]
        norms = snap["layer_norms"]
        for n in layer_names:
            share_traj[n].append(float(shares.get(n, 0.0)))
            layer_g[n].append(float(norms.get(n, 0.0)))

        if ug["reject"] and top_share_at_reject is None and not burn:
            top_share_at_reject = top_share_layers(shares, k=3)

        step_mse = update_online_rfperm(mse_state, Xc, yc, burn_in=burn, alpha=alpha)
        mse_traj.append(float(step_mse["T"]))
        mse_p.append(float(step_mse["p"]))
        mse_rej.append(int(step_mse["reject"]))

        if t == 0:
            mmd_v, po_v = 0.0, 0.0
        else:
            mmd_v = rbf_mmd2(Xp, Xc, rng=np.random.default_rng(seed + t))
            po_v = batch_po(Xp, yp, Xc, yc, seed + 17 * t)
        sm = mmd_mon.step(mmd_v, burn_in=burn, alpha=alpha)
        sp = po_mon.step(po_v, burn_in=burn, alpha=alpha)
        mmd_traj.append(float(sm["T"]))
        mmd_p.append(float(sm["p"]))
        mmd_rej.append(int(sm["reject"]))
        po_traj.append(float(sp["T"]))
        po_p.append(float(sp["p"]))
        po_rej.append(int(sp["reject"]))

        Xp, yp = Xc, yc

    after = n_burn
    t_grad = first_reject_index(grad_state.reject_hist, after=after)
    t_mse = first_reject_index(mse_state.reject_hist, after=after)
    t_mmd = first_reject_index(mmd_mon.reject_hist, after=after)
    t_po = first_reject_index(po_mon.reject_hist, after=after)

    return {
        "dataset": name,
        "n_batches": len(stream),
        "n_burn": n_burn,
        "alpha": alpha,
        "protocol": "unfrozen_l2_single_stream",
        "first_reject": {
            "grad": t_grad,
            "mse": t_mse,
            "mmd": t_mmd,
            "po": t_po,
        },
        "lead_time": {
            "grad_minus_mse": lead_time(t_grad, t_mse),
            "grad_minus_mmd": lead_time(t_grad, t_mmd),
            "grad_minus_po": lead_time(t_grad, t_po),
        },
        "top_share_at_reject": top_share_at_reject,
        "_states_for_shift": {
            "grad": list(grad_state.reject_hist),
            "mse": list(mse_state.reject_hist),
            "mmd": list(mmd_mon.reject_hist),
            "po": list(po_mon.reject_hist),
        },
        "duty": {
            "grad": float(np.mean(grad_state.reject_hist[after:]))
            if after < len(grad_state.reject_hist)
            else 0.0,
            "mse": float(np.mean(mse_state.reject_hist[after:]))
            if after < len(mse_state.reject_hist)
            else 0.0,
            "mmd": float(np.mean(mmd_mon.reject_hist[after:]))
            if after < len(mmd_mon.reject_hist)
            else 0.0,
            "po": float(np.mean(po_mon.reject_hist[after:]))
            if after < len(po_mon.reject_hist)
            else 0.0,
        },
        "traj": {
            "train_mse": train_mse,
            "mse_T": mse_traj,
            "mse_p": mse_p,
            "mse_rej": mse_rej,
            "mmd_T": mmd_traj,
            "mmd_p": mmd_p,
            "mmd_rej": mmd_rej,
            "po_T": po_traj,
            "po_p": po_p,
            "po_rej": po_rej,
            "grad_p": grad_p,
            "grad_g": grad_g,
            "grad_rej": grad_rej,
            "share": share_traj,
            "layer_g": layer_g,
        },
    }


def plot_dataset(res: dict, out_dir: Path) -> Path:
    name = res["dataset"]
    traj = res["traj"]
    n_burn = res["n_burn"]
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.5), sharex=True)

    ax = axes[0]
    ax.plot(traj["grad_p"], label="Grad p (unfrozen ℓ₂)", color="#BF616A", lw=1.6)
    ax.plot(traj["mse_p"], label="MSE p", color="#5E81AC", lw=1.6, ls="--")
    ax.plot(traj["mmd_p"], label="MMD p", color="#B48EAD", lw=1.2, ls=":")
    ax.plot(traj["po_p"], label="PO p", color="#D08770", lw=1.2, ls="-.")
    ax.axvline(n_burn - 0.5, color="gray", ls="--", lw=0.9, label="end burn")
    fr = res["first_reject"]
    if fr.get("grad") is not None:
        ax.axvline(fr["grad"], color="#BF616A", lw=1.0, alpha=0.7)
    if fr.get("mse") is not None:
        ax.axvline(fr["mse"], color="#5E81AC", lw=1.0, alpha=0.7)
    ax.set_ylabel("p-value")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(ncol=3, fontsize=7, loc="upper right")
    ax.set_title(f"{name}: single-stream Grad-OnlineRFPerm vs MSE / MMD / PO")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    g = np.asarray(traj["grad_g"], float)
    ax.plot(g / (g.max() + 1e-12), label="‖∇_U‖₂ (rel)", color="#BF616A", lw=1.5)
    for n, color in zip(["fc1", "fc2", "fc3"], ["#EBCB8B", "#A3BE8C", "#88C0D0"]):
        s = np.asarray(traj["share"][n], float)
        ax.plot(s, label=f"share {n}", color=color, lw=1.1, alpha=0.85)
    ax.plot(np.asarray(traj["train_mse"], float), label="serving MSE", color="#4C566A", lw=1.0)
    ax.set_ylabel("energy / share")
    ax.legend(ncol=4, fontsize=7)
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    for key, color, lab in [
        ("mse_rej", "#5E81AC", "MSE reject"),
        ("mmd_rej", "#B48EAD", "MMD reject"),
        ("po_rej", "#D08770", "PO reject"),
    ]:
        ax.step(range(len(traj[key])), traj[key], where="post", label=lab, color=color, lw=1.2)
    ax.step(
        range(len(traj["grad_rej"])),
        traj["grad_rej"],
        where="post",
        label="Grad reject (1 stream)",
        color="#BF616A",
        lw=1.5,
    )
    ax.set_ylabel("reject")
    ax.set_xlabel("batch t")
    ax.set_ylim(-0.05, 1.15)
    lt = res["lead_time"]["grad_minus_mse"]
    ax.set_title(f"rejects  |  lead(grad−mse)={lt}")
    ax.legend(ncol=4, fontsize=7)
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    path = out_dir / f"{name}_grad_rfperm.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def plot_lead_summary(all_res: Dict[str, dict], out_dir: Path) -> Path:
    names = list(all_res.keys())
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    xs = np.arange(len(names))
    vals = []
    for n in names:
        lt = all_res[n]["lead_time"]["grad_minus_mse"]
        vals.append(0.0 if lt is None else float(lt))
    colors = ["#A3BE8C" if v < 0 else ("#EBCB8B" if v == 0 else "#BF616A") for v in vals]
    ax.bar(xs, vals, color=colors, width=0.55)
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set_xticks(xs)
    ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel("lead time  t_grad − t_mse  (neg = earlier)")
    ax.set_title("Grad-OnlineRFPerm lead vs MSE break")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    path = out_dir / "lead_time_summary.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def write_report(all_res: Dict[str, dict], out_dir: Path, meta: dict) -> Path:
    lines = [
        "# Grad-OnlineRFPerm (`param.grad.norm`) — Sep17 MVP datasets",
        "",
        "**One** OnlineRFPerm stream on unfrozen global grad energy",
        "`g_t = ||∇_{θ_U} L||_2` (not per-layer tests; no cross-param multiplicity).",
        "Layer shares are diagnostics only (freeze-depth after global reject).",
        "",
        "Lead time: `t_grad − t_mse` (negative ⇒ Grad rejects **before** MSE break).",
        "",
        f"- batch_size={meta['batch_size']}, n_batches={meta['n_batches']}, "
        f"n_burn={meta['n_burn']}, alpha={meta['alpha']}, seed={meta['seed']}, "
        f"protocol={meta.get('protocol')}",
        "",
        "| dataset | t_grad | t_mse | t_mmd | t_po | lead(g−mse) | mse duty | top share @ reject |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for name, r in all_res.items():
        fr = r["first_reject"]
        lt = r["lead_time"]

        def _f(x):
            return "—" if x is None else str(x)

        top = r.get("top_share_at_reject") or []
        top_s = ", ".join(f"{k}:{v:.2f}" for k, v in top) if top else "—"
        lines.append(
            f"| `{name}` | {_f(fr['grad'])} | {_f(fr['mse'])} | {_f(fr['mmd'])} | "
            f"{_f(fr['po'])} | {_f(lt['grad_minus_mse'])} | {r['duty']['mse']:.2f} | {top_s} |"
        )
    syn = all_res.get("synthetic")
    if syn and "detection_delay" in syn:
        d = syn["detection_delay"]
        lines += [
            "",
            f"### Synthetic controlled shift @ batch {syn.get('shift_batch')}",
            "",
            f"Detection delay after shift: Grad=`{d.get('grad')}`, MSE=`{d.get('mse')}`, "
            f"MMD=`{d.get('mmd')}`, PO=`{d.get('po')}` "
            f"(post-shift lead grad−mse=`{syn['lead_time'].get('grad_minus_mse_postshift')}`).",
            "",
        ]
    lines += [
        "",
        "## Method (no multiple testing across params)",
        "",
        "```",
        "θ_U = {params with requires_grad=True}",
        "g_t = ||∇_{θ_U} L(batch; f_ref)||_2     # ONE scalar, ℓ₂ of full unfrozen grad",
        "T_t = g_t - e_ref",
        "p_t = rank/EWMA(T_t vs history); alpha-investing FDR → reject  # ONE stream",
        "# NOT: mean of per-layer norms as separate tests",
        "# NOT: OnlineRFPerm per layer then any(reject)",
        "shares_ℓ = ||g_ℓ|| / g_t                 # diagnostic only after reject",
        "```",
        "",
    ]
    path = out_dir / "AGOD_grad_rfperm_monitor.md"
    path.write_text("\n".join(lines) + "\n")
    docs = ROOT / "docs" / "agod" / "AGOD_grad_rfperm_monitor.md"
    docs.parent.mkdir(parents=True, exist_ok=True)
    docs.write_text(path.read_text())
    return path


def _attach_shift_delay(res: dict, meta: dict) -> None:
    t_shift = meta.get("shift_batch")
    if t_shift is None:
        res.pop("_states_for_shift", None)
        return
    st = res.pop("_states_for_shift", {})
    tg = first_reject_index(st.get("grad", []), after=int(t_shift))
    delay = {"grad": None if tg is None else int(tg - t_shift)}
    for key in ("mse", "mmd", "po"):
        t0 = first_reject_index(st.get(key, []), after=int(t_shift))
        delay[key] = None if t0 is None else int(t0 - t_shift)
    res["shift_batch"] = int(t_shift)
    res["detection_delay"] = delay
    res["lead_time"]["grad_minus_mse_postshift"] = lead_time(
        tg, first_reject_index(st.get("mse", []), after=int(t_shift))
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["synthetic", "covertype", "bank", "electricity", "eeg"],
    )
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--n-batches", type=int, default=48)
    ap.add_argument("--n-burn", type=int, default=8)
    ap.add_argument("--max-n", type=int, default=12000)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--seeds",
        type=int,
        nargs="*",
        default=None,
        help="If provided, run multi-seed aggregate instead of single --seed",
    )
    ap.add_argument("--train-steps", type=int, default=12)
    ap.add_argument("--lr", type=float, default=1e-2)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results" / "agod_grad_rfperm")
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.seeds:
        _run_multiseed(args)
        return
    all_res: Dict[str, dict] = {}

    LOADERS = _loaders(args.batch_size)
    for name in args.datasets:
        if name not in LOADERS:
            raise SystemExit(f"unknown dataset {name}; choose from {list(LOADERS)}")
        print(f"=== {name} ===", flush=True)
        X, y, meta = LOADERS[name](args.max_n, args.seed)
        print(f"  loaded n={meta['n']} d={meta['d']} y_unique={meta.get('y_unique')}", flush=True)
        need = args.batch_size * args.n_batches
        if meta["n"] < need:
            # shrink batches
            nb = max(args.n_burn + 4, meta["n"] // args.batch_size)
            print(f"  shrink n_batches {args.n_batches} → {nb}", flush=True)
            n_batches = nb
        else:
            n_batches = args.n_batches
        stream = make_stream(X, y, args.batch_size, n_batches)
        res = run_dataset(
            name,
            stream,
            n_burn=min(args.n_burn, max(2, n_batches // 4)),
            alpha=args.alpha,
            seed=args.seed,
            train_steps=args.train_steps,
            lr=args.lr,
            use_share=True,
        )
        res["meta"] = meta
        _attach_shift_delay(res, meta)
        all_res[name] = res
        fr, lt = res["first_reject"], res["lead_time"]
        extra = ""
        if "detection_delay" in res:
            extra = f" | delay@shift={res['detection_delay']}"
        top = res.get("top_share_at_reject")
        top_s = "" if not top else f" | top_share={top}"
        print(
            f"  first: grad={fr['grad']} mse={fr['mse']} mmd={fr['mmd']} po={fr['po']} "
            f"| lead(g-mse)={lt['grad_minus_mse']}{extra}{top_s}",
            flush=True,
        )
        plot_dataset(res, out_dir)

    plot_lead_summary(all_res, out_dir)
    meta = {
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "n_burn": args.n_burn,
        "alpha": args.alpha,
        "seed": args.seed,
        "train_steps": args.train_steps,
        "protocol": "unfrozen_l2_single_stream",
    }
    report = write_report(all_res, out_dir, meta)
    # compact JSON (drop long traj for summary file; full traj kept per-ds)
    summary = {
        "meta": meta,
        "datasets": {
            k: {
                "meta": v.get("meta"),
                "first_reject": v["first_reject"],
                "lead_time": v["lead_time"],
                "duty": v["duty"],
                **(
                    {
                        "shift_batch": v["shift_batch"],
                        "detection_delay": v["detection_delay"],
                    }
                    if "detection_delay" in v
                    else {}
                ),
            }
            for k, v in all_res.items()
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    # full results (can be large)
    full = {k: {kk: vv for kk, vv in v.items() if kk != "traj"} | {"traj_keys": list(v["traj"].keys())} for k, v in all_res.items()}
    # keep traj in separate files
    for k, v in all_res.items():
        (out_dir / f"{k}_traj.json").write_text(json.dumps(v["traj"]))
    (out_dir / "results_meta.json").write_text(json.dumps(full, indent=2, default=str))
    print(f"wrote {report}", flush=True)
    print(f"summary → {out_dir / 'summary.json'}", flush=True)


def _run_multiseed(args: argparse.Namespace) -> None:
    from collections import defaultdict

    out_dir = args.out_dir
    LOADERS = _loaders(args.batch_size)
    seeds = list(args.seeds)
    agg = {d: defaultdict(list) for d in args.datasets}
    all_runs: Dict[str, dict] = {}
    seed0_res: Dict[str, dict] = {}

    for seed in seeds:
        print(f"===== seed {seed} =====", flush=True)
        all_runs[str(seed)] = {}
        for name in args.datasets:
            X, y, meta = LOADERS[name](args.max_n, seed)
            need = args.batch_size * args.n_batches
            n_batches = (
                args.n_batches
                if meta["n"] >= need
                else max(args.n_burn + 4, meta["n"] // args.batch_size)
            )
            stream = make_stream(X, y, args.batch_size, n_batches)
            res = run_dataset(
                name,
                stream,
                n_burn=min(args.n_burn, max(2, n_batches // 4)),
                alpha=args.alpha,
                seed=seed,
                train_steps=args.train_steps,
                lr=args.lr,
                use_share=True,
            )
            res["meta"] = meta
            _attach_shift_delay(res, meta)
            fr, lt = res["first_reject"], res["lead_time"]
            print(
                f"  {name}: grad={fr['grad']} mse={fr['mse']} lead={lt['grad_minus_mse']}",
                flush=True,
            )
            for k in ("grad_minus_mse", "grad_minus_mmd", "grad_minus_po"):
                v = lt.get(k)
                if v is not None:
                    agg[name][k].append(v)
            if lt.get("grad_minus_mse") is not None:
                agg[name]["grad_earlier"].append(int(lt["grad_minus_mse"] < 0))
            all_runs[str(seed)][name] = {
                "first_reject": fr,
                "lead_time": lt,
                "detection_delay": res.get("detection_delay"),
                "duty": res["duty"],
                "top_share_at_reject": res.get("top_share_at_reject"),
            }
            if seed == seeds[0]:
                seed0_res[name] = res
                plot_dataset(res, out_dir)

    if seed0_res:
        plot_lead_summary(seed0_res, out_dir)

    lines = [
        "# Grad-OnlineRFPerm (`param.grad.norm`) — Sep17 MVP datasets",
        "",
        "**Single stream:** `g_t = ||∇_{θ_U} L||_2` over all unfrozen params,",
        "one OnlineRFPerm (no per-layer / per-param multiple testing).",
        "Layer shares = diagnostics only after global reject.",
        "",
        f"Lead: `t_grad − t_mse` (negative ⇒ Grad earlier). seeds={seeds}, "
        f"batch={args.batch_size}, n_batches={args.n_batches}, n_burn={args.n_burn}, alpha={args.alpha}",
        "",
        "| dataset | mean lead(g−mse) | median | P(earlier) | P(≤0) |",
        "|---|---:|---:|---:|---:|",
    ]
    summary_ds = {}
    for name in args.datasets:
        a = agg[name]
        gm = np.asarray(a["grad_minus_mse"], float)
        pe = float(np.mean(a["grad_earlier"])) if a["grad_earlier"] else float("nan")
        p_le = float(np.mean(gm <= 0)) if len(gm) else float("nan")
        lines.append(
            f"| `{name}` | {gm.mean():+.2f} | {np.median(gm):+.1f} | {pe:.0%} | {p_le:.0%} |"
        )
        summary_ds[name] = {
            "lead_grad_mse_mean": float(gm.mean()),
            "lead_grad_mse_median": float(np.median(gm)),
            "lead_grad_mse_std": float(gm.std()),
            "p_grad_earlier": pe,
            "p_grad_not_later": p_le,
            "leads_grad": [float(x) for x in gm],
        }
    all_g = np.concatenate([np.asarray(agg[d]["grad_minus_mse"], float) for d in args.datasets])
    lines += [
        "",
        f"**Overall** ({len(all_g)} runs): mean lead(g−mse)={all_g.mean():+.2f}, "
        f"P(Grad earlier)={np.mean(all_g < 0):.0%}, P(≤0)={np.mean(all_g <= 0):.0%}.",
        "",
        "## Method (engineering口径)",
        "",
        "```",
        "θ_U = unfrozen params (requires_grad=True)",
        "g_t = ||∇_{θ_U} L||_2 = sqrt(Σ_i ||∇θ_i||²)   # ONE scalar — not mean of norms",
        "OnlineRFPerm once on T_t = g_t - e_ref",
        "shares_ℓ = ||g_ℓ|| / g_t   # diagnostic ranking only, no FDR per layer",
        "```",
        "",
    ]
    text = "\n".join(lines) + "\n"
    (out_dir / "AGOD_grad_rfperm_monitor.md").write_text(text)
    docs = ROOT / "docs" / "agod" / "AGOD_grad_rfperm_monitor.md"
    docs.parent.mkdir(parents=True, exist_ok=True)
    docs.write_text(text)
    payload = {
        "meta": {
            "seeds": seeds,
            "batch_size": args.batch_size,
            "n_batches": args.n_batches,
            "n_burn": args.n_burn,
            "alpha": args.alpha,
            "protocol": "unfrozen_l2_single_stream",
        },
        "by_dataset": summary_ds,
        "overall": {
            "mean_lead_grad_mse": float(all_g.mean()),
            "p_grad_earlier": float(np.mean(all_g < 0)),
            "p_grad_not_later": float(np.mean(all_g <= 0)),
            "n_runs": int(len(all_g)),
        },
        "runs": all_runs,
    }
    (out_dir / "multiseed_summary.json").write_text(json.dumps(payload, indent=2))
    (out_dir / "summary.json").write_text(
        json.dumps({"overall": payload["overall"], "by_dataset": summary_ds}, indent=2)
    )
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    names = list(args.datasets)
    vals = [summary_ds[n]["lead_grad_mse_mean"] for n in names]
    colors = ["#A3BE8C" if v < 0 else ("#EBCB8B" if v == 0 else "#BF616A") for v in vals]
    ax.bar(np.arange(len(names)), vals, color=colors, width=0.55)
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel("mean lead  t_grad − t_mse")
    ax.set_title(f"Grad-OnlineRFPerm mean lead ({len(seeds)} seeds, 1 stream)")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "lead_time_summary.png", dpi=140)
    plt.close(fig)
    print(text, flush=True)


if __name__ == "__main__":
    main()
