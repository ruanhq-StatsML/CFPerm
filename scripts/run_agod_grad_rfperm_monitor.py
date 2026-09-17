#!/usr/bin/env python3
"""Grad-OnlineRFPerm monitor: ``param.grad.norm()`` vs MSE / MMD / PO break.

Continuous-time OnlineRFPerm on per-layer gradient energy (and relative
shares), compared to serving-MSE OnlineRFPerm, consecutive-batch RBF-MMD²,
and window PO-risk — reporting lead time ``t_grad - t_mse``.

Datasets (Sep17 MVP): synthetic, Covertype, bank-marketing, electricity,
eeg-eye-state.

  PYTHONPATH=. python3 scripts/run_agod_grad_rfperm_monitor.py \\
    --datasets synthetic covertype bank electricity eeg \\
    --batch-size 128 --n-batches 48 --n-burn 8
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
    earliest_layer_reject,
    first_reject_index,
    init_grad_rfperm,
    layer_grad_norms,
    lead_time,
    relative_grad_shares,
    update_grad_rfperm,
)
from agod.online_rfperm import OnlineRFPermState, rank_pvalue, update_online_rfperm
from agod.online_rfperm import online_fdr_step

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


def load_synthetic(max_n: int, seed: int, *, shift_at: float = 0.45) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Controlled concept drift: mean/scale change after ``shift_at`` fraction."""
    rng = np.random.default_rng(seed)
    d = 16
    n = max_n
    X = rng.normal(size=(n, d)).astype(np.float32)
    y = np.zeros(n, dtype=np.float64)
    t_cut = int(n * shift_at)
    w0 = rng.normal(size=d)
    w0 /= np.linalg.norm(w0) + 1e-9
    w1 = rng.normal(size=d)
    w1 /= np.linalg.norm(w1) + 1e-9
    # pre: mild linear; post: rotated + amplified + covariate shift
    y[:t_cut] = X[:t_cut] @ w0 + rng.normal(0, 0.25, t_cut)
    X[t_cut:] = X[t_cut:] * 1.35 + rng.normal(0.4, 0.15, size=(n - t_cut, d)).astype(np.float32)
    y[t_cut:] = X[t_cut:] @ w1 * 1.8 + rng.normal(0, 0.55, n - t_cut)
    X = StandardScaler().fit_transform(X).astype(np.float32)
    meta = {
        "n": n,
        "d": d,
        "y_unique": -1,
        "shift_frac": shift_at,
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


LOADERS = {
    "synthetic": load_synthetic,
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
    torch.manual_seed(seed)
    d = stream[0][0].shape[1]
    model = StreamMLP(d)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    layer_names = ["fc1", "fc2", "fc3"]
    grad_states = init_grad_rfperm(layer_names)
    share_states = init_grad_rfperm([f"{n}_share" for n in layer_names])

    # MSE OnlineRFPerm on fixed f_ref (sklearn RF) — serving degradation
    X0, y0 = stream[0]
    from agod.online_rfperm import fit_online_rfperm

    mse_state = fit_online_rfperm(X0, y0, seed=seed)
    mmd_mon = ScalarRFPerm(e_ref=0.0)
    po_mon = ScalarRFPerm(e_ref=0.0)

    mse_traj, mmd_traj, po_traj = [], [], []
    mse_p, mmd_p, po_p = [], [], []
    mse_rej, mmd_rej, po_rej = [], [], []
    grad_p: Dict[str, List[float]] = {n: [] for n in layer_names}
    share_p: Dict[str, List[float]] = {n: [] for n in layer_names}
    grad_g: Dict[str, List[float]] = {n: [] for n in layer_names}
    train_mse: List[float] = []

    Xp, yp = X0, y0
    for t in range(len(stream)):
        Xc, yc = stream[t]
        burn = t < n_burn

        # train one batch; capture grads from last step
        xt = torch.from_numpy(np.asarray(Xc, np.float32))
        yt = torch.from_numpy(np.asarray(yc, np.float32))
        model.train()
        last_norms: Dict[str, float] = {}
        for step_i in range(train_steps):
            pred = model(xt)
            loss = ((pred - yt) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            if step_i == train_steps - 1:
                last_norms = layer_grad_norms(model)
            opt.step()
        train_mse.append(float(loss.detach().item()))

        # relative shares
        shares = relative_grad_shares(last_norms) if last_norms else {n: 0.0 for n in layer_names}

        for n in layer_names:
            g = float(last_norms.get(n, 0.0))
            s = float(shares.get(n, 0.0))
            ug = update_grad_rfperm(grad_states[n], g, burn_in=burn, alpha=alpha)
            us = update_grad_rfperm(
                share_states[f"{n}_share"], s, burn_in=burn, alpha=alpha, use_relative_to_ref=True
            )
            grad_p[n].append(float(ug["p"]))
            share_p[n].append(float(us["p"]))
            grad_g[n].append(g)

        # MSE OnlineRFPerm (skip t=0 as ref already fitted; still burn-compatible)
        step_mse = update_online_rfperm(mse_state, Xc, yc, burn_in=burn, alpha=alpha)
        mse_traj.append(float(step_mse["T"]))
        mse_p.append(float(step_mse["p"]))
        mse_rej.append(int(step_mse["reject"]))

        # MMD / PO vs previous batch
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
    grad_first = earliest_layer_reject(grad_states, after=after)
    share_first = earliest_layer_reject(share_states, after=after)
    t_mse = first_reject_index(mse_state.reject_hist, after=after)
    t_mmd = first_reject_index(mmd_mon.reject_hist, after=after)
    t_po = first_reject_index(po_mon.reject_hist, after=after)
    t_grad = grad_first.get("__any__")
    t_share = share_first.get("__any__")

    return {
        "dataset": name,
        "n_batches": len(stream),
        "n_burn": n_burn,
        "alpha": alpha,
        "use_share": use_share,
        "first_reject": {
            "grad_any": t_grad,
            "grad_share_any": t_share,
            "grad_per_layer": {k: v for k, v in grad_first.items() if k != "__any__"},
            "share_per_layer": {k: v for k, v in share_first.items() if k != "__any__"},
            "mse": t_mse,
            "mmd": t_mmd,
            "po": t_po,
        },
        "lead_time": {
            "grad_minus_mse": lead_time(t_grad, t_mse),
            "grad_minus_mmd": lead_time(t_grad, t_mmd),
            "grad_minus_po": lead_time(t_grad, t_po),
            "share_minus_mse": lead_time(t_share, t_mse),
        },
        "duty": {
            "grad": {
                n: float(np.mean(grad_states[n].reject_hist[after:])) if after < len(grad_states[n].reject_hist) else 0.0
                for n in layer_names
            },
            "mse": float(np.mean(mse_state.reject_hist[after:])) if after < len(mse_state.reject_hist) else 0.0,
            "mmd": float(np.mean(mmd_mon.reject_hist[after:])) if after < len(mmd_mon.reject_hist) else 0.0,
            "po": float(np.mean(po_mon.reject_hist[after:])) if after < len(po_mon.reject_hist) else 0.0,
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
            "share_p": share_p,
            "grad_g": grad_g,
            "grad_rej": {n: list(grad_states[n].reject_hist) for n in layer_names},
        },
    }


def plot_dataset(res: dict, out_dir: Path) -> Path:
    name = res["dataset"]
    traj = res["traj"]
    n_burn = res["n_burn"]
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.5), sharex=True)

    ax = axes[0]
    for n, color in zip(["fc1", "fc2", "fc3"], ["#BF616A", "#EBCB8B", "#A3BE8C"]):
        ax.plot(traj["grad_p"][n], label=f"grad p {n}", color=color, lw=1.4)
    ax.plot(traj["mse_p"], label="MSE p", color="#5E81AC", lw=1.6, ls="--")
    ax.plot(traj["mmd_p"], label="MMD p", color="#B48EAD", lw=1.2, ls=":")
    ax.plot(traj["po_p"], label="PO p", color="#D08770", lw=1.2, ls="-.")
    ax.axvline(n_burn - 0.5, color="gray", ls="--", lw=0.9, label="end burn")
    fr = res["first_reject"]
    if fr["grad_any"] is not None:
        ax.axvline(fr["grad_any"], color="#BF616A", lw=1.0, alpha=0.7)
    if fr["mse"] is not None:
        ax.axvline(fr["mse"], color="#5E81AC", lw=1.0, alpha=0.7)
    ax.set_ylabel("p-value")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(ncol=4, fontsize=7, loc="upper right")
    ax.set_title(f"{name}: Grad-OnlineRFPerm vs MSE / MMD / PO")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    for n, color in zip(["fc1", "fc2", "fc3"], ["#BF616A", "#EBCB8B", "#A3BE8C"]):
        g = np.asarray(traj["grad_g"][n], float)
        ax.plot(g / (g.max() + 1e-12), label=f"‖g‖ {n} (rel)", color=color, lw=1.3)
    ax.plot(np.asarray(traj["train_mse"], float), label="train MSE", color="#4C566A", lw=1.2)
    ax.set_ylabel("grad / MSE")
    ax.legend(ncol=4, fontsize=7)
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    for key, color, lab in [
        ("mse_rej", "#5E81AC", "MSE reject"),
        ("mmd_rej", "#B48EAD", "MMD reject"),
        ("po_rej", "#D08770", "PO reject"),
    ]:
        ax.step(range(len(traj[key])), traj[key], where="post", label=lab, color=color, lw=1.2)
    # any-layer grad reject
    any_g = np.maximum.reduce([traj["grad_rej"][n] for n in ["fc1", "fc2", "fc3"]])
    ax.step(range(len(any_g)), any_g, where="post", label="Grad any reject", color="#BF616A", lw=1.5)
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
        "Continuous-time OnlineRFPerm on **per-layer gradient energy** (and relative",
        "shares), vs serving-MSE OnlineRFPerm / consecutive-batch RBF-MMD² / window PO.",
        "",
        "Lead time: `t_grad − t_mse` (negative ⇒ Grad rejects **before** MSE break).",
        "",
        f"- batch_size={meta['batch_size']}, n_batches={meta['n_batches']}, "
        f"n_burn={meta['n_burn']}, alpha={meta['alpha']}, seed={meta['seed']}",
        "",
        "| dataset | t_grad | t_mse | t_mmd | t_po | lead(g−mse) | lead(share−mse) | mse duty |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, r in all_res.items():
        fr = r["first_reject"]
        lt = r["lead_time"]

        def _f(x):
            return "—" if x is None else str(x)

        lines.append(
            f"| `{name}` | {_f(fr['grad_any'])} | {_f(fr['mse'])} | {_f(fr['mmd'])} | "
            f"{_f(fr['po'])} | {_f(lt['grad_minus_mse'])} | {_f(lt['share_minus_mse'])} | "
            f"{r['duty']['mse']:.2f} |"
        )
    lines += [
        "",
        "## Per-layer first Grad reject",
        "",
        "| dataset | fc1 | fc2 | fc3 |",
        "|---|---:|---:|---:|",
    ]
    for name, r in all_res.items():
        pl = r["first_reject"]["grad_per_layer"]

        def _f(x):
            return "—" if x is None else str(x)

        lines.append(f"| `{name}` | {_f(pl.get('fc1'))} | {_f(pl.get('fc2'))} | {_f(pl.get('fc3'))} |")
    lines += [
        "",
        "## Method",
        "",
        "```",
        "for each stream batch t:",
        "  train MLP a few Adam steps; capture ||grad||_2 per Linear layer",
        "  T_grad = ||g_t|| - e_ref   (burn-in sets e_ref)",
        "  p = rank/EWMA vs T history; alpha-investing FDR → reject",
        "  also: MSE-OnlineRFPerm(f_ref), MMD²(X_{t-1},X_t), PO window gap",
        "lead = t_first_grad_reject - t_first_mse_reject",
        "```",
        "",
        "Freeze hint: earliest rejecting layer (often early `fc1` under covariate",
        "shift, later `fc3` under concept shift) can guide back-prop depth.",
        "",
    ]
    path = out_dir / "AGOD_grad_rfperm_monitor.md"
    path.write_text("\n".join(lines) + "\n")
    # also mirror under docs/agod
    docs = ROOT / "docs" / "agod" / "AGOD_grad_rfperm_monitor.md"
    docs.parent.mkdir(parents=True, exist_ok=True)
    docs.write_text(path.read_text())
    return path


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
    ap.add_argument("--train-steps", type=int, default=12)
    ap.add_argument("--lr", type=float, default=1e-2)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "results" / "agod_grad_rfperm")
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    all_res: Dict[str, dict] = {}

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
        all_res[name] = res
        fr, lt = res["first_reject"], res["lead_time"]
        print(
            f"  first: grad={fr['grad_any']} mse={fr['mse']} mmd={fr['mmd']} po={fr['po']} "
            f"| lead(g-mse)={lt['grad_minus_mse']} share-mse={lt['share_minus_mse']}",
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


if __name__ == "__main__":
    main()
