#!/usr/bin/env python3
"""Continuous-batch PO-risk IPTW: uniform / prop / sqrt / inv → MSE @ 10 batches.

Green (no DRE/DGA):
  uniform : w = 1
  prop    : w ∝ PO
  sqrt    : w ∝ √PO     # soft upweight — w_i = sqrt(PO(X_i,Y_i,T_i=1))
  inv     : w ∝ 1/PO

Eval: fit batch t with weights → MSE on batch t+1.

  PYTHONPATH=. python3 scripts/run_agod_po_iptw_mse.py \\
    --dataset affec --n-batches 10 --batch-size 256
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

from agod.po_iptw import instance_po_risk, po_iptw_weights


def batch_po_risk(X0, y0, X1, y1, *, seed: int) -> float:
    """Window PO: ref RF absolute-error gap on cur (no DRE)."""
    if len(X0) < 16 or len(X1) < 16:
        return 0.0
    rf = RandomForestRegressor(
        n_estimators=20, max_depth=4, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    rf.fit(X0, y0)
    e0 = float(np.mean(np.abs(y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(y1 - rf.predict(X1))))
    return max(e1 - e0, 0.0)


def load_affec(root: Path, max_n: int, seed: int) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    cache = root / "results/affec_fsds/affec_fsds_xyw_cache.npz"
    if not cache.is_file():
        return None
    z = np.load(cache, allow_pickle=True)
    X, Y = z["X"].astype(np.float32), z["Y"].astype(np.float64)
    n = min(max_n, len(X))
    idx = np.random.default_rng(seed).choice(len(X), size=n, replace=False)
    return X[idx], Y[idx]


def load_synth(n: int, d: int, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, d)).astype(np.float32)
    y = np.zeros(n, float)
    chunk = max(n // 10, 1)
    for b in range(10):
        lo, hi = b * chunk, n if b == 9 else (b + 1) * chunk
        w = rng.normal(size=d)
        w /= np.linalg.norm(w) + 1e-9
        noise = 0.3 + 0.4 * (b % 3)
        y[lo:hi] = X[lo:hi] @ w * (1.0 + 0.15 * b) + rng.normal(0, noise, hi - lo)
    return X, y


def make_stream(X: np.ndarray, y: np.ndarray, batch: int) -> List[Tuple[np.ndarray, np.ndarray]]:
    return [(X[i : i + batch], y[i : i + batch]) for i in range(0, len(X) - batch + 1, batch)]


class MLP(nn.Module):
    def __init__(self, d: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d, 64), nn.ReLU(), nn.Linear(64, 32), nn.ReLU(), nn.Linear(32, 1)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


def fit_rf(X, y, w, seed) -> RandomForestRegressor:
    rf = RandomForestRegressor(
        n_estimators=40, max_depth=6, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    rf.fit(X, y, sample_weight=w)
    return rf


def fit_mlp(X, y, w, seed, steps: int = 80) -> MLP:
    torch.manual_seed(seed)
    m = MLP(X.shape[1])
    opt = torch.optim.Adam(m.parameters(), lr=1e-2)
    xt = torch.from_numpy(np.asarray(X, np.float32))
    yt = torch.from_numpy(np.asarray(y, np.float32))
    wt = torch.from_numpy(np.asarray(w, np.float32))
    m.train()
    for _ in range(steps):
        loss = (wt * (m(xt) - yt) ** 2).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
    m.eval()
    return m


@torch.no_grad()
def pred_mlp(m: MLP, X) -> np.ndarray:
    return m(torch.from_numpy(np.asarray(X, np.float32))).numpy()


def _pack(xs: List[float]) -> dict:
    a = np.asarray(xs, float)
    return {
        "mean": float(a.mean()),
        "std": float(a.std()),
        "p50": float(np.median(a)),
        "p90": float(np.percentile(a, 90)),
        "traj": [float(v) for v in a],
    }


def run_mode(stream, mode: str, model_kind: str, seed: int) -> dict:
    mse_cur: List[float] = []
    mse_next: List[float] = []
    batch_po: List[float] = []
    X0, y0 = stream[0]
    probe = fit_rf(X0, y0, np.ones(len(y0)), seed)

    for t in range(1, len(stream)):
        Xp, yp = stream[t - 1]
        Xc, yc = stream[t]
        po_b = batch_po_risk(Xp, yp, Xc, yc, seed=seed + t)
        batch_po.append(po_b)

        po_row = instance_po_risk(yc, probe.predict(Xc), batch_po=po_b, mix=0.5)
        w = po_iptw_weights(po_row, mode=mode)  # type: ignore[arg-type]

        if model_kind == "rf":
            model = fit_rf(Xc, yc, w, seed + 17 * t)
            pred_c = model.predict(Xc)
            mse_cur.append(float(mean_squared_error(yc, pred_c)))
            if t + 1 < len(stream):
                Xn, yn = stream[t + 1]
                mse_next.append(float(mean_squared_error(yn, model.predict(Xn))))
        else:
            model = fit_mlp(Xc, yc, w, seed + 17 * t)
            pred_c = pred_mlp(model, Xc)
            mse_cur.append(float(mean_squared_error(yc, pred_c)))
            if t + 1 < len(stream):
                Xn, yn = stream[t + 1]
                mse_next.append(float(mean_squared_error(yn, pred_mlp(model, Xn))))

        probe = fit_rf(Xc, yc, w, seed + 31 * t)

    return {
        "mode": mode,
        "batch_po": batch_po,
        "mse_cur": _pack(mse_cur),
        "mse_next": _pack(mse_next),
    }


def report_md(results: Dict[str, dict]) -> str:
    lines = [
        "# PO-risk IPTW continuous-batch MSE",
        "",
        "| mode | MSE_next mean | std | p50 | p90 | MSE_cur mean |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for m, r in results.items():
        n, c = r["mse_next"], r["mse_cur"]
        lines.append(
            f"| `{m}` | {n['mean']:.4f} | {n['std']:.4f} | {n['p50']:.4f} | "
            f"{n['p90']:.4f} | {c['mean']:.4f} |"
        )
    best = min(results, key=lambda k: results[k]["mse_next"]["mean"])
    lines += [
        "",
        f"**Best next-batch MSE mean:** `{best}`",
        "",
        "```python",
        "w = po / po.mean()              # prop",
        "w = np.sqrt(po) / mean          # sqrt  ← soft high-PO upweight",
        "w = (1/po) / (1/po).mean()      # inv",
        "rf.fit(X, y, sample_weight=w)",
        "```",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--dataset", choices=["affec", "synth"], default="affec")
    ap.add_argument("--n-batches", type=int, default=10)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--max-n", type=int, default=3000)
    ap.add_argument("--model", choices=["rf", "mlp"], default="rf")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_po_iptw"))
    args = ap.parse_args()

    if args.dataset == "affec":
        packed = load_affec(args.root, args.max_n, args.seed)
        if packed is None:
            print("[warn] no affec cache → synth")
            X, y = load_synth(args.max_n, 32, args.seed)
            ds = "synth_fallback"
        else:
            X, y = packed
            ds = "affec"
    else:
        X, y = load_synth(args.max_n, 32, args.seed)
        ds = "synth"

    need = args.n_batches * args.batch_size
    if len(X) < need:
        raise SystemExit(f"need>={need} rows, got {len(X)}")
    stream = make_stream(X[:need], y[:need], args.batch_size)
    print(f"dataset={ds} batches={len(stream)} bs={args.batch_size} model={args.model}")

    results: Dict[str, dict] = {}
    for mode in ("uniform", "prop", "sqrt", "inv"):
        print(f"  [{mode}] ...", flush=True)
        results[mode] = run_mode(stream, mode, args.model, args.seed)
        n = results[mode]["mse_next"]
        print(f"    next MSE mean={n['mean']:.4f} std={n['std']:.4f} p90={n['p90']:.4f}")

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "dataset": ds,
        "model": args.model,
        "n_batches": len(stream),
        "results": results,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report_md(results)
    (args.out / "PO_IPTW_MSE_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_po_iptw_mse.md").write_text(md, encoding="utf-8")
    print(md)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
