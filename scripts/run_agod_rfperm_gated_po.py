#!/usr/bin/env python3
"""OnlineRFPerm-gated √PO IPTW on stream packs (PDF Algorithm 1).

Modes:
  uniform      w=1
  sqrt         always w∝√PO
  sqrt_gated   OnlineRFPerm reject → √PO else uniform
  dre          logistic density-ratio

Two held-out packs by default: stocks_MSFT, stocks_IWM.

  PYTHONPATH=. python3 scripts/run_agod_rfperm_gated_po.py \\
    --datasets stocks_MSFT stocks_IWM --batch-size 100 --n-batches 40
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error

from agod.online_rfperm import fit_online_rfperm, gate_blend, update_online_rfperm
from agod.po_iptw import dre_weights, instance_po_risk, po_iptw_weights
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

MODES = ("uniform", "sqrt", "sqrt_gated", "dre")
COLORS = {
    "uniform": "#4C566A",
    "sqrt": "#A3BE8C",
    "sqrt_gated": "#88C0D0",
    "dre": "#BF616A",
}


def ensure_loaders() -> None:
    if "stocks_MSFT" not in LOADERS:
        LOADERS["stocks_MSFT"] = lambda root, max_n=20000: load_stocks(root, "MSFT", max_n)
    if "stocks_IWM" not in LOADERS:
        LOADERS["stocks_IWM"] = lambda root, max_n=20000: load_stocks(root, "IWM", max_n)


def pca_fit(X, d, seed):
    if X.shape[1] <= d:
        return X.astype(np.float32)
    return PCA(n_components=d, random_state=seed).fit_transform(X).astype(np.float32)


def fit_rf(X, y, w, seed):
    m = RandomForestRegressor(
        n_estimators=40, max_depth=8, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    m.fit(X, y, sample_weight=w)
    return m


def batch_po(X0, y0, X1, y1, seed):
    if len(X0) < 16 or len(X1) < 16:
        return 0.0
    rf = RandomForestRegressor(
        n_estimators=20, max_depth=4, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    rf.fit(X0, y0)
    e0 = float(np.mean(np.abs(y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(y1 - rf.predict(X1))))
    return max(e1 - e0, 0.0)


def make_stream(X, y, bs, n_batches):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, len(X), bs)]


def run_mode(stream, mode: str, seed: int, *, n_burn: int, alpha: float, soft_gate: bool) -> dict:
    mse_next, mae_next = [], []
    gate_on, p_traj, T_traj = [], [], []
    X0, y0 = stream[0]
    probe = fit_rf(X0, y0, np.ones(len(y0)), seed)
    # OnlineRFPerm: f_ref on first batch (= production model)
    rfperm = fit_online_rfperm(X0, y0, seed=seed)

    for t in range(1, len(stream)):
        Xp, yp = stream[t - 1]
        Xc, yc = stream[t]
        burn = t <= n_burn
        step = update_online_rfperm(
            rfperm, Xc, yc, burn_in=burn, alpha=alpha, ewma=True, fdr="alpha_investing"
        )
        p_traj.append(float(step["p"]))
        T_traj.append(float(step["T"]))
        gate_on.append(int(step["reject"]))

        po_b = batch_po(Xp, yp, Xc, yc, seed + t)
        if mode == "dre":
            w = dre_weights(Xp, Xc, seed=seed + t)
        elif mode == "uniform":
            w = np.ones(len(yc))
        else:
            pred = probe.predict(Xc)
            po_row = instance_po_risk(yc, pred, batch_po=po_b, mix=0.5)
            w_sqrt = po_iptw_weights(po_row, mode="sqrt")
            if mode == "sqrt":
                w = w_sqrt
            else:  # sqrt_gated
                w = gate_blend(
                    bool(step["reject"]),
                    w_sqrt,
                    soft=soft_gate,
                    p=float(step["p"]),
                    alpha=alpha,
                )
                w = w / (w.mean() + 1e-8)

        model = fit_rf(Xc, yc, w, seed + 17 * t)
        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            pn = model.predict(Xn)
            mse_next.append(float(mean_squared_error(yn, pn)))
            mae_next.append(float(mean_absolute_error(yn, pn)))
        probe = fit_rf(Xc, yc, w, seed + 31 * t)

    mse = np.asarray(mse_next, float)
    mae = np.asarray(mae_next, float)
    duty = float(np.mean(gate_on)) if gate_on else 0.0
    return {
        "mode": mode,
        "mse_next": mse_next,
        "mae_next": mae_next,
        "mse_mean": float(mse.mean()) if len(mse) else float("nan"),
        "mse_std": float(mse.std()) if len(mse) else float("nan"),
        "mae_mean": float(mae.mean()) if len(mae) else float("nan"),
        "cum_mse": float(mse.sum()) if len(mse) else float("nan"),
        "gate_duty": duty if mode == "sqrt_gated" else (1.0 if mode == "sqrt" else 0.0),
        "p_traj": p_traj,
        "T_traj": T_traj,
        "gate_on": gate_on,
    }


def regrets(results: Dict[str, dict]) -> Dict[str, dict]:
    base_u = np.asarray(results["uniform"]["mse_next"], float)
    base_d = np.asarray(results["dre"]["mse_next"], float)
    out = {}
    for mode, r in results.items():
        m = np.asarray(r["mse_next"], float)
        ex_u = m - base_u[: len(m)]
        ex_d = m - base_d[: len(m)]
        out[mode] = {
            "cum_regret_vs_uniform": float(ex_u.sum()),
            "cum_regret_vs_dre": float(ex_d.sum()),
            "win_rate_vs_uniform": float(np.mean(ex_u < 0)),
            "win_rate_vs_dre": float(np.mean(ex_d < 0)),
            "cum_excess_traj_vs_dre": [float(v) for v in np.cumsum(ex_d)],
            "cum_excess_traj_vs_uniform": [float(v) for v in np.cumsum(ex_u)],
        }
    return out


def plot_all(all_res: dict, out_dir: Path) -> List[Path]:
    paths = []
    ds = list(all_res.keys())
    # relative MSE bars
    fig, ax = plt.subplots(figsize=(max(7, 1.8 * len(ds)), 4.5))
    x = np.arange(len(ds))
    w = 0.18
    for i, mode in enumerate(MODES):
        vals = []
        for d in ds:
            u = all_res[d]["results"]["uniform"]["mse_mean"]
            vals.append(all_res[d]["results"][mode]["mse_mean"] / (u + 1e-12))
        ax.bar(x + (i - 1.5) * w, vals, w, label=mode, color=COLORS[mode])
    ax.axhline(1.0, color="black", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("Next-MSE / uniform")
    ax.set_title("OnlineRFPerm-gated √PO vs always-√PO / DRE")
    ax.legend(ncol=4, fontsize=8)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "mse_rel_bars.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    # gate duty + p for gated
    fig, axes = plt.subplots(1, len(ds), figsize=(4.2 * len(ds), 3.6), squeeze=False)
    for j, d in enumerate(ds):
        ax = axes[0, j]
        g = all_res[d]["results"]["sqrt_gated"]
        ax.plot(g["p_traj"], label="p_t", color="#5E81AC", lw=1.5)
        ax.step(
            range(len(g["gate_on"])),
            g["gate_on"],
            where="post",
            label="gate ON",
            color="#A3BE8C",
            alpha=0.8,
        )
        ax.set_title(f"{d} (duty={g['gate_duty']:.2f})")
        ax.set_xlabel("batch t")
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    fig.suptitle("OnlineRFPerm p-values & gate", y=1.02)
    fig.tight_layout()
    p = out_dir / "gate_p_traj.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    paths.append(p)

    # cum regret gated vs dre
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for d in ds:
        ax.plot(
            all_res[d]["regret"]["sqrt_gated"]["cum_excess_traj_vs_dre"],
            label=f"{d} gated",
            lw=1.8,
        )
        ax.plot(
            all_res[d]["regret"]["sqrt"]["cum_excess_traj_vs_dre"],
            label=f"{d} always-√",
            lw=1.2,
            ls="--",
            alpha=0.7,
        )
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("stream step")
    ax.set_ylabel("cum excess MSE vs DRE (↓ better)")
    ax.set_title("Gated vs always-√PO cumulative regret vs DRE")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = out_dir / "cum_regret_gated_vs_dre.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)
    return paths


def report(all_res: dict) -> str:
    lines = [
        "# OnlineRFPerm-gated √PO — held-out stream packs",
        "",
        "Gate = PDF OnlineRFPerm (rank/EWMA p + alpha-investing FDR).",
        "Reweight only when reject: `w=√PO`; else uniform.",
        "",
        "Primary metric = next-MSE on **significant (reject) batches only**.",
        "",
        "| dataset | n_sig | unif | always-√ | **gated-√** | dre | duty | best |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for d, blob in all_res.items():
        res = blob["results"]
        means = {m: res[m].get("mse_mean_sig", res[m]["mse_mean"]) for m in MODES}
        best = min(means, key=means.get)
        n_sig = res["sqrt_gated"].get("n_significant", "?")
        lines.append(
            f"| `{d}` | {n_sig} | {means['uniform']:.4g} | {means['sqrt']:.4g} | "
            f"**{means['sqrt_gated']:.4g}** | {means['dre']:.4g} | "
            f"{res['sqrt_gated']['gate_duty']:.2f} | `{best}` |"
        )
    lines += [
        "",
        "```python",
        "# OnlineRFPerm gate (PDF Alg.1)",
        "T = MSE(f_ref, batch) - E_ref",
        "p = rank/EWMA vs historical T; online FDR → reject",
        "w = sqrt(PO) if reject else 1",
        "metric = mean(mse_next[reject])  # drop non-significant",
        "```",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ensure_loaders()
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--datasets", nargs="+", default=["stocks_MSFT", "stocks_IWM"])
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--n-burn", type=int, default=5)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--soft-gate", action="store_true")
    ap.add_argument("--max-n", type=int, default=8000)
    ap.add_argument("--pca-d", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_rfperm_gated"))
    args = ap.parse_args()

    need = args.batch_size * args.n_batches
    all_res: Dict[str, dict] = {}
    skipped: List[str] = []

    for name in args.datasets:
        loader = LOADERS.get(name)
        if loader is None:
            print(f"[skip] unknown {name}")
            skipped.append(name)
            continue
        packed = loader(args.root, max_n=max(need, args.max_n))
        if packed is None:
            print(f"[skip] {name}: missing")
            skipped.append(name)
            continue
        X, y, meta = packed
        if len(X) < need:
            print(f"[skip] {name}: need {need}, got {len(X)}")
            skipped.append(name)
            continue
        X = pca_fit(X[:need], args.pca_d, args.seed)
        y = y[:need]
        stream = make_stream(X, y, args.batch_size, args.n_batches)
        print(
            f"=== {name} n={len(X)} d={X.shape[1]} batches={len(stream)} burn={args.n_burn} ===",
            flush=True,
        )
        results = {}
        for mode in MODES:
            print(f"  [{mode}]", flush=True)
            results[mode] = run_mode(
                stream,
                mode,
                args.seed,
                n_burn=args.n_burn,
                alpha=args.alpha,
                soft_gate=args.soft_gate,
            )
            print(
                f"    mse={results[mode]['mse_mean']:.6g} "
                f"duty={results[mode]['gate_duty']:.2f}"
            )
        annotate_results_with_sig(results, preferred_gate_modes=("sqrt_gated",))
        all_res[name] = {"results": results, "regret": regrets(results), "meta": meta}

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "n_burn": args.n_burn,
        "alpha": args.alpha,
        "soft_gate": args.soft_gate,
        "modes": list(MODES),
        "skipped": skipped,
        "datasets": all_res,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_res)
    (args.out / "RFPERM_GATED_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_rfperm_gated_po.md").write_text(md, encoding="utf-8")
    plots = plot_all(all_res, args.out) if all_res else []
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
