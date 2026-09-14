#!/usr/bin/env python3
"""OnlineRFPerm gate → post-hoc T=0/T=1 PO re-fit → √PO weights.

Elaboration
-----------
1. Uniform is the default (should look slightly best overall).
2. OnlineRFPerm: small p / FDR reject ⇒ batch is *significantly* different.
3. Only then post-hoc reweight:
     T=0 : recent control window (batch just before the pair)
     T=1 : previous batch ∪ current batch
   Re-fit μ0 (and optional μ1); PO_i = |Y−μ0(X)| (+ μ-gap blend) on T=1;
   take √PO weights on the *current* slice; fit the downstream RF.

Modes: uniform | sqrt (always, stale probe residual) | sqrt_gated
       | sqrt_gated_refit (RFPerm + T=0/T=1 re-fit) | dre

  PYTHONPATH=. python3 scripts/run_agod_rfperm_po_refit.py \\
    --datasets metro_interstate beijing_pm25 stocks_AAPL waymo_proxy
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
from agod.po_refit import current_batch_sqrt_po_weights
from agod.stream_packs import LOADERS, load_stocks

MODES = ("uniform", "sqrt", "sqrt_gated", "sqrt_gated_refit", "dre")
COLORS = {
    "uniform": "#4C566A",
    "sqrt": "#A3BE8C",
    "sqrt_gated": "#88C0D0",
    "sqrt_gated_refit": "#5E81AC",
    "dre": "#BF616A",
}


def ensure_loaders() -> None:
    for t in ("MSFT", "IWM", "AAPL"):
        key = f"stocks_{t}"
        if key not in LOADERS:
            LOADERS[key] = (lambda ticker: (lambda root, max_n=20000: load_stocks(root, ticker, max_n)))(t)


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


def run_mode(
    stream,
    mode: str,
    seed: int,
    *,
    n_burn: int,
    alpha: float,
    soft_gate: bool,
    n_control: int,
) -> dict:
    mse_next, mae_next = [], []
    gate_on, p_traj, T_traj, refit_on = [], [], [], []
    X0, y0 = stream[0]
    probe = fit_rf(X0, y0, np.ones(len(y0)), seed)
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
        rejected = bool(step["reject"])
        gate_on.append(int(rejected))
        did_refit = 0

        if mode == "dre":
            w = dre_weights(Xp, Xc, seed=seed + t)
        elif mode == "uniform":
            w = np.ones(len(yc))
        elif mode == "sqrt":
            po_b = batch_po(Xp, yp, Xc, yc, seed + t)
            pred = probe.predict(Xc)
            po_row = instance_po_risk(yc, pred, batch_po=po_b, mix=0.5)
            w = po_iptw_weights(po_row, mode="sqrt")
        elif mode == "sqrt_gated":
            po_b = batch_po(Xp, yp, Xc, yc, seed + t)
            pred = probe.predict(Xc)
            po_row = instance_po_risk(yc, pred, batch_po=po_b, mix=0.5)
            w_sqrt = po_iptw_weights(po_row, mode="sqrt")
            w = gate_blend(rejected, w_sqrt, soft=soft_gate, p=float(step["p"]), alpha=alpha)
            w = w / (w.mean() + 1e-8)
        else:  # sqrt_gated_refit
            if rejected and t >= 1:
                w = current_batch_sqrt_po_weights(
                    stream, t, seed=seed + 7 * t, n_control=n_control
                )
                did_refit = 1
            else:
                w = np.ones(len(yc))
        refit_on.append(did_refit)

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
    refit_duty = float(np.mean(refit_on)) if refit_on else 0.0
    return {
        "mode": mode,
        "mse_next": mse_next,
        "mae_next": mae_next,
        "mse_mean": float(mse.mean()) if len(mse) else float("nan"),
        "mse_std": float(mse.std()) if len(mse) else float("nan"),
        "mae_mean": float(mae.mean()) if len(mae) else float("nan"),
        "cum_mse": float(mse.sum()) if len(mse) else float("nan"),
        "gate_duty": duty,
        "refit_duty": refit_duty if mode == "sqrt_gated_refit" else duty,
        "p_traj": p_traj,
        "T_traj": T_traj,
        "gate_on": gate_on,
        "refit_on": refit_on,
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
    fig, ax = plt.subplots(figsize=(max(8, 1.9 * len(ds)), 4.6))
    x = np.arange(len(ds))
    w = 0.15
    for i, mode in enumerate(MODES):
        vals = [
            all_res[d]["results"][mode]["mse_mean"]
            / (all_res[d]["results"]["uniform"]["mse_mean"] + 1e-12)
            for d in ds
        ]
        ax.bar(x + (i - 2) * w, vals, w, label=mode, color=COLORS[mode])
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("Next-MSE / uniform")
    ax.set_title("RFPerm + T0/T1 PO re-fit vs always-√ / DRE")
    ax.legend(ncol=3, fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "mse_rel_bars.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    for d in ds:
        ax.plot(
            all_res[d]["regret"]["sqrt_gated_refit"]["cum_excess_traj_vs_dre"],
            label=f"{d} refit",
            lw=1.8,
        )
        ax.plot(
            all_res[d]["regret"]["sqrt_gated"]["cum_excess_traj_vs_dre"],
            label=f"{d} gated",
            lw=1.1,
            ls="--",
            alpha=0.75,
        )
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("stream step")
    ax.set_ylabel("cum excess MSE vs DRE (↓ better)")
    ax.set_title("T0/T1 PO re-fit vs probe-gated √PO (vs DRE)")
    ax.legend(fontsize=6, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p = out_dir / "cum_regret_refit_vs_dre.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    fig, axes = plt.subplots(1, len(ds), figsize=(4.0 * len(ds), 3.5), squeeze=False)
    for j, d in enumerate(ds):
        ax = axes[0, j]
        g = all_res[d]["results"]["sqrt_gated_refit"]
        ax.plot(g["p_traj"], color="#5E81AC", lw=1.3, label="p_t")
        ax.step(range(len(g["refit_on"])), g["refit_on"], where="post", color="#A3BE8C", label="refit ON")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"{d} duty={g['refit_duty']:.2f}")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
    fig.suptitle("OnlineRFPerm p + T0/T1 PO re-fit duty", y=1.02)
    fig.tight_layout()
    p = out_dir / "refit_duty_traj.png"
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    paths.append(p)
    return paths


def report(all_res: dict) -> str:
    lines = [
        "# OnlineRFPerm → T=0/T=1 PO re-fit → √PO (post-hoc)",
        "",
        "Default = **uniform**. On significant OnlineRFPerm reject only:",
        "",
        "- `T=0`: recent control batch(es) before the pair",
        "- `T=1`: previous ∪ current batch",
        "- re-fit μ0 (optional μ1); `PO=|Y−μ0(X)|`; `w=√PO` on current batch",
        "",
        "| dataset | unif | always-√ | gated | **refit** | dre | refit-duty | best |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    beats_dre = beats_always = beats_gated = 0
    n = 0
    for d, blob in all_res.items():
        res = blob["results"]
        means = {m: res[m]["mse_mean"] for m in MODES}
        best = min(means, key=means.get)
        n += 1
        if means["sqrt_gated_refit"] < means["dre"]:
            beats_dre += 1
        if means["sqrt_gated_refit"] < means["sqrt"]:
            beats_always += 1
        if means["sqrt_gated_refit"] < means["sqrt_gated"]:
            beats_gated += 1
        lines.append(
            f"| `{d}` | {means['uniform']:.4g} | {means['sqrt']:.4g} | "
            f"{means['sqrt_gated']:.4g} | **{means['sqrt_gated_refit']:.4g}** | "
            f"{means['dre']:.4g} | {res['sqrt_gated_refit']['refit_duty']:.2f} | `{best}` |"
        )
    lines += [
        "",
        f"**refit vs DRE:** `{beats_dre}/{n}`",
        f"**refit vs always-√:** `{beats_always}/{n}`",
        f"**refit vs probe-gated:** `{beats_gated}/{n}`",
        "",
        "Expectation: uniform slightly best overall; refit fires sparsely on",
        "clearly-shifted batches and should beat always-√ / DRE there.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ensure_loaders()
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=[
            "metro_interstate",
            "beijing_pm25",
            "stocks_AAPL",
            "waymo_proxy",
            "stocks_MSFT",
            "stocks_IWM",
        ],
    )
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--n-burn", type=int, default=5)
    ap.add_argument("--n-control", type=int, default=1)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--soft-gate", action="store_true")
    ap.add_argument("--max-n", type=int, default=8000)
    ap.add_argument("--pca-d", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_rfperm_po_refit"))
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
        if packed is None or len(packed[0]) < need:
            print(f"[skip] {name}")
            skipped.append(name)
            continue
        X, y, meta = packed
        X = pca_fit(X[:need], args.pca_d, args.seed)
        y = y[:need]
        stream = make_stream(X, y, args.batch_size, args.n_batches)
        print(
            f"=== {name} n={len(X)} d={X.shape[1]} batches={len(stream)} "
            f"burn={args.n_burn} control={args.n_control} ===",
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
                n_control=args.n_control,
            )
            print(
                f"    mse={results[mode]['mse_mean']:.6g} "
                f"gate={results[mode]['gate_duty']:.2f} "
                f"refit={results[mode]['refit_duty']:.2f}"
            )
        all_res[name] = {"results": results, "regret": regrets(results), "meta": meta}

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "n_burn": args.n_burn,
        "n_control": args.n_control,
        "alpha": args.alpha,
        "modes": list(MODES),
        "skipped": skipped,
        "datasets": all_res,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_res)
    (args.out / "RFPERM_PO_REFIT_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_rfperm_po_refit.md").write_text(md, encoding="utf-8")
    plots = plot_all(all_res, args.out) if all_res else []
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
