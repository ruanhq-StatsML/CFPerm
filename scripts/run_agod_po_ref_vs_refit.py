#!/usr/bin/env python3
"""Compare reference PO-risk vs OnlineRFPerm-triggered PO-learner re-fit.

Eval logic (deliberately blunt)
------------------------------
1. OnlineRFPerm opens only on significant batches.
2. On those steps, three PO scores on the *current* batch:
     ref_po   = |Y − f_ref(X)|     frozen burn-in reference
     probe_po = |Y − probe(X)|     rolling prev-batch probe
     refit_po = |Y − μ0_R(X)|      μ0 re-trained on recent control R
3. Quality: Spearman / top-20% overlap of PO vs oracle |Y − μ_oracle|,
   where μ_oracle is fit on recent∪current (uniform) — a proxy for
   "which rows are actually hard".
4. Downstream: IPTW w∝√PO only on reject; else w=1. Report next-MSE on
   **significant batches only** (non-reject ≡ uniform → drop them).

  PYTHONPATH=. python3 scripts/run_agod_po_ref_vs_refit.py \\
    --datasets metro_interstate beijing_pm25 stocks_AAPL waymo_proxy stocks_MSFT stocks_IWM
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

from agod.hard_rank_metrics import aggregate_hard_rank
from agod.online_rfperm import fit_online_rfperm, update_online_rfperm
from agod.po_iptw import dre_weights
from agod.po_ref_vs_refit import (
    po_quality_vs_truth,
    probe_po,
    reference_po_from_fref,
    refit_po_current,
    weights_from_po,
)
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

MODES = ("uniform", "ref_po", "probe_po", "refit_po", "dre")
COLORS = {
    "uniform": "#4C566A",
    "ref_po": "#EBCB8B",
    "probe_po": "#88C0D0",
    "refit_po": "#5E81AC",
    "dre": "#BF616A",
}


def ensure_loaders() -> None:
    for t in ("MSFT", "IWM", "AAPL"):
        key = f"stocks_{t}"
        if key not in LOADERS:
            LOADERS[key] = (
                lambda ticker: (lambda root, max_n=20000: load_stocks(root, ticker, max_n))
            )(t)


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


def make_stream(X, y, bs, n_batches):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, len(X), bs)]


def oracle_hardness(stream, t: int, seed: int) -> np.ndarray:
    """Proxy truth: |y − μ| on current; μ fit on recent∪current (uniform)."""
    Xp, yp = stream[t - 1]
    Xc, yc = stream[t]
    X = np.vstack([Xp, Xc])
    y = np.concatenate([np.asarray(yp, float), np.asarray(yc, float)])
    mu = fit_rf(X, y, np.ones(len(y)), seed)
    return np.abs(np.asarray(yc, float) - mu.predict(Xc))


def _agg_q(qs: List[dict], key: str) -> float:
    vals = [q[key] for q in qs if key in q and q[key] == q[key]]
    return float(np.mean(vals)) if vals else float("nan")


def _pack_q(qs: List[dict], po_means: List[float]) -> dict:
    """Aggregate full hard-rank metric suite over reject batches."""
    agg = aggregate_hard_rank(qs) if qs else {}
    return {
        "spearman": agg.get("spearman", float("nan")),
        "pearson": agg.get("pearson", float("nan")),
        "precision_at_k": agg.get("precision_at_k", agg.get("topk_overlap", float("nan"))),
        "topk_overlap": agg.get("topk_overlap", agg.get("precision_at_k", float("nan"))),
        "lift_at_k": agg.get("lift_at_k", float("nan")),
        "ndcg_at_k": agg.get("ndcg_at_k", float("nan")),
        "auroc_topk": agg.get("auroc_topk", float("nan")),
        "precision_at_10pct": agg.get("precision_at_10pct", float("nan")),
        "precision_at_20pct": agg.get("precision_at_20pct", float("nan")),
        "precision_at_30pct": agg.get("precision_at_30pct", float("nan")),
        "po_mean": float(np.mean(po_means)) if po_means else float("nan"),
        "n_batches": float(len(qs)),
    }


def run_mode(
    stream,
    mode: str,
    seed: int,
    *,
    n_burn: int,
    alpha: float,
    n_control: int,
    window_mode: str,
) -> dict:
    mse_next, mae_next = [], []
    gate_on, p_traj, T_traj = [], [], []
    q_ref, q_probe, q_refit = [], [], []
    po_mean_ref, po_mean_probe, po_mean_refit = [], [], []

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
        rejected = bool(step["reject"])
        gate_on.append(int(rejected))
        p_traj.append(float(step["p"]))
        T_traj.append(float(step["T"]))

        po_ref = reference_po_from_fref(rfperm.f_ref, Xc, yc)
        po_prb = probe_po(probe, Xc, yc, batch_po=None, mix=0.0)
        po_rft = refit_po_current(
            stream,
            t,
            seed=seed + 7 * t,
            n_recent=n_control,
            window_mode=window_mode,
        )

        if rejected:
            truth = oracle_hardness(stream, t, seed + 99 * t)
            q_ref.append(po_quality_vs_truth(po_ref, truth))
            q_probe.append(po_quality_vs_truth(po_prb, truth))
            q_refit.append(po_quality_vs_truth(po_rft, truth))
            po_mean_ref.append(float(np.mean(po_ref)))
            po_mean_probe.append(float(np.mean(po_prb)))
            po_mean_refit.append(float(np.mean(po_rft)))

        if mode == "dre":
            w = dre_weights(Xp, Xc, seed=seed + t)
        elif mode == "uniform":
            w = np.ones(len(yc))
        elif mode == "ref_po":
            w = weights_from_po(po_ref, "sqrt") if rejected else np.ones(len(yc))
        elif mode == "probe_po":
            w = weights_from_po(po_prb, "sqrt") if rejected else np.ones(len(yc))
        elif mode == "refit_po":
            w = weights_from_po(po_rft, "sqrt") if rejected else np.ones(len(yc))
        else:
            raise ValueError(mode)

        model = fit_rf(Xc, yc, w, seed + 17 * t)
        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            pn = model.predict(Xn)
            mse_next.append(float(mean_squared_error(yn, pn)))
            mae_next.append(float(mean_absolute_error(yn, pn)))
        probe = fit_rf(Xc, yc, np.ones(len(yc)), seed + 31 * t)

    mse = np.asarray(mse_next, float)
    mae = np.asarray(mae_next, float)
    return {
        "mode": mode,
        "mse_next": mse_next,
        "mae_next": mae_next,
        "mse_mean": float(mse.mean()) if len(mse) else float("nan"),
        "mse_std": float(mse.std()) if len(mse) else float("nan"),
        "mae_mean": float(mae.mean()) if len(mae) else float("nan"),
        "gate_duty": float(np.mean(gate_on)) if gate_on else 0.0,
        "gate_on": gate_on,
        "p_traj": p_traj,
        "T_traj": T_traj,
        "po_quality": {
            "n_reject": len(q_ref),
            "ref": _pack_q(q_ref, po_mean_ref),
            "probe": _pack_q(q_probe, po_mean_probe),
            "refit": _pack_q(q_refit, po_mean_refit),
        },
    }


def plot_all(all_res: dict, out_dir: Path) -> List[Path]:
    paths = []
    ds = list(all_res.keys())
    x = np.arange(len(ds))

    fig, ax = plt.subplots(figsize=(max(8, 1.7 * len(ds)), 4.6))
    w = 0.15
    for i, mode in enumerate(MODES):
        vals = [
            all_res[d]["results"][mode].get(
                "rel_mse_vs_uniform_sig",
                all_res[d]["results"][mode]["mse_mean"]
                / (all_res[d]["results"]["uniform"]["mse_mean"] + 1e-12),
            )
            for d in ds
        ]
        ax.bar(x + (i - 2) * w, vals, w, label=mode, color=COLORS[mode])
    ax.axhline(1.0, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("Next-MSE / uniform (sig-only)")
    ax.set_title("Reference / probe / re-fit PO — reject steps only")
    ax.legend(ncol=3, fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "mse_rel_sig_ref_vs_refit.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(max(8, 1.5 * len(ds)), 4.4))
    w = 0.25
    for i, src in enumerate(("ref", "probe", "refit")):
        vals = [all_res[d]["results"]["refit_po"]["po_quality"][src]["spearman"] for d in ds]
        ax.bar(x + (i - 1) * w, vals, w, label=src)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("Spearman(PO, oracle |resid|)")
    ax.set_title("PO ranking quality on reject batches (↑ better)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "po_spearman_ref_vs_refit.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(max(8, 1.5 * len(ds)), 4.4))
    for i, src in enumerate(("ref", "probe", "refit")):
        vals = [
            all_res[d]["results"]["refit_po"]["po_quality"][src].get(
                "precision_at_k",
                all_res[d]["results"]["refit_po"]["po_quality"][src].get("topk_overlap", np.nan),
            )
            for d in ds
        ]
        ax.bar(x + (i - 1) * w, vals, w, label=src)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("Precision@20% (hard-row recovery)")
    ax.set_title("Top-hard Precision@k on reject batches (↑ better)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "po_precision_at_k_ref_vs_refit.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    fig, ax = plt.subplots(figsize=(max(8, 1.5 * len(ds)), 4.4))
    for i, src in enumerate(("ref", "probe", "refit")):
        vals = [all_res[d]["results"]["refit_po"]["po_quality"][src].get("auroc_topk", np.nan) for d in ds]
        ax.bar(x + (i - 1) * w, vals, w, label=src)
    ax.axhline(0.5, color="k", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(ds, rotation=15, ha="right")
    ax.set_ylabel("AUROC (truth top-20% = positive)")
    ax.set_title("Hard-class AUROC on reject batches (↑ better)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    p = out_dir / "po_auroc_ref_vs_refit.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)
    return paths


def report(all_res: dict) -> str:
    lines = [
        "# Hard-sample ranking: reference PO vs re-fit PO-learner",
        "",
        "Primary claim = **PO ranks hard OOD rows**. Downstream MSE is secondary.",
        "",
        "On each OnlineRFPerm reject batch:",
        "",
        "1. `truth_i = |Y_i − μ_oracle(X_i)|` with μ_oracle fit on recent∪current (diagnostic).",
        "2. Score `ref_po` / `probe_po` / `refit_po` on the same rows.",
        "3. Measure ranking quality (Spearman, Precision@20%, Lift, NDCG, AUROC).",
        "",
        "## Hard-row ranking (reject batches only) — primary",
        "",
        "| dataset | spearman ref→probe→**refit** | P@20% ref→probe→**refit** | AUROC ref→probe→**refit** | Lift@20% **refit** | NDCG **refit** |",
        "|---|---|---|---|---:|---:|",
    ]
    for d, blob in all_res.items():
        q = blob["results"]["refit_po"]["po_quality"]
        lines.append(
            f"| `{d}` | "
            f"{q['ref']['spearman']:.2f}→{q['probe']['spearman']:.2f}→**{q['refit']['spearman']:.2f}** | "
            f"{q['ref']['precision_at_k']:.2f}→{q['probe']['precision_at_k']:.2f}→**{q['refit']['precision_at_k']:.2f}** | "
            f"{q['ref']['auroc_topk']:.2f}→{q['probe']['auroc_topk']:.2f}→**{q['refit']['auroc_topk']:.2f}** | "
            f"{q['refit']['lift_at_k']:.2f} | {q['refit']['ndcg_at_k']:.2f} |"
        )

    # ranking wins: who has best spearman / precision
    sp_wins = {"ref": 0, "probe": 0, "refit": 0}
    pk_wins = {"ref": 0, "probe": 0, "refit": 0}
    for d, blob in all_res.items():
        q = blob["results"]["refit_po"]["po_quality"]
        sp_wins[max(("ref", "probe", "refit"), key=lambda s: q[s]["spearman"])] += 1
        pk_wins[max(("ref", "probe", "refit"), key=lambda s: q[s]["precision_at_k"])] += 1

    lines += [
        "",
        f"**Best Spearman wins:** ref={sp_wins['ref']}, probe={sp_wins['probe']}, **refit={sp_wins['refit']}**",
        f"**Best Precision@20% wins:** ref={pk_wins['ref']}, probe={pk_wins['probe']}, **refit={pk_wins['refit']}**",
        "",
        "### Metric definitions (top 20% = hard)",
        "",
        "- **Spearman**: `corr(rank(PO), rank(truth))` — full order concordance.",
        "- **Precision@k**: `|Top_k(PO) ∩ Top_k(truth)| / k` — hard-set recovery.",
        "- **Lift@k**: Precision@k / (k/n) — vs random (1.0 = chance).",
        "- **NDCG@k**: graded by truth hardness — rewards ordering the *hardest* first.",
        "- **AUROC**: truth top-20% as positive class — threshold-free hard detection.",
        "",
        "## Downstream next-MSE (sig-only) — secondary",
        "",
        "| dataset | n_sig | unif | ref_po | probe_po | **refit_po** | dre | best |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    wins = {m: 0 for m in MODES}
    refit_beats_ref = refit_beats_probe = 0
    n = 0
    for d, blob in all_res.items():
        res = blob["results"]
        means = {m: res[m].get("mse_mean_sig", res[m]["mse_mean"]) for m in MODES}
        best = min(means, key=means.get)
        wins[best] += 1
        n += 1
        if means["refit_po"] < means["ref_po"]:
            refit_beats_ref += 1
        if means["refit_po"] < means["probe_po"]:
            refit_beats_probe += 1
        n_sig = res["refit_po"].get("n_significant", "?")
        lines.append(
            f"| `{d}` | {n_sig} | {means['uniform']:.4g} | {means['ref_po']:.4g} | "
            f"{means['probe_po']:.4g} | **{means['refit_po']:.4g}** | "
            f"{means['dre']:.4g} | `{best}` |"
        )
    lines += [
        "",
        f"**MSE wins (sig-only):** " + ", ".join(f"`{m}`={wins[m]}" for m in MODES),
        f"**refit_po < ref_po (MSE):** `{refit_beats_ref}/{n}`",
        f"**refit_po < probe_po (MSE):** `{refit_beats_probe}/{n}`",
        "",
        "See `docs/agod/AGOD_hard_rank_eval.md` for the full protocol.",
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
    ap.add_argument(
        "--window-mode",
        choices=("recent_ood", "prev_cur"),
        default="recent_ood",
    )
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--max-n", type=int, default=8000)
    ap.add_argument("--pca-d", type=int, default=32)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_po_ref_vs_refit"))
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
            f"window={args.window_mode} recent={args.n_control} ===",
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
                n_control=args.n_control,
                window_mode=args.window_mode,
            )
            print(
                f"    mse={results[mode]['mse_mean']:.6g} "
                f"gate={results[mode]['gate_duty']:.2f}"
            )
        annotate_results_with_sig(results, preferred_gate_modes=("refit_po", "ref_po"))
        q = results["refit_po"]["po_quality"]
        print(
            f"  [PO quality @reject] spearman "
            f"ref={q['ref']['spearman']:.3f} probe={q['probe']['spearman']:.3f} "
            f"refit={q['refit']['spearman']:.3f} | "
            f"topk ref={q['ref']['topk_overlap']:.3f} "
            f"probe={q['probe']['topk_overlap']:.3f} "
            f"refit={q['refit']['topk_overlap']:.3f}"
        )
        print(
            f"  [sig-only MSE] unif={results['uniform']['mse_mean_sig']:.6g} "
            f"ref={results['ref_po']['mse_mean_sig']:.6g} "
            f"probe={results['probe_po']['mse_mean_sig']:.6g} "
            f"refit={results['refit_po']['mse_mean_sig']:.6g}"
        )
        all_res[name] = {"results": results, "meta": meta}

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "n_burn": args.n_burn,
        "n_control": args.n_control,
        "window_mode": args.window_mode,
        "modes": list(MODES),
        "skipped": skipped,
        "datasets": {
            k: {
                "meta": v["meta"],
                "results": {
                    m: {
                        kk: vv
                        for kk, vv in r.items()
                        if kk
                        not in (
                            "mse_next",
                            "mae_next",
                            "p_traj",
                            "T_traj",
                            "gate_on",
                        )
                    }
                    for m, r in v["results"].items()
                },
            }
            for k, v in all_res.items()
        },
        "note": (
            "Frozen f_ref PO vs rolling probe vs recent-control PO re-fit after "
            "OnlineRFPerm reject; sig-only next-MSE + PO ranking quality."
        ),
    }
    # keep full trajectories in a separate heavy dump
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (args.out / "summary_full.json").write_text(
        json.dumps(
            {
                "datasets": all_res,
                "skipped": skipped,
                "window_mode": args.window_mode,
                "n_control": args.n_control,
            },
            indent=2,
            default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o,
        ),
        encoding="utf-8",
    )
    md = report(all_res)
    (args.out / "PO_REF_VS_REFIT_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_po_ref_vs_refit.md").write_text(md, encoding="utf-8")
    plots = plot_all(all_res, args.out) if all_res else []
    print(md)
    print("plots:", [str(p) for p in plots])


if __name__ == "__main__":
    main()
