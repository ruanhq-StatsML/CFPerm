#!/usr/bin/env python3
"""Obs-PO hard-reweight v4: hard-support + beijing-class packMSE gate.

Logic we buy
------------
1. **Hard-rank always** — obs PO ranks hard rows (Spearman / P@20%).
   Mechanism that matches: ``hard_support`` (top-k boost scaled by λ).
2. **Pack MSE only under beijing-class drift** — chase all-row next-MSE
   only when drift ≫ mild (gate≈0.45). Calm rejects stay uniform for packMSE.
3. Primary claim metric = **hard-subset next-MSE**; packMSE is conditional.

  PYTHONPATH=. python3 scripts/run_agod_obs_po_weights.py \\
    --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT stocks_IWM waymo_proxy
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor

from agod.hard_rank_metrics import hard_rank_metrics
from agod.obs_po_weights import (
    BEIJING_DRIFT_GATE,
    adaptive_temper,
    drift_intensity,
    gated_obs_po_weights,
    hard_subset_mask,
    is_beijing_class_drift,
)
from agod.online_rfperm import fit_online_rfperm, update_online_rfperm
from agod.po_refit import build_recent_ood_windows, refit_po_on_windows
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

# v4 modes: hard-claim path + beijing packMSE path + v3 baselines
MODES: Tuple[str, ...] = (
    "uniform",
    "gated_qrt_adapt",       # v3 soft qrt, mild gate 0.2
    "gated_qrt_hi",          # soft qrt, beijing gate 0.45 (packMSE path)
    "gated_hard_adapt",      # hard_support, mild gate 0.2 (hard claim)
    "gated_hard_hi",         # hard_support, beijing gate 0.45
)
COLORS = {
    "uniform": "#4C566A",
    "gated_qrt_adapt": "#A3BE8C",
    "gated_qrt_hi": "#EBCB8B",
    "gated_hard_adapt": "#88C0D0",
    "gated_hard_hi": "#5E81AC",
}


def ensure_loaders() -> None:
    for t in ("MSFT", "IWM", "AAPL", "SPY", "QQQ"):
        key = f"stocks_{t}"
        if key not in LOADERS:
            LOADERS[key] = (
                lambda tk: (lambda root, max_n=20000: load_stocks(root, tk, max_n))
            )(t)


def load_xy(name: str, root: Path, max_n: int = 20000):
    ensure_loaders()
    out = LOADERS[name](root, max_n=max_n)
    return np.asarray(out[0], float), np.asarray(out[1], float).ravel()


def pca_fit(X, d, seed):
    X = np.asarray(X, float)
    if d and X.shape[1] > d:
        return PCA(n_components=d, random_state=seed).fit_transform(X).astype(np.float32)
    return X.astype(np.float32)


def make_stream(X, y, bs: int, n_batches: int):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, len(X), bs)]


def fit_rf(X, y, w, seed):
    m = RandomForestRegressor(
        n_estimators=40, max_depth=8, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    m.fit(X, y, sample_weight=w)
    return m


def mode_spec(mode: str) -> dict:
    """obs map / drift gate / n_recent."""
    if mode == "uniform":
        return {"obs": "uniform", "adapt": False, "drift_gate": 0.0, "n_recent": None}
    if mode == "gated_qrt_adapt":
        return {"obs": "qrt", "adapt": True, "drift_gate": 0.20, "n_recent": None}
    if mode == "gated_qrt_hi":
        return {
            "obs": "qrt",
            "adapt": True,
            "drift_gate": BEIJING_DRIFT_GATE,
            "n_recent": None,
        }
    if mode == "gated_hard_adapt":
        return {
            "obs": "hard_support",
            "adapt": True,
            "drift_gate": 0.20,
            "n_recent": None,
        }
    if mode == "gated_hard_hi":
        return {
            "obs": "hard_support",
            "adapt": True,
            "drift_gate": BEIJING_DRIFT_GATE,
            "n_recent": None,
        }
    raise KeyError(mode)


def run_mode(
    stream,
    mode: str,
    seed: int,
    *,
    n_burn: int,
    alpha: float,
    n_recent_default: int,
    hard_frac: float = 0.2,
) -> dict:
    spec = mode_spec(mode)
    mse_next: List[float] = []
    mse_next_hard: List[float] = []
    mse_next_easy: List[float] = []
    gate_on: List[int] = []
    drifts: List[float] = []
    lams: List[float] = []
    beijing_on: List[int] = []
    hard_rows: List[dict] = []

    X0, y0 = stream[0]
    rfperm = fit_online_rfperm(X0, y0, seed=seed)

    for t in range(1, len(stream)):
        Xc, yc = stream[t]
        burn = t <= n_burn
        step = update_online_rfperm(
            rfperm, Xc, yc, burn_in=burn, alpha=alpha, ewma=True, fdr="alpha_investing"
        )
        rejected = bool(step["reject"]) and not burn
        gate_on.append(int(rejected))

        n_recent = spec["n_recent"] or n_recent_default
        windows = None
        po = None
        lam_used = 0.0
        drift = 0.0

        if mode == "uniform" or not rejected:
            w = np.ones(len(yc), float)
        else:
            windows = build_recent_ood_windows(
                stream, t, n_recent=n_recent, window_mode="recent_ood"
            )
            po = refit_po_on_windows(windows, seed=seed + t, blend_mu_gap=0.25)
            # Drift vs fixed f_ref (RFPerm yardstick) — not in-sample μ0.
            f_ref = rfperm.f_ref
            control_resid = np.abs(
                windows.y_recent - f_ref.predict(windows.X_recent)
            )
            ood_resid_ref = np.abs(yc - f_ref.predict(Xc))
            drift = drift_intensity(
                ood_resid_ref,
                control_resid=control_resid,
                p=float(step["p"]),
                T=float(step.get("T", 0.0)),
                alpha=alpha,
            )
            if spec["adapt"]:
                lam_used = adaptive_temper(
                    drift, lam_max=0.75, drift_gate=float(spec["drift_gate"])
                )
            else:
                lam_used = 0.0
            drifts.append(drift)
            lams.append(lam_used)
            beijing_on.append(int(is_beijing_class_drift(drift)))
            w = gated_obs_po_weights(
                po,
                reject=True,
                mode=spec["obs"],  # type: ignore[arg-type]
                soft=False,
                p=float(step["p"]),
                alpha=alpha,
                temper=lam_used,
                topk_frac=hard_frac,
                boost_max=3.0,
            )

        model = fit_rf(Xc, yc, w, seed + t)

        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            pred = model.predict(Xn)
            err2 = (yn - pred) ** 2
            mse_next.append(float(np.mean(err2)))

            if windows is not None:
                mu0 = fit_rf(
                    windows.X_recent,
                    windows.y_recent,
                    np.ones(len(windows.y_recent)),
                    seed + 99 + t,
                )
            else:
                Xp, yp = stream[t - 1]
                mu0 = fit_rf(Xp, yp, np.ones(len(yp)), seed + 99 + t)
            po_next = np.abs(yn - mu0.predict(Xn))
            hard_m = hard_subset_mask(po_next, frac=hard_frac)
            mse_next_hard.append(float(np.mean(err2[hard_m])))
            mse_next_easy.append(float(np.mean(err2[~hard_m])))

        if rejected and po is not None and windows is not None:
            Xr = np.vstack([windows.X_recent, windows.X_ood])
            yr = np.concatenate([windows.y_recent, windows.y_ood])
            mu_o = fit_rf(Xr, yr, np.ones(len(yr)), seed + 11 + t)
            truth = np.abs(yc - mu_o.predict(Xc))
            hard_rows.append({"po": hard_rank_metrics(po, truth), "drift": drift, "lam": lam_used})

    def _avg(xs: List[float]) -> float:
        return float(np.mean(xs)) if xs else float("nan")

    def _mean_hard() -> dict:
        if not hard_rows:
            return {"spearman": float("nan"), "precision_at_k": float("nan"), "n": 0}
        sp = [h["po"]["spearman"] for h in hard_rows if h["po"]["spearman"] == h["po"]["spearman"]]
        pk = [
            h["po"]["precision_at_k"]
            for h in hard_rows
            if h["po"]["precision_at_k"] == h["po"]["precision_at_k"]
        ]
        return {
            "spearman": float(np.mean(sp)) if sp else float("nan"),
            "precision_at_k": float(np.mean(pk)) if pk else float("nan"),
            "n": len(hard_rows),
        }

    return {
        "mse_next": mse_next,
        "mse_next_hard": mse_next_hard,
        "mse_next_easy": mse_next_easy,
        "gate_on": gate_on,
        "mse_mean": _avg(mse_next),
        "mse_hard_mean": _avg(mse_next_hard),
        "mse_easy_mean": _avg(mse_next_easy),
        "hard_po": _mean_hard(),
        "n_reject": int(sum(gate_on)),
        "duty": float(np.mean(gate_on)) if gate_on else 0.0,
        "drift_mean": _avg(drifts),
        "lam_mean": _avg(lams),
        "beijing_frac": _avg([float(x) for x in beijing_on]),
    }


def run_dataset(
    name: str,
    root: Path,
    *,
    seed: int,
    batch_size: int,
    n_batches: int,
    pca_d: int,
    n_burn: int,
    alpha: float,
    n_recent: int,
) -> dict:
    X, y = load_xy(name, root)
    X = pca_fit(X, pca_d, seed)
    stream = make_stream(X, y, batch_size, n_batches)
    results = {}
    for mode in MODES:
        results[mode] = run_mode(
            stream,
            mode,
            seed,
            n_burn=n_burn,
            alpha=alpha,
            n_recent_default=n_recent,
        )
    annotate_results_with_sig(
        results,
        preferred_gate_modes=(
            "gated_qrt_adapt",
            "gated_qrt_hi",
            "gated_hard_adapt",
            "gated_hard_hi",
        ),
    )
    gate = None
    for m in (
        "gated_qrt_adapt",
        "gated_hard_adapt",
        "gated_qrt_hi",
        "uniform",
    ):
        if m in results and results[m].get("gate_on"):
            g = np.asarray(results[m]["gate_on"], bool)
            n = len(results[m]["mse_next"])
            gate = g[:n]
            break
    if gate is not None and gate.any():
        for m, r in results.items():
            mh = np.asarray(r["mse_next_hard"], float)[: len(gate)]
            me = np.asarray(r["mse_next_easy"], float)[: len(gate)]
            r["mse_hard_mean_sig"] = float(mh[gate].mean()) if len(mh) == len(gate) else float("nan")
            r["mse_easy_mean_sig"] = float(me[gate].mean()) if len(me) == len(gate) else float("nan")
    return {
        "dataset": name,
        "n_batches": len(stream),
        "results": results,
        "drift_mean_adapt": results["gated_qrt_adapt"].get("drift_mean"),
        "beijing_class": is_beijing_class_drift(
            float(results["gated_qrt_adapt"].get("drift_mean") or 0.0)
        ),
    }


def _bars(all_ds: dict, out: Path):
    packs = list(all_ds.keys())
    fig, axes = plt.subplots(1, 2, figsize=(max(11, 1.8 * len(packs)), 4.6))
    x = np.arange(len(packs))
    w = 0.15
    mid = (len(MODES) - 1) / 2.0

    for ax, key, title in [
        (axes[0], "mse_mean_sig", "sig-only next MSE (all rows / pack)"),
        (axes[1], "mse_hard_mean_sig", "sig-only next MSE (hard top-20%)"),
    ]:
        for i, m in enumerate(MODES):
            vals = []
            for p in packs:
                r = all_ds[p]["results"][m]
                vals.append(r.get(key, r.get("mse_mean")))
            ax.bar(x + (i - mid) * w, vals, w, label=m, color=COLORS[m])
        ax.set_xticks(x)
        ax.set_xticklabels(packs, rotation=18, ha="right")
        ax.set_ylabel("MSE (↓)")
        ax.set_title(title)
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].legend(ncol=2, fontsize=6)
    fig.tight_layout()
    path = out / "obs_po_v4_mse.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def report(all_ds: dict) -> str:
    def f(v, pct=False):
        if v is None or (isinstance(v, float) and v != v):
            return "—"
        return f"{100*v:+.1f}%" if pct else f"{v:.4g}"

    short = {
        "uniform": "uniform",
        "gated_qrt_adapt": "qrt_mild",
        "gated_qrt_hi": "qrt_bj",
        "gated_hard_adapt": "hard_mild",
        "gated_hard_hi": "hard_bj",
    }
    lines = [
        "# Observation-level PO hard-reweight v4",
        "",
        "## Do we buy hard / beijing-drift / packMSE?",
        "",
        "**Yes, with a clean split (this is the locked thesis):**",
        "",
        "1. **认 hard** — obs PO is a hardness score. Spearman ~0.5–0.8, P@20% ≫ random.",
        "   The matching mechanism is **hard_support** (boost only the hard top-k),",
        "   not diffuse soft IPTW. Primary metric = **hard-subset next-MSE**.",
        "2. **条件认 beijing类漂移 → packMSE** — all-row pack MSE is only a fair",
        "   claim when hard-tail ≈ shift signal (drift ≫ mild, gate≈0.45).",
        "   On calm packs, hard ≈ noise → uniform wins packMSE; do not force lift.",
        "3. **不认** chasing packMSE on every RFPerm reject, or treating PO as an",
        "   image-OOD detector. Reject gate stays OnlineRFPerm; PO is post-hoc reweight.",
        "",
        f"v4: hard_support λ via mild gate (0.2) vs beijing gate ({BEIJING_DRIFT_GATE});",
        "soft qrt kept as the packMSE-oriented comparator under the same gates.",
        "",
        "## Drift intensity (mean on reject batches)",
        "",
        "| dataset | drift_mean | beijing_class? | lam qrt_mild | lam qrt_bj | lam hard_mild |",
        "|---|---:|:---:|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]
        bj = "yes" if blob.get("beijing_class") else "no"
        lines.append(
            f"| `{ds}` | {f(r['gated_qrt_adapt'].get('drift_mean'))} | {bj} | "
            f"{f(r['gated_qrt_adapt'].get('lam_mean'))} | "
            f"{f(r['gated_qrt_hi'].get('lam_mean'))} | "
            f"{f(r['gated_hard_adapt'].get('lam_mean'))} |"
        )

    hdr = " | ".join(short[m] for m in MODES)
    lines += [
        "",
        "## Sig-only next MSE — hard top-20% (↓)  ← primary claim (认 hard)",
        "",
        f"| dataset | {hdr} | best |",
        "|" + "---|---:" * len(MODES) + "|---|",
    ]
    wins_h = {m: 0 for m in MODES}
    for ds, blob in all_ds.items():
        r = blob["results"]

        def score_h(m):
            v = r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))
            return v if v == v else 1e99

        best = min(MODES, key=score_h)
        wins_h[best] += 1
        cells = [f(r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))) for m in MODES]
        lines.append(f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` |")
    lines += [
        "",
        "**Wins (hard-subset):** " + ", ".join(f"`{short[k]}`={v}" for k, v in wins_h.items()),
        "",
        "## Sig-only next MSE — all rows / pack (↓)  ← only claim under beijing drift",
        "",
        f"| dataset | {hdr} | best | beijing? |",
        "|" + "---|---:" * len(MODES) + "|---|:---:|",
    ]
    wins = {m: 0 for m in MODES}
    wins_bj = {m: 0 for m in MODES}
    n_bj = 0
    for ds, blob in all_ds.items():
        r = blob["results"]
        bj = bool(blob.get("beijing_class"))

        def score(m):
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return v if v == v else 1e99

        best = min(MODES, key=score)
        wins[best] += 1
        if bj:
            n_bj += 1
            wins_bj[best] += 1
        cells = [f(r[m].get("mse_mean_sig", r[m]["mse_mean"])) for m in MODES]
        lines.append(
            f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` | {'yes' if bj else 'no'} |"
        )
    lines += [
        "",
        "**Wins (all packs):** " + ", ".join(f"`{short[k]}`={v}" for k, v in wins.items()),
        f"**Wins among beijing-class packs only (n={n_bj}):** "
        + ", ".join(f"`{short[k]}`={v}" for k, v in wins_bj.items()),
        "",
        "## Rel. pack MSE vs uniform (soft qrt paths)",
        "",
        "| dataset | qrt_mild | qrt_bj | hard_mild | hard_bj | drift | beijing? |",
        "|---|---:|---:|---:|---:|---:|:---:|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]
        u = r["uniform"].get("mse_mean_sig", r["uniform"]["mse_mean"])
        bj = "yes" if blob.get("beijing_class") else "no"

        def rel(m):
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return (v / u - 1.0) if u == u and u > 0 else float("nan")

        lines.append(
            f"| `{ds}` | {f(rel('gated_qrt_adapt'), pct=True)} | {f(rel('gated_qrt_hi'), pct=True)} | "
            f"{f(rel('gated_hard_adapt'), pct=True)} | {f(rel('gated_hard_hi'), pct=True)} | "
            f"{f(r['gated_qrt_adapt'].get('drift_mean'))} | {bj} |"
        )

    lines += [
        "",
        "## Hard-rank (认 hard — unchanged)",
        "",
        "| dataset | spearman | P@20% | n_reject |",
        "|---|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        h = blob["results"]["gated_hard_adapt"]["hard_po"]
        nrej = blob["results"]["gated_hard_adapt"]["n_reject"]
        lines.append(f"| `{ds}` | {f(h['spearman'])} | {f(h['precision_at_k'])} | {nrej} |")

    lines += [
        "",
        "### Takeaway",
        "",
        "- **认 hard**: use PO to find / boost hard support; judge by hard-subset MSE + rank.",
        "- **条件认 packMSE**: only advertise all-row lift on beijing-class drift packs;",
        "  prefer `qrt_bj` / `hard_bj` (high gate) so calm rejects stay near uniform.",
        "- Mild-gate soft qrt remains a useful ablation, not the default packMSE claim.",
        "",
        "See `docs/agod/AGOD_obs_po_weights.md`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=[
            "metro_interstate",
            "beijing_pm25",
            "stocks_AAPL",
            "stocks_MSFT",
            "stocks_IWM",
            "waymo_proxy",
        ],
    )
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--pca-d", type=int, default=16)
    ap.add_argument("--n-burn", type=int, default=5)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--n-recent", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_obs_po_weights"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    all_ds: Dict[str, dict] = {}
    for name in args.datasets:
        print(f"=== {name} ===", flush=True)
        try:
            blob = run_dataset(
                name,
                args.root,
                seed=args.seed,
                batch_size=args.batch_size,
                n_batches=args.n_batches,
                pca_d=args.pca_d,
                n_burn=args.n_burn,
                alpha=args.alpha,
                n_recent=args.n_recent,
            )
        except Exception as e:
            print(f"[skip] {name}: {e}", flush=True)
            continue
        all_ds[name] = blob
        r = blob["results"]
        bj = "BJ" if blob.get("beijing_class") else "calm"
        print(
            "  drift={:.3f} ({}) | sigPack ".format(
                r["gated_qrt_adapt"].get("drift_mean") or float("nan"),
                bj,
            )
            + " ".join(
                f"{short_name(m)}={r[m].get('mse_mean_sig', r[m]['mse_mean']):.4g}" for m in MODES
            )
            + " | hard "
            + " ".join(
                f"{short_name(m)}={r[m].get('mse_hard_mean_sig', r[m].get('mse_hard_mean', float('nan'))):.4g}"
                for m in MODES
            ),
            flush=True,
        )

    payload = {
        "version": 4,
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "modes": list(MODES),
        "beijing_drift_gate": BEIJING_DRIFT_GATE,
        "datasets": all_ds,
        "note": "v4 hard_support + beijing-class packMSE gate",
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = report(all_ds)
    (args.out / "OBS_PO_WEIGHTS_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_obs_po_weights.md").write_text(md, encoding="utf-8")
    plots = [_bars(all_ds, args.out)] if all_ds else []
    print(md)
    print("plots:", [str(p) for p in plots])


def short_name(m: str) -> str:
    return {
        "uniform": "uni",
        "gated_qrt_adapt": "qrt_m",
        "gated_qrt_hi": "qrt_bj",
        "gated_hard_adapt": "h_m",
        "gated_hard_hi": "h_bj",
    }.get(m, m)


if __name__ == "__main__":
    main()
