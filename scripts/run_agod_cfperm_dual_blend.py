#!/usr/bin/env python3
"""CFPerm-gated dual/blend v8.1 — multi-seed + Jaccard(CFPerm, RFPerm).

Stack
-----
L0  CFPerm DRPerm(recent vs current) → reject?
L1  intensity (p, T, PO-gap)         → λ, beijing?
L2  shape: hard_support / soft CV / blend(mix)
L3  policy: dual | hard_m | blend_50

  PYTHONPATH=. python3 scripts/run_agod_cfperm_dual_blend.py \\
    --datasets metro_interstate beijing_pm25 stocks_AAPL stocks_MSFT stocks_IWM waymo_proxy \\
    --n-batches 40 --batch-size 256 --n-perm 39 --seeds 0 1 2 --jaccard-rfperm
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

from agod.cfperm_gate import (
    cfperm_batch_test,
    cfperm_intensity_to_temper,
    eval_cfperm_size_power,
    synthetic_shift_trial,
)
from agod.hard_rank_metrics import hard_rank_metrics
from agod.obs_po_cv import SOFT_POWERS, cv_select_power
from agod.obs_po_weights import (
    blend_hard_qrt_weights,
    gated_obs_po_weights,
    hard_subset_mask,
)
from agod.online_rfperm import fit_online_rfperm, update_online_rfperm
from agod.po_refit import build_recent_ood_windows, refit_po_on_windows
from agod.sig_batch_metrics import annotate_results_with_sig
from agod.stream_packs import LOADERS, load_stocks

DEFAULT_DATASETS = (
    "metro_interstate",
    "beijing_pm25",
    "stocks_AAPL",
    "stocks_MSFT",
    "stocks_IWM",
    "waymo_proxy",
)
MODES: Tuple[str, ...] = ("uniform", "hard_m", "dual", "blend_50")
COLORS = {
    "uniform": "#4C566A",
    "hard_m": "#88C0D0",
    "dual": "#5E81AC",
    "blend_50": "#D08770",
}


def ensure_loaders() -> None:
    for t in ("MSFT", "IWM", "AAPL", "SPY", "QQQ"):
        key = f"stocks_{t}"
        if key not in LOADERS:

            def _mk(tk: str):
                return lambda root, max_n=20000: load_stocks(root, tk, max_n)

            LOADERS[key] = _mk(t)


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


def short_name(m: str) -> str:
    return {"uniform": "uni", "hard_m": "hard_m", "dual": "dual", "blend_50": "b50"}[m]


def precompute_gates(
    stream,
    *,
    seed: int,
    n_burn: int,
    alpha: float,
    n_recent: int,
    n_perm: int,
    risk: str,
) -> List[dict]:
    cache: List[dict] = []
    for t in range(1, len(stream)):
        Xc, yc = stream[t]
        windows = build_recent_ood_windows(
            stream, t, n_recent=n_recent, window_mode="recent_ood"
        )
        if t <= n_burn:
            cache.append(
                {
                    "t": t,
                    "burn": True,
                    "reject": False,
                    "intensity": 0.0,
                    "lam": 0.0,
                    "is_bj": False,
                    "windows": windows,
                }
            )
            continue
        gate = cfperm_batch_test(
            windows.X_recent,
            windows.y_recent,
            Xc,
            yc,
            risk=risk,  # type: ignore[arg-type]
            n_perm=n_perm,
            alpha=alpha,
            seed=seed + 17 * t,
            e_mode="known",
        )
        # Store intensity only; λ / beijing remapped later (threshold scan).
        cache.append(
            {
                "t": t,
                "burn": False,
                "reject": bool(gate.reject),
                "intensity": float(gate.intensity),
                "lam": 0.0,
                "is_bj": False,
                "p_value": float(gate.p_value),
                "statistic": float(gate.statistic),
                "windows": windows,
            }
        )
    return cache


def remap_gate_temper(
    gate_cache: List[dict],
    *,
    beijing_gate: float = 0.25,
    temper_gate: float = 0.20,
    lam_max: float = 0.75,
) -> List[dict]:
    """Recompute (λ, is_beijing) from cached intensity without re-running CFPerm."""
    out: List[dict] = []
    for step in gate_cache:
        step = dict(step)
        if not step.get("burn", False):
            lam, is_bj = cfperm_intensity_to_temper(
                float(step["intensity"]),
                lam_max=lam_max,
                gate=temper_gate,
                beijing_gate=beijing_gate,
            )
            step["lam"] = float(lam)
            step["is_bj"] = bool(is_bj)
        out.append(step)
    return out


def rfperm_reject_mask(
    stream, *, n_burn: int, alpha: float, seed: int
) -> np.ndarray:
    X0, y0 = stream[0]
    state = fit_online_rfperm(X0, y0, seed=seed)
    out: List[int] = []
    for t in range(1, len(stream)):
        Xc, yc = stream[t]
        burn = t <= n_burn
        step = update_online_rfperm(
            state, Xc, yc, burn_in=burn, alpha=alpha, ewma=True, fdr="alpha_investing"
        )
        out.append(int(bool(step["reject"]) and not burn))
    return np.asarray(out, int)


def jaccard_binary(a, b) -> float:
    a = np.asarray(a, bool)
    b = np.asarray(b, bool)
    inter = int(np.sum(a & b))
    union = int(np.sum(a | b))
    return float(inter / union) if union else float("nan")


def run_mode(
    stream,
    mode: str,
    seed: int,
    gate_cache: List[dict],
    *,
    hard_frac: float,
    cv_folds: int,
) -> dict:
    mse_next: List[float] = []
    mse_hard: List[float] = []
    gate_on: List[int] = []
    intensities: List[float] = []
    lams: List[float] = []
    beijing_on: List[int] = []  # among rejects only
    beijing_all: List[int] = []  # all post-burn steps
    hard_rows: List[dict] = []
    fams: List[str] = []

    for step in gate_cache:
        t = int(step["t"])
        Xc, yc = stream[t]
        windows = step["windows"]

        if step["burn"]:
            gate_on.append(0)
            w = np.ones(len(yc), float)
        else:
            rejected = bool(step["reject"])
            gate_on.append(int(rejected))
            intensities.append(float(step["intensity"]))
            lam = float(step["lam"])
            is_bj = bool(step["is_bj"])
            lams.append(lam)
            if rejected:
                beijing_on.append(int(is_bj))
            else:
                # keep all-step beijing for diagnostics too
                pass
            beijing_all.append(int(is_bj))

            if mode == "uniform" or not rejected:
                w = np.ones(len(yc), float)
            else:
                po = refit_po_on_windows(windows, seed=seed + t, blend_mu_gap=0.25)
                use_hard = mode == "hard_m" or (
                    mode in ("dual", "blend_50") and not is_bj
                )
                if use_hard:
                    w = gated_obs_po_weights(
                        po,
                        reject=True,
                        mode="hard_support",
                        temper=max(lam, 0.25),
                        topk_frac=hard_frac,
                        boost_max=3.0,
                    )
                    fams.append("hard_support")
                elif mode == "dual":
                    cap = max(lam, 0.15)
                    sel = cv_select_power(
                        Xc,
                        yc,
                        po,
                        powers=SOFT_POWERS,
                        n_folds=cv_folds,
                        seed=seed + t,
                        temper_cap=cap,
                        objective="all",
                    )
                    w = np.asarray(sel["weights"], float)
                    fams.append(f"cv_PO^{float(sel['power']):g}")
                else:
                    w = blend_hard_qrt_weights(
                        po, lam=max(lam, 0.25), mix=0.5, topk_frac=hard_frac
                    )
                    fams.append("blend_50")

                Xr = np.vstack([windows.X_recent, windows.X_ood])
                yr = np.concatenate([windows.y_recent, windows.y_ood])
                mu_o = fit_rf(Xr, yr, np.ones(len(yr)), seed + 11 + t)
                truth = np.abs(yc - mu_o.predict(Xc))
                hard_rows.append({"po": hard_rank_metrics(po, truth)})

        model = fit_rf(Xc, yc, w, seed + t)
        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            err2 = (yn - model.predict(Xn)) ** 2
            mse_next.append(float(np.mean(err2)))
            Xp, yp = stream[t - 1]
            mu0 = fit_rf(Xp, yp, np.ones(len(yp)), seed + 99 + t)
            po_n = np.abs(yn - mu0.predict(Xn))
            hm = hard_subset_mask(po_n, frac=hard_frac)
            mse_hard.append(float(np.mean(err2[hm])))

    def _avg(xs):
        return float(np.mean(xs)) if xs else float("nan")

    fam_mode = "—"
    if fams:
        vals, counts = np.unique(np.asarray(fams), return_counts=True)
        fam_mode = str(vals[int(np.argmax(counts))])

    hard_po = {"spearman": float("nan"), "precision_at_k": float("nan"), "n": 0}
    if hard_rows:
        sp = [
            h["po"]["spearman"]
            for h in hard_rows
            if h["po"]["spearman"] == h["po"]["spearman"]
        ]
        pk = [
            h["po"]["precision_at_k"]
            for h in hard_rows
            if h["po"]["precision_at_k"] == h["po"]["precision_at_k"]
        ]
        hard_po = {
            "spearman": float(np.mean(sp)) if sp else float("nan"),
            "precision_at_k": float(np.mean(pk)) if pk else float("nan"),
            "n": len(hard_rows),
        }

    return {
        "mse_next": mse_next,
        "mse_next_hard": mse_hard,
        "gate_on": gate_on,
        "mse_mean": _avg(mse_next),
        "mse_hard_mean": _avg(mse_hard),
        "n_reject": int(sum(gate_on)),
        "duty": float(np.mean(gate_on)) if gate_on else 0.0,
        "intensity_mean": _avg(intensities),
        "lam_mean": _avg(lams),
        "beijing_frac": _avg([float(x) for x in beijing_on]),  # among rejects
        "beijing_frac_all": _avg([float(x) for x in beijing_all]),
        "family_mode": fam_mode,
        "hard_po": hard_po,
    }


def run_dataset(name, root, **kw) -> dict:
    pca_d = kw.pop("pca_d")
    batch_size = kw.pop("batch_size")
    n_batches = kw.pop("n_batches")
    seed = kw["seed"]
    hard_frac = kw.pop("hard_frac")
    cv_folds = kw.pop("cv_folds")
    jaccard_rfperm = bool(kw.pop("jaccard_rfperm", False))
    beijing_gate = float(kw.pop("beijing_gate", 0.25))
    temper_gate = float(kw.pop("temper_gate", 0.20))
    lam_max = float(kw.pop("lam_max", 0.75))
    modes = tuple(kw.pop("modes", MODES))

    X, y = load_xy(name, root)
    X = pca_fit(X, pca_d, seed)
    stream = make_stream(X, y, batch_size, n_batches)
    raw_cache = precompute_gates(
        stream,
        seed=seed,
        n_burn=kw["n_burn"],
        alpha=kw["alpha"],
        n_recent=kw["n_recent"],
        n_perm=kw["n_perm"],
        risk=kw["risk"],
    )
    gate_cache = remap_gate_temper(
        raw_cache,
        beijing_gate=beijing_gate,
        temper_gate=temper_gate,
        lam_max=lam_max,
    )
    results = {
        mode: run_mode(
            stream, mode, seed, gate_cache, hard_frac=hard_frac, cv_folds=cv_folds
        )
        for mode in modes
    }
    # pad missing modes for annotate if scanning with subset
    if "uniform" not in results:
        raise ValueError("uniform mode required for sig annotate")
    annotate_results_with_sig(
        results,
        preferred_gate_modes=tuple(
            m for m in ("hard_m", "dual", "blend_50") if m in results
        ),
    )
    gate = None
    for m in ("dual", "hard_m", "uniform"):
        if m in results and results[m].get("gate_on"):
            g = np.asarray(results[m]["gate_on"], bool)
            gate = g[: len(results[m]["mse_next"])]
            break
    if gate is not None and gate.any():
        for r in results.values():
            mh = np.asarray(r["mse_next_hard"], float)[: len(gate)]
            if len(mh) == len(gate):
                r["mse_hard_mean_sig"] = float(mh[gate].mean())

    jaccard = float("nan")
    rfperm_duty = float("nan")
    if jaccard_rfperm:
        cf = np.asarray(results["uniform"]["gate_on"], int)
        rf = rfperm_reject_mask(
            stream, n_burn=kw["n_burn"], alpha=kw["alpha"], seed=seed
        )
        n = min(len(cf), len(rf))
        jaccard = jaccard_binary(cf[:n], rf[:n])
        rfperm_duty = float(np.mean(rf[:n])) if n else float("nan")

    # intensity percentiles on reject steps (for threshold choice)
    inten_rej = [
        float(s["intensity"])
        for s in gate_cache
        if (not s.get("burn")) and s.get("reject")
    ]
    inten_q = {}
    if inten_rej:
        qs = [0.5, 0.75, 0.9]
        vals = np.quantile(inten_rej, qs)
        inten_q = {f"q{int(100*q)}": float(v) for q, v in zip(qs, vals)}
        inten_q["mean"] = float(np.mean(inten_rej))
        inten_q["max"] = float(np.max(inten_rej))
        inten_q["n_reject"] = len(inten_rej)

    dual = results.get("dual", results.get("hard_m", results["uniform"]))
    return {
        "dataset": name,
        "results": results,
        "duty": dual["duty"],
        "intensity_mean": dual["intensity_mean"],
        "beijing_frac": dual["beijing_frac"],
        "beijing_frac_all": dual.get("beijing_frac_all", float("nan")),
        "jaccard_cfperm_rfperm": jaccard,
        "rfperm_duty": rfperm_duty,
        "beijing_gate": beijing_gate,
        "intensity_reject_quantiles": inten_q,
        "seed": seed,
        "_raw_gate_cache": raw_cache,  # for scan reuse
        "_stream": stream,
    }


def _nanmean(xs: List[float]) -> float:
    arr = np.asarray(xs, float)
    arr = arr[np.isfinite(arr)]
    return float(arr.mean()) if len(arr) else float("nan")


def average_seed_blobs(blobs: List[dict]) -> dict:
    if len(blobs) == 1:
        out = {k: v for k, v in blobs[0].items() if not k.startswith("_")}
        out["n_seeds"] = 1
        out["seeds"] = [blobs[0].get("seed")]
        return out
    name = blobs[0]["dataset"]
    modes = list(blobs[0]["results"].keys())
    results: Dict[str, dict] = {}
    for mode in modes:
        keys = [
            "mse_mean",
            "mse_hard_mean",
            "mse_mean_sig",
            "mse_hard_mean_sig",
            "duty",
            "intensity_mean",
            "lam_mean",
            "beijing_frac",
            "beijing_frac_all",
            "n_reject",
            "rel_mse_vs_uniform_sig",
        ]
        agg = {
            k: _nanmean([b["results"][mode].get(k, float("nan")) for b in blobs])
            for k in keys
        }
        fams = [b["results"][mode].get("family_mode", "—") for b in blobs]
        vals, counts = np.unique(np.asarray(fams), return_counts=True)
        agg["family_mode"] = str(vals[int(np.argmax(counts))])
        agg["gate_on"] = blobs[0]["results"][mode].get("gate_on", [])
        agg["mse_next"] = blobs[0]["results"][mode].get("mse_next", [])
        sp = [
            b["results"][mode].get("hard_po", {}).get("spearman", float("nan"))
            for b in blobs
        ]
        pk = [
            b["results"][mode].get("hard_po", {}).get("precision_at_k", float("nan"))
            for b in blobs
        ]
        agg["hard_po"] = {
            "spearman": _nanmean(sp),
            "precision_at_k": _nanmean(pk),
            "n": int(
                np.nansum(
                    [b["results"][mode].get("hard_po", {}).get("n", 0) for b in blobs]
                )
            ),
        }
        results[mode] = agg
    return {
        "dataset": name,
        "results": results,
        "duty": results.get("dual", results[modes[0]])["duty"],
        "intensity_mean": results.get("dual", results[modes[0]])["intensity_mean"],
        "beijing_frac": results.get("dual", results[modes[0]])["beijing_frac"],
        "beijing_frac_all": results.get("dual", results[modes[0]]).get(
            "beijing_frac_all", float("nan")
        ),
        "jaccard_cfperm_rfperm": _nanmean(
            [b.get("jaccard_cfperm_rfperm", float("nan")) for b in blobs]
        ),
        "rfperm_duty": _nanmean([b.get("rfperm_duty", float("nan")) for b in blobs]),
        "beijing_gate": blobs[0].get("beijing_gate"),
        "n_seeds": len(blobs),
        "seeds": [b.get("seed") for b in blobs],
    }


def report(all_ds: dict, synth: dict | None = None) -> str:
    def f(v, pct=False):
        if v is None or (isinstance(v, float) and v != v):
            return "—"
        return f"{100 * v:+.1f}%" if pct else f"{v:.4g}"

    n_seeds = max((blob.get("n_seeds", 1) for blob in all_ds.values()), default=1)
    lines = [
        "# CFPerm-gated dual/blend (v8.3 full multi-seed (beijing_gate=0.25))",
        "",
        f"_Mean over **{n_seeds}** seed(s). L0 = CFPerm DRPerm (`e_mode=known`)._",
        "",
        "## Gate = CFPerm subset (not OnlineRFPerm)",
        "",
        "| piece | role |",
        "|---|---|",
        "| **DRPerm** (`risk=dr`) | L0 batch shift: PO-risk + permute-W |",
        "| **RRPerm** (`risk=rr`) | optional L0 via R-risk |",
        "| **CFPerm-VIMP** | post-hoc feature attribution (not stream gate) |",
        "",
        "L1 intensity = `0.55·po_gap + 0.30·p_strength + 0.15·T_strength`.",
        "dual: mild → hard_support; beijing (intensity>beijing_gate) → soft CV.",
        "blend_50: mild → hard; beijing → hard/qrt mix=0.5.",
        "",
    ]
    if synth:
        lines += [
            "## Synthetic gate eval (how we estimate/evaluate L0)",
            "",
            f"- **size** (null reject) = **{f(synth['size'])}**",
            f"- **power** (alt reject) = **{f(synth['power'])}**",
            f"- mean p null/alt = {f(synth['mean_p_null'])} / {f(synth['mean_p_alt'])}",
            f"- mean T null/alt = {f(synth['mean_stat_null'])} / {f(synth['mean_stat_alt'])}",
            "",
            "DGP: concept drift; evaluate size≈α, power↑.",
            "",
        ]
    lines += [
        "## Stream packs (CFPerm duty / intensity)",
        "",
        "| dataset | duty | intensitȳ | beijing_frac | dual fam |",
        "|---|---:|---:|---:|---|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]["dual"]
        lines.append(
            f"| `{ds}` | {f(r['duty'])} | {f(r['intensity_mean'])} | "
            f"{f(r['beijing_frac'])} | `{r['family_mode']}` |"
        )

    # Use modes present in results (scan may drop blend_50).
    modes = list(MODES)
    if all_ds:
        modes = [m for m in MODES if m in next(iter(all_ds.values()))["results"]]

    short = {m: short_name(m) for m in modes}
    hdr = " | ".join(short[m] for m in modes)
    lines += [
        "",
        "## Sig-only hard top-20% next-MSE (↓)",
        "",
        f"| dataset | {hdr} | best |",
        "|" + "---|---:" * len(modes) + "|---|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]

        def score_h(m, r=r):
            v = r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))
            return v if v == v else 1e99

        best = min(modes, key=score_h)
        cells = [
            f(r[m].get("mse_hard_mean_sig", r[m].get("mse_hard_mean"))) for m in modes
        ]
        lines.append(f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` |")

    lines += [
        "",
        "## Sig-only pack next-MSE (↓)",
        "",
        f"| dataset | {hdr} | best |",
        "|" + "---|---:" * len(modes) + "|---|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]

        def score(m, r=r):
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return v if v == v else 1e99

        best = min(modes, key=score)
        cells = [f(r[m].get("mse_mean_sig", r[m]["mse_mean"])) for m in modes]
        lines.append(f"| `{ds}` | " + " | ".join(cells) + f" | `{short[best]}` |")

    lines += [
        "",
        "## Rel. pack MSE vs uniform (CFPerm-sig)",
        "",
        "| dataset | hard_m | dual | b50 | duty |",
        "|---|---:|---:|---:|---:|",
    ]
    for ds, blob in all_ds.items():
        r = blob["results"]
        u = r["uniform"].get("mse_mean_sig", r["uniform"]["mse_mean"])

        def rel(m):
            if m not in r:
                return float("nan")
            v = r[m].get("mse_mean_sig", r[m]["mse_mean"])
            return (v / u - 1.0) if u == u and u > 0 else float("nan")

        lines.append(
            f"| `{ds}` | {f(rel('hard_m'), pct=True)} | {f(rel('dual'), pct=True)} | "
            f"{f(rel('blend_50'), pct=True)} | {f(r['dual']['duty'])} |"
        )

    if any(
        np.isfinite(blob.get("jaccard_cfperm_rfperm", float("nan")))
        for blob in all_ds.values()
    ):
        lines += [
            "",
            "## Jaccard(CFPerm reject, OnlineRFPerm reject)",
            "",
            "| dataset | Jaccard | CFPerm duty | RFPerm duty |",
            "|---|---:|---:|---:|",
        ]
        for ds, blob in all_ds.items():
            lines.append(
                f"| `{ds}` | {f(blob.get('jaccard_cfperm_rfperm'))} | "
                f"{f(blob['results']['dual']['duty'])} | "
                f"{f(blob.get('rfperm_duty'))} |"
            )

    lines += [
        "",
        "### Estimate / evaluate checklist",
        "",
        "1. **L0 estimate**: DRPerm on (recent, current); get p, T, reject.",
        "2. **L0 evaluate**: synthetic size/power; stream duty (selective).",
        "3. **L1 estimate**: intensity → λ / beijing.",
        "4. **L2/L3 evaluate**: sig-only hard-subset + pack MSE; calm ≈ uniform.",
        "5. **Ablation**: Jaccard(CFPerm, RFPerm) reject sets.",
        "",
        "See `docs/agod/AGOD_cfperm_dual_blend.md`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--n-batches", type=int, default=40)
    ap.add_argument("--pca-d", type=int, default=12)
    ap.add_argument("--n-burn", type=int, default=4)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--n-recent", type=int, default=1)
    ap.add_argument("--cv-folds", type=int, default=3)
    ap.add_argument("--n-perm", type=int, default=39)
    ap.add_argument("--risk", choices=["dr", "rr"], default="dr")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--seeds", type=int, nargs="+", default=None)
    ap.add_argument("--jaccard-rfperm", action="store_true")
    ap.add_argument("--skip-synth", action="store_true")
    ap.add_argument(
        "--beijing-gate",
        type=float,
        default=0.25,
        help="Intensity threshold for beijing soft path (v8.1 used 0.45 → rarely fired).",
    )
    ap.add_argument("--temper-gate", type=float, default=0.20)
    ap.add_argument("--lam-max", type=float, default=0.75)
    ap.add_argument(
        "--scan-beijing-gates",
        type=float,
        nargs="+",
        default=None,
        help="If set, reuse CFPerm cache and evaluate dual at each beijing gate.",
    )
    ap.add_argument("--out", type=Path, default=Path("results/agod_cfperm_dual_blend"))
    args = ap.parse_args()
    seeds = list(args.seeds) if args.seeds else [args.seed]

    min_p = 1.0 / (args.n_perm + 1)
    if min_p >= args.alpha:
        need = int(np.ceil(1.0 / args.alpha - 1.0))
        raise SystemExit(
            f"n_perm={args.n_perm} → min p={min_p:.4f} >= alpha={args.alpha}; "
            f"use n_perm>={need}."
        )
    args.out.mkdir(parents=True, exist_ok=True)

    synth = None
    if not args.skip_synth:
        print("=== synthetic CFPerm size/power ===", flush=True)
        trials = []
        for i in range(10):
            trials.append(
                synthetic_shift_trial(
                    shift=0.0, seed=i, n_perm=args.n_perm, n0=64, n1=64, p=6
                )
            )
        for i in range(10):
            trials.append(
                synthetic_shift_trial(
                    shift=1.5, seed=100 + i, n_perm=args.n_perm, n0=64, n1=64, p=6
                )
            )
        s = eval_cfperm_size_power(trials)
        synth = {
            "size": s.size,
            "power": s.power,
            "mean_p_null": s.mean_p_null,
            "mean_p_alt": s.mean_p_alt,
            "mean_stat_null": s.mean_stat_null,
            "mean_stat_alt": s.mean_stat_alt,
        }
        print(f"  size={s.size:.3f} power={s.power:.3f}", flush=True)

    all_ds: Dict[str, dict] = {}
    scan_rows: List[dict] = []
    scan_gates = list(args.scan_beijing_gates) if args.scan_beijing_gates else None

    for name in args.datasets:
        print(f"=== {name} seeds={seeds} ===", flush=True)
        seed_blobs: List[dict] = []
        for sd in seeds:
            try:
                blob = run_dataset(
                    name,
                    args.root,
                    seed=sd,
                    batch_size=args.batch_size,
                    n_batches=args.n_batches,
                    pca_d=args.pca_d,
                    n_burn=args.n_burn,
                    alpha=args.alpha,
                    n_recent=args.n_recent,
                    hard_frac=0.2,
                    cv_folds=args.cv_folds,
                    n_perm=args.n_perm,
                    risk=args.risk,
                    jaccard_rfperm=args.jaccard_rfperm,
                    beijing_gate=args.beijing_gate,
                    temper_gate=args.temper_gate,
                    lam_max=args.lam_max,
                    modes=("uniform", "hard_m", "dual", "blend_50")
                    if not scan_gates
                    else ("uniform", "hard_m", "dual"),
                )
            except Exception as e:
                print(f"  [skip] seed={sd}: {e}", flush=True)
                continue
            seed_blobs.append(blob)
            print(
                "  seed={}: duty={:.2f} inten={:.3f} bj_rej={:.2f} fam={} q={}".format(
                    sd,
                    blob["results"]["dual"]["duty"],
                    blob["results"]["dual"]["intensity_mean"],
                    blob["results"]["dual"].get("beijing_frac", float("nan")),
                    blob["results"]["dual"].get("family_mode"),
                    blob.get("intensity_reject_quantiles"),
                ),
                flush=True,
            )

            # Threshold scan: reuse CFPerm cache, rematerialize dual only
            if scan_gates and "_raw_gate_cache" in blob:
                stream = blob["_stream"]
                raw = blob["_raw_gate_cache"]
                for bg in scan_gates:
                    gcache = remap_gate_temper(
                        raw,
                        beijing_gate=float(bg),
                        temper_gate=args.temper_gate,
                        lam_max=args.lam_max,
                    )
                    dual_r = run_mode(
                        stream,
                        "dual",
                        sd,
                        gcache,
                        hard_frac=0.2,
                        cv_folds=args.cv_folds,
                    )
                    uni_r = run_mode(
                        stream,
                        "uniform",
                        sd,
                        gcache,
                        hard_frac=0.2,
                        cv_folds=args.cv_folds,
                    )
                    # sig pack MSE on CFPerm rejects
                    gmask = np.asarray(dual_r["gate_on"], bool)
                    n = min(len(gmask), len(dual_r["mse_next"]), len(uni_r["mse_next"]))
                    gmask = gmask[:n]
                    if gmask.any():
                        dual_sig = float(np.mean(np.asarray(dual_r["mse_next"])[:n][gmask]))
                        uni_sig = float(np.mean(np.asarray(uni_r["mse_next"])[:n][gmask]))
                        rel = dual_sig / uni_sig - 1.0 if uni_sig > 0 else float("nan")
                    else:
                        dual_sig = uni_sig = rel = float("nan")
                    scan_rows.append(
                        {
                            "dataset": name,
                            "seed": sd,
                            "beijing_gate": float(bg),
                            "duty": dual_r["duty"],
                            "beijing_frac_reject": dual_r["beijing_frac"],
                            "family_mode": dual_r["family_mode"],
                            "dual_mse_sig": dual_sig,
                            "uni_mse_sig": uni_sig,
                            "rel_vs_uni": rel,
                        }
                    )
                    print(
                        f"    scan bj_gate={bg:.2f}: bj_rej={dual_r['beijing_frac']:.2f} "
                        f"fam={dual_r['family_mode']} rel={rel:+.1%}"
                        if rel == rel
                        else f"    scan bj_gate={bg:.2f}: bj_rej={dual_r['beijing_frac']:.2f} fam={dual_r['family_mode']} rel=—",
                        flush=True,
                    )

        if not seed_blobs:
            continue
        blob = average_seed_blobs(seed_blobs)
        all_ds[name] = blob
        r = blob["results"]
        modes_present = [m for m in MODES if m in r]
        print(
            "  mean duty={:.2f} inten={:.3f} bj_rej={:.2f} | pack ".format(
                r["dual"]["duty"],
                r["dual"]["intensity_mean"],
                r["dual"].get("beijing_frac", float("nan")),
            )
            + " ".join(
                f"{short_name(m)}={r[m].get('mse_mean_sig', r[m]['mse_mean']):.4g}"
                for m in modes_present
            ),
            flush=True,
        )

    # Summarize beijing-gate scan
    scan_summary = None
    if scan_rows:
        scan_summary = {}
        for bg in sorted(set(r["beijing_gate"] for r in scan_rows)):
            rows = [r for r in scan_rows if r["beijing_gate"] == bg]
            scan_summary[str(bg)] = {
                "beijing_frac_reject_mean": _nanmean(
                    [r["beijing_frac_reject"] for r in rows]
                ),
                "rel_vs_uni_mean": _nanmean([r["rel_vs_uni"] for r in rows]),
                "n": len(rows),
                "by_dataset": {},
            }
            for ds in sorted(set(r["dataset"] for r in rows)):
                drows = [r for r in rows if r["dataset"] == ds]
                scan_summary[str(bg)]["by_dataset"][ds] = {
                    "beijing_frac_reject": _nanmean(
                        [r["beijing_frac_reject"] for r in drows]
                    ),
                    "rel_vs_uni": _nanmean([r["rel_vs_uni"] for r in drows]),
                    "family_modes": sorted(
                        {r["family_mode"] for r in drows if r["family_mode"]}
                    ),
                }
        print("=== beijing_gate scan summary ===", flush=True)
        for bg, s in scan_summary.items():
            print(
                f"  gate={bg}: bj_rej̄={s['beijing_frac_reject_mean']:.3f} "
                f"rel̄={s['rel_vs_uni_mean']:+.1%}"
                if s["rel_vs_uni_mean"] == s["rel_vs_uni_mean"]
                else f"  gate={bg}: bj_rej̄={s['beijing_frac_reject_mean']:.3f} rel̄=—",
                flush=True,
            )
        (args.out / "beijing_gate_scan.json").write_text(
            json.dumps({"rows": scan_rows, "summary": scan_summary}, indent=2),
            encoding="utf-8",
        )

    md = report(all_ds, synth)
    if scan_summary:
        md += "\n## Beijing-gate scan (dual soft path)\n\n"
        md += "| beijing_gate | bj_frac@reject̄ | dual rel pack vs unī |\n|---|---:|---:|\n"
        for bg, s in scan_summary.items():
            rel = s["rel_vs_uni_mean"]
            md += (
                f"| {bg} | {s['beijing_frac_reject_mean']:.3f} | "
                + (f"{100*rel:+.1f}%" if rel == rel else "—")
                + " |\n"
            )
        md += (
            "\n_Lower beijing_gate → more soft-CV dual path. "
            "v8.1 default 0.45 rarely fired; retune toward reject-intensity quantiles._\n"
        )
    (args.out / "CFPERM_DUAL_BLEND_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_cfperm_dual_blend.md").write_text(md, encoding="utf-8")

    slim = {}
    for ds, blob in all_ds.items():
        slim[ds] = {
            "duty": blob["duty"],
            "intensity_mean": blob["intensity_mean"],
            "beijing_frac": blob["beijing_frac"],
            "beijing_frac_all": blob.get("beijing_frac_all"),
            "beijing_gate": blob.get("beijing_gate", args.beijing_gate),
            "jaccard_cfperm_rfperm": blob.get("jaccard_cfperm_rfperm"),
            "rfperm_duty": blob.get("rfperm_duty"),
            "n_seeds": blob.get("n_seeds", 1),
            "seeds": blob.get("seeds"),
            "results": {
                m: {
                    k: v
                    for k, v in blob["results"][m].items()
                    if k not in ("mse_next", "mse_next_hard", "gate_on")
                }
                for m in blob["results"]
            },
        }
    payload = {
        "version": "8.3",
        "gate": "cfperm",
        "risk": args.risk,
        "n_perm": args.n_perm,
        "seeds": seeds,
        "beijing_gate": args.beijing_gate,
        "temper_gate": args.temper_gate,
        "lam_max": args.lam_max,
        "jaccard_rfperm": args.jaccard_rfperm,
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "modes": list(MODES),
        "synth": synth,
        "beijing_gate_scan": scan_summary,
        "datasets": slim,
    }
    (args.out / "summary.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8"
    )

    if all_ds:
        packs = list(all_ds.keys())
        modes_plot = [m for m in MODES if m in next(iter(all_ds.values()))["results"]]
        fig, ax = plt.subplots(figsize=(max(8, 1.6 * len(packs)), 4.2))
        x = np.arange(len(packs))
        w = 0.18
        mid = (len(modes_plot) - 1) / 2.0
        for i, m in enumerate(modes_plot):
            vals = [
                all_ds[p]["results"][m].get(
                    "mse_mean_sig", all_ds[p]["results"][m]["mse_mean"]
                )
                for p in packs
            ]
            ax.bar(x + (i - mid) * w, vals, w, label=short_name(m), color=COLORS[m])
        ax.set_xticks(x)
        ax.set_xticklabels(packs, rotation=15, ha="right")
        ax.set_ylabel("sig pack MSE (↓)")
        ax.set_title(
            f"CFPerm dual/blend (seeds={seeds}, beijing_gate={args.beijing_gate})"
        )
        ax.legend(fontsize=8)
        ax.grid(True, axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(args.out / "cfperm_dual_blend_mse.png", dpi=140)
        plt.close(fig)
    print(md)


if __name__ == "__main__":
    main()
