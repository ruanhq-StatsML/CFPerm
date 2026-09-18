#!/usr/bin/env python3
"""SAFFRON / ADDIS / α-investing on OnlineRFPerm streams — manuscript form.

Report form (OnlinePermOOB / OnlineRFPerm.pdf):

  first1 / first2 / first3   — first k consecutive rejects (1-based end idx)
  alarm P25 / median / P75   — distribution of *all* alarm times on the trail
  day DetRate                — e.g. NYC-taxi: n_alarm_days / n_days

  PYTHONPATH=. python3 scripts/run_saffron_addis_realdata.py \\
    --datasets all --seeds 0 1 2 --addis-lambdas 0.25 0.35 0.40
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
import pandas as pd
from online_fdr import Addis, Saffron
from sklearn.preprocessing import StandardScaler

from agod.first_k_metrics import (
    aggregate_day_forms,
    first_k_form,
    fmt_num,
)
from agod.online_rfperm import batch_T, fit_online_rfperm, rank_pvalue
from agod.stream_packs import LOADERS as PACK_LOADERS

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------


def load_synthetic(max_n: int, seed: int, *, shift_batch: int = 20, bs: int = 128):
    rng = np.random.default_rng(seed)
    d = 16
    n = max_n
    X = rng.normal(size=(n, d)).astype(np.float32)
    y = np.zeros(n, dtype=np.float64)
    t_cut = int(shift_batch * bs)
    t_cut = max(bs * 4, min(t_cut, n - bs * 4))
    w0 = rng.normal(size=d)
    w0 /= np.linalg.norm(w0) + 1e-9
    w1 = rng.normal(size=d)
    w1 /= np.linalg.norm(w1) + 1e-9
    y[:t_cut] = X[:t_cut] @ w0 + rng.normal(0, 0.2, t_cut)
    X[t_cut:] = X[t_cut:] * 1.6 + 0.75
    y[t_cut:] = X[t_cut:] @ w1 * 2.2 + rng.normal(0, 0.45, n - t_cut)
    mu = X[:t_cut].mean(0, keepdims=True)
    sd = X[:t_cut].std(0, keepdims=True) + 1e-6
    X = ((X - mu) / sd).astype(np.float32)
    return X, y, {"shift_batch": int(t_cut // bs), "name": "synthetic"}


def load_openml(name: str, max_n: int, seed: int):
    from sklearn.datasets import fetch_openml, fetch_covtype
    from sklearn.preprocessing import LabelEncoder

    if name == "covertype":
        bun = fetch_covtype()
        X = np.asarray(bun.data, float)
        y = (bun.target == 1).astype(float)
    else:
        bun = fetch_openml(name, version=1, as_frame=True, parser="auto")
        df = bun.data.copy()
        for c in df.columns:
            if df[c].dtype.kind in "OSUb":
                df[c] = LabelEncoder().fit_transform(df[c].astype(str).fillna("NA"))
            else:
                df[c] = df[c].astype(float)
        X = df.to_numpy(float)
        y = bun.target
        if getattr(y, "dtype", None) is not None and y.dtype.kind in "OSU":
            y = LabelEncoder().fit_transform(y.astype(str)).astype(float)
        else:
            y = np.asarray(y, float)
            if np.unique(y).size > 20:
                y = (y > np.median(y)).astype(float)
    rng = np.random.default_rng(seed)
    n = min(max_n, len(X))
    start = int(rng.integers(0, max(1, len(X) - n)))
    X = X[start : start + n]
    y = y[start : start + n]
    mask = np.isfinite(X).all(1) & np.isfinite(y)
    X, y = X[mask], y[mask]
    X = StandardScaler().fit_transform(X).astype(np.float32)
    return X, y, {"name": name, "shift_batch": None}


def load_pack(name: str, max_n: int, seed: int):
    if name not in PACK_LOADERS:
        raise KeyError(name)
    X, y, meta = PACK_LOADERS[name](ROOT, max_n=max_n)
    # time-ordered packs: do not shuffle; scale on prefix
    n = len(X)
    fit_n = max(64, n // 5)
    sc = StandardScaler().fit(X[:fit_n])
    X = sc.transform(X).astype(np.float32)
    meta = dict(meta)
    meta["shift_batch"] = None
    meta["seed"] = seed
    return X, y, meta


OPENML_LOADERS = {
    "synthetic": lambda max_n, seed, bs: load_synthetic(max_n, seed, bs=bs),
    "covertype": lambda max_n, seed, bs: load_openml("covertype", max_n, seed),
    "bank": lambda max_n, seed, bs: load_openml("bank-marketing", max_n, seed),
    "electricity": lambda max_n, seed, bs: load_openml("electricity", max_n, seed),
    "eeg": lambda max_n, seed, bs: load_openml("eeg-eye-state", max_n, seed),
    "adult": lambda max_n, seed, bs: load_openml("adult", max_n, seed),
}

PACK_NAMES = [
    "metro_interstate",
    "beijing_pm25",
    "nyc_taxi",
    "stocks_AAPL",
    "stocks_MSFT",
    "stocks_IWM",
    "waymo_proxy",
]

ALL_DATASETS = list(OPENML_LOADERS.keys()) + PACK_NAMES


def load_any(name: str, max_n: int, seed: int, bs: int):
    if name in OPENML_LOADERS:
        return OPENML_LOADERS[name](max_n, seed, bs)
    return load_pack(name, max_n, seed)


def make_stream(X, y, bs: int, n_batches: int):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    if len(X) < need:
        raise ValueError(f"need {need}, got {len(X)}")
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, need, bs)]


# ---------------------------------------------------------------------------
# p-stream + FDR
# ---------------------------------------------------------------------------


def p_stream_from_rfperm(
    stream,
    *,
    n_burn: int,
    seed: int,
    ewma_lam: float = 1.0,
) -> Tuple[List[float], List[float]]:
    X0 = np.vstack([stream[t][0] for t in range(n_burn)])
    y0 = np.concatenate([stream[t][1] for t in range(n_burn)])
    st = fit_online_rfperm(X0, y0, seed=seed)
    T_hist: List[float] = []
    for t in range(n_burn):
        T = batch_T(st, stream[t][0], stream[t][1])
        T_hist.append(T)
    ps: List[float] = []
    Ts: List[float] = list(T_hist)
    for t in range(n_burn, len(stream)):
        T = batch_T(st, stream[t][0], stream[t][1])
        p = rank_pvalue(T, T_hist, ewma=True, lam=ewma_lam)
        ps.append(float(p))
        T_hist.append(T)
        Ts.append(T)
    return Ts, ps


def rejects_from_procedure(ps: Sequence[float], name: str, **kw) -> List[int]:
    if name == "saffron":
        proc = Saffron(alpha=kw["alpha"], wealth=kw["wealth"], lambda_=kw["lambda_"])
        return [int(proc.test_one(float(p))) for p in ps]
    if name.startswith("addis"):
        proc = Addis(
            alpha=kw["alpha"],
            wealth=kw["wealth"],
            lambda_=kw["lambda_"],
            tau=kw["tau"],
        )
        return [int(proc.test_one(float(p))) for p in ps]
    if name == "alpha_investing":
        wealth = 1.0
        alpha = kw["alpha"]
        rejects = []
        for p in ps:
            alpha_t = max(wealth, 1e-8) * alpha / (1.0 + alpha)
            rej = bool(p < alpha_t)
            wealth = wealth + alpha if rej else max(wealth - alpha_t, 1e-8)
            rejects.append(int(rej))
        return rejects
    raise ValueError(name)


def apply_grace(rejects: List[int], grace: int) -> List[int]:
    scored = list(rejects)
    for i in range(min(max(0, grace), len(scored))):
        scored[i] = 0
    return scored


def run_one_stream(
    X,
    y,
    *,
    bs: int,
    n_batches: int,
    n_burn: int,
    seed: int,
    procedures: List[Tuple[str, dict]],
    alpha: float,
    wealth: float,
    addis_tau: float,
    saffron_lambda: float,
    grace: int,
    shift_batch: Optional[int],
) -> List[Dict]:
    need = bs * n_batches
    if len(X) < need:
        nb = max(n_burn + 4, len(X) // bs)
    else:
        nb = n_batches
    stream = make_stream(X, y, bs, nb)
    n_burn_eff = min(n_burn, max(2, nb // 4))
    _, ps = p_stream_from_rfperm(stream, n_burn=n_burn_eff, seed=seed)
    shift_rel = None if shift_batch is None else int(shift_batch - n_burn_eff)
    rows = []
    for pname, extra in procedures:
        kw = {
            "alpha": alpha,
            "wealth": wealth,
            "tau": addis_tau,
            "lambda_": saffron_lambda,
        }
        kw.update(extra)
        rejects = rejects_from_procedure(
            ps, "addis" if pname.startswith("addis") else pname, **kw
        )
        scored = apply_grace(rejects, grace)
        form = first_k_form(scored)
        # also raw (no grace) for SUM / day alarms
        form_raw = first_k_form(rejects)
        delay = None
        if shift_rel is not None and np.isfinite(form["first1"]):
            # first1 is 1-based end idx on post-burn trail
            delay = float(form["first1"] - 1 - shift_rel)
        rows.append(
            {
                "procedure": pname,
                "lambda": extra.get("lambda_", None),
                "n_burn": n_burn_eff,
                "n_batches": nb,
                "n_monitor": len(rejects),
                "shift_batch": shift_batch,
                "grace": grace,
                "delay_vs_shift": delay,
                **{f"g_{k}": v for k, v in form.items()},
                **{f"raw_{k}": v for k, v in form_raw.items()},
            }
        )
    return rows


def day_slices(
    X: np.ndarray,
    y: np.ndarray,
    day_id: np.ndarray,
    *,
    min_len: int,
    max_days: Optional[int] = None,
) -> List[Tuple[np.ndarray, np.ndarray]]:
    out = []
    for d in np.unique(day_id):
        m = day_id == d
        if int(m.sum()) < min_len:
            continue
        out.append((X[m], y[m]))
        if max_days is not None and len(out) >= max_days:
            break
    return out


def run_day_level(
    X,
    y,
    day_id,
    *,
    bs: int,
    n_burn: int,
    seed: int,
    procedures: List[Tuple[str, dict]],
    alpha: float,
    wealth: float,
    addis_tau: float,
    saffron_lambda: float,
    grace: int,
    max_days: Optional[int],
) -> List[Dict]:
    """Each calendar day = one stream unit → DetRate + first1 quantiles."""
    slices = day_slices(X, y, day_id, min_len=max(bs * (n_burn + 2), bs * 4), max_days=max_days)
    rows = []
    for pname, extra in procedures:
        day_forms = []
        for di, (Xd, yd) in enumerate(slices):
            nb = len(Xd) // bs
            if nb < n_burn + 2:
                continue
            sub = run_one_stream(
                Xd,
                yd,
                bs=bs,
                n_batches=nb,
                n_burn=n_burn,
                seed=seed + di,
                procedures=[(pname, extra)],
                alpha=alpha,
                wealth=wealth,
                addis_tau=addis_tau,
                saffron_lambda=saffron_lambda,
                grace=grace,
                shift_batch=None,
            )
            if not sub:
                continue
            r = sub[0]
            day_forms.append(
                {
                    "SUM": r["g_SUM"],
                    "first1": r["g_first1"],
                    "first2": r["g_first2"],
                    "first3": r["g_first3"],
                }
            )
        agg = aggregate_day_forms(day_forms)
        rows.append(
            {
                "procedure": pname,
                "lambda": extra.get("lambda_", None),
                "mode": "day_level",
                **agg,
            }
        )
    return rows


def build_procedures(saffron_lambda: float, addis_lambdas: Sequence[float]):
    return [
        ("saffron", {"lambda_": saffron_lambda}),
        *[(f"addis_lam{lam:g}", {"lambda_": lam}) for lam in addis_lambdas],
        ("alpha_investing", {}),
    ]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["all"],
        help="dataset names or 'all'",
    )
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--n-batches", type=int, default=48)
    ap.add_argument("--n-burn", type=int, default=8)
    ap.add_argument("--max-n", type=int, default=12000)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--saffron-lambda", type=float, default=0.5)
    ap.add_argument(
        "--addis-lambdas",
        type=float,
        nargs="+",
        default=[0.25, 0.35, 0.40],
    )
    ap.add_argument("--addis-tau", type=float, default=0.5)
    ap.add_argument("--grace", type=int, default=4)
    ap.add_argument(
        "--day-batch-size",
        type=int,
        default=8,
        help="Batch size inside each calendar-day stream (NYC/metro/beijing)",
    )
    ap.add_argument("--day-burn", type=int, default=2)
    ap.add_argument("--max-days", type=int, default=None)
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "saffron_addis_realdata",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    datasets = ALL_DATASETS if "all" in args.datasets else list(args.datasets)
    wealth = args.alpha / 2.0
    procedures = build_procedures(args.saffron_lambda, args.addis_lambdas)

    stream_rows: List[Dict] = []
    day_rows: List[Dict] = []

    for ds in datasets:
        print(f"######## dataset={ds} ########", flush=True)
        for seed in args.seeds:
            print(f"=== {ds} seed={seed} ===", flush=True)
            try:
                X, y, meta = load_any(ds, args.max_n, seed, args.batch_size)
            except Exception as e:
                print(f"  SKIP load {ds}: {e}", flush=True)
                continue
            shift_batch = meta.get("shift_batch")
            try:
                rows = run_one_stream(
                    X,
                    y,
                    bs=args.batch_size,
                    n_batches=args.n_batches,
                    n_burn=args.n_burn,
                    seed=seed,
                    procedures=procedures,
                    alpha=args.alpha,
                    wealth=wealth,
                    addis_tau=args.addis_tau,
                    saffron_lambda=args.saffron_lambda,
                    grace=args.grace,
                    shift_batch=shift_batch,
                )
            except Exception as e:
                print(f"  SKIP stream {ds}: {e}", flush=True)
                continue
            for r in rows:
                stream_rows.append({"dataset": ds, "seed": seed, "mode": "stream", **r})
                print(
                    f"  {r['procedure']}: first1={fmt_num(r['g_first1'])} "
                    f"first2={fmt_num(r['g_first2'])} first3={fmt_num(r['g_first3'])} "
                    f"P25/med/P75={fmt_num(r['g_alarm_p25'])}/"
                    f"{fmt_num(r['g_alarm_median'])}/{fmt_num(r['g_alarm_p75'])} "
                    f"SUM={r['g_SUM']}",
                    flush=True,
                )

            # day-level for packs with day_id (NYC / metro / beijing)
            day_id = meta.get("day_id")
            if day_id is not None and seed == args.seeds[0]:
                print(f"  -- day-level ({meta.get('n_days')} calendar days) --", flush=True)
                try:
                    drows = run_day_level(
                        X,
                        y,
                        np.asarray(day_id),
                        bs=args.day_batch_size,
                        n_burn=args.day_burn,
                        seed=seed,
                        procedures=procedures,
                        alpha=args.alpha,
                        wealth=wealth,
                        addis_tau=args.addis_tau,
                        saffron_lambda=args.saffron_lambda,
                        grace=max(0, args.grace // 2),
                        max_days=args.max_days,
                    )
                    for r in drows:
                        day_rows.append({"dataset": ds, **r})
                        print(
                            f"  day {r['procedure']}: "
                            f"{r['n_alarm_days']}/{r['n_days']} days alarmed "
                            f"(DetRate={r['det_rate']:.3f}) "
                            f"first1 med={fmt_num(r['first1_median'])} "
                            f"P25/P75={fmt_num(r['first1_p25'])}/{fmt_num(r['first1_p75'])}",
                            flush=True,
                        )
                except Exception as e:
                    print(f"  SKIP day-level {ds}: {e}", flush=True)

    df = pd.DataFrame(stream_rows)
    df.to_csv(args.out_dir / "saffron_addis_runs.csv", index=False)

    # ---- manuscript form tables ----
    form_cols = [
        "dataset",
        "procedure",
        "g_first1",
        "g_first2",
        "g_first3",
        "g_alarm_p25",
        "g_alarm_median",
        "g_alarm_p75",
        "g_alarm_mean",
        "g_alarm_std",
        "g_SUM",
        "g_n_trail",
    ]
    if len(df):
        form_agg = (
            df.groupby(["dataset", "procedure"], as_index=False)
            .agg(
                first1=("g_first1", "median"),
                first2=("g_first2", "median"),
                first3=("g_first3", "median"),
                alarm_p25=("g_alarm_p25", "median"),
                alarm_median=("g_alarm_median", "median"),
                alarm_p75=("g_alarm_p75", "median"),
                alarm_mean=("g_alarm_mean", "mean"),
                alarm_std=("g_alarm_std", "mean"),
                SUM=("g_SUM", "mean"),
                n_trail=("g_n_trail", "first"),
                lambda_=("lambda", "first"),
            )
            .sort_values(["dataset", "procedure"])
        )
        form_agg.to_csv(args.out_dir / "saffron_addis_firstk_form.csv", index=False)
    else:
        form_agg = pd.DataFrame()

    day_df = pd.DataFrame(day_rows)
    if len(day_df):
        day_df.to_csv(args.out_dir / "saffron_addis_day_detrate.csv", index=False)

    overall = (
        df.groupby("procedure", as_index=False)
        .agg(
            first1_med=("g_first1", "median"),
            first2_med=("g_first2", "median"),
            first3_med=("g_first3", "median"),
            alarm_p25=("g_alarm_p25", "median"),
            alarm_median=("g_alarm_median", "median"),
            alarm_p75=("g_alarm_p75", "median"),
            mean_SUM=("g_SUM", "mean"),
            lambda_=("lambda", "first"),
        )
        .sort_values("procedure")
        if len(df)
        else pd.DataFrame()
    )
    if len(overall):
        overall.to_csv(args.out_dir / "saffron_addis_overall.csv", index=False)

    # plot: first1 + alarm median across procedures
    if len(overall):
        procs_order = ["saffron"] + [
            f"addis_lam{lam:g}" for lam in args.addis_lambdas
        ] + ["alpha_investing"]
        procs_order = [p for p in procs_order if p in set(overall["procedure"])]
        o = overall.set_index("procedure").loc[procs_order]
        fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
        ax = axes[0]
        ax.bar(range(len(o)), o["first1_med"], color="#4C78A8")
        ax.set_xticks(range(len(o)))
        ax.set_xticklabels(procs_order, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("median first1")
        ax.set_title("first1 (first reject)")
        ax = axes[1]
        ax.bar(range(len(o)), o["alarm_median"], color="#F58518", label="median")
        ax.plot(range(len(o)), o["alarm_p25"], "v", color="#54A24B", label="P25")
        ax.plot(range(len(o)), o["alarm_p75"], "^", color="#E45756", label="P75")
        ax.set_xticks(range(len(o)))
        ax.set_xticklabels(procs_order, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("alarm time")
        ax.set_title("Alarm-time P25 / median / P75")
        ax.legend(fontsize=7)
        ax = axes[2]
        ax.bar(range(len(o)), o["mean_SUM"], color="#B279A2")
        ax.set_xticks(range(len(o)))
        ax.set_xticklabels(procs_order, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("mean SUM")
        ax.set_title("Total alarms (SUM)")
        fig.suptitle(
            f"OnlineRFPerm → SAFFRON/ADDIS  first1/2/3 form  (α={args.alpha})",
            fontsize=11,
        )
        fig.tight_layout()
        fig.savefig(args.out_dir / "saffron_addis_compare.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    # markdown report
    md = [
        "# SAFFRON / ADDIS — manuscript first1/2/3 form",
        "",
        "Counting rule (`agod/first_k_metrics.py`, OnlinePermOOB):",
        "```python",
        "first_k = first 1-based end index of k consecutive rejects",
        "P25 / median / P75 = quantiles of *all* alarm times on the trail",
        "DetRate = n_alarm_days / n_days   # NYC-taxi / metro / beijing",
        "```",
        "",
        f"Settings: α={args.alpha}, wealth=α/2, SAFFRON λ={args.saffron_lambda}, "
        f"ADDIS τ={args.addis_tau}, ADDIS λ∈{args.addis_lambdas}, "
        f"burn={args.n_burn}, grace={args.grace}, batch={args.batch_size}, "
        f"seeds={args.seeds}.",
        "",
        f"Datasets: {', '.join(datasets)}.",
        "",
        "## Stream form — first1 / first2 / first3 + alarm P25/median/P75",
        "",
        "| dataset | procedure | first1 | first2 | first3 | P25 | median | P75 | SUM |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    if len(form_agg):
        for _, r in form_agg.iterrows():
            md.append(
                f"| {r['dataset']} | {r['procedure']} | "
                f"{fmt_num(r['first1'])} | {fmt_num(r['first2'])} | {fmt_num(r['first3'])} | "
                f"{fmt_num(r['alarm_p25'])} | {fmt_num(r['alarm_median'])} | "
                f"{fmt_num(r['alarm_p75'])} | {fmt_num(r['SUM'], 1)} |"
            )
    md += ["", "## Day-level DetRate (NYC-taxi / metro / beijing)", ""]
    if len(day_df):
        md += [
            "| dataset | procedure | n_days | n_alarm_days | DetRate | "
            "first1 med | first1 P25 | first1 P75 | first1 mean±std |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for _, r in day_df.iterrows():
            md.append(
                f"| {r['dataset']} | {r['procedure']} | {r['n_days']} | "
                f"{r['n_alarm_days']} | {r['det_rate']:.3f} | "
                f"{fmt_num(r['first1_median'])} | {fmt_num(r['first1_p25'])} | "
                f"{fmt_num(r['first1_p75'])} | "
                f"{fmt_num(r['first1_mean'], 2)}±{fmt_num(r['first1_std'], 2)} |"
            )
        nyc = day_df[day_df["dataset"] == "nyc_taxi"]
        if len(nyc):
            md += [
                "",
                "### NYC-taxi headline",
                "",
            ]
            for _, r in nyc.iterrows():
                md.append(
                    f"- **{r['procedure']}**: **{int(r['n_alarm_days'])}/{int(r['n_days'])}** "
                    f"days alarmed (DetRate={r['det_rate']:.3f})"
                )
    else:
        md.append("_no day-level rows_")

    md += [
        "",
        "## Overall (median across datasets×seeds)",
        "```",
        overall.to_string(index=False) if len(overall) else "(empty)",
        "```",
        "",
        "## Takeaway",
        "- Real-data readout is **first1/2/3 + running alarm-time range (P25/median/P75)**, not overall AR.",
        "- Day packs report **how many calendar days fired** (NYC-taxi DetRate).",
        "- Raising ADDIS λ → more conservative (smaller `(τ−λ)` spend).",
        "",
    ]
    (args.out_dir / "SAFFRON_ADDIS_REPORT.md").write_text("\n".join(md))
    (args.out_dir / "meta.json").write_text(
        json.dumps(
            {
                "alpha": args.alpha,
                "saffron_lambda": args.saffron_lambda,
                "addis_lambdas": args.addis_lambdas,
                "addis_tau": args.addis_tau,
                "datasets": datasets,
                "seeds": args.seeds,
                "form": "first1/2/3 + alarm P25/median/P75 + day DetRate",
            },
            indent=2,
        )
    )
    print("wrote", args.out_dir)
    if len(overall):
        print(overall.to_string(index=False))
    if len(day_df):
        print(day_df.to_string(index=False))


if __name__ == "__main__":
    main()
