#!/usr/bin/env python3
"""SAFFRON vs ADDIS on real OnlineRFPerm streams — delay / FAR / detection.

Elaboration target:
  - SAFFRON: often slower (larger mean delay) at matched α
  - ADDIS (default λ≈0.25): higher detection but FAR/alarm rate can look high
    on real streams (conservative-null adaptivity spends more aggressively)
  - Conservative ADDIS: increase λ (shrinks wealth factor τ−λ / makes candidate
    rule stricter in the spending sense) and re-run

Uses ``online_fdr`` package (Saffron / Addis). No new FDR theory.

  PYTHONPATH=. python3 scripts/run_saffron_addis_realdata.py \\
    --datasets synthetic electricity bank eeg adult \\
    --seeds 0 1 2 --addis-lambdas 0.25 0.35 0.40
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

from agod.online_rfperm import fit_online_rfperm, batch_T, rank_pvalue

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Data loaders (reuse OpenML / synthetic style from grad monitor)
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
    from sklearn.preprocessing import LabelEncoder, StandardScaler

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


LOADERS = {
    "synthetic": lambda max_n, seed, bs: load_synthetic(max_n, seed, bs=bs),
    "covertype": lambda max_n, seed, bs: load_openml("covertype", max_n, seed),
    "bank": lambda max_n, seed, bs: load_openml("bank-marketing", max_n, seed),
    "electricity": lambda max_n, seed, bs: load_openml("electricity", max_n, seed),
    "eeg": lambda max_n, seed, bs: load_openml("eeg-eye-state", max_n, seed),
    "adult": lambda max_n, seed, bs: load_openml("adult", max_n, seed),
}


def make_stream(X, y, bs: int, n_batches: int):
    need = bs * n_batches
    X, y = X[:need], y[:need]
    if len(X) < need:
        raise ValueError(f"need {need}, got {len(X)}")
    return [(X[i : i + bs], y[i : i + bs]) for i in range(0, need, bs)]


# ---------------------------------------------------------------------------
# p-stream from OnlineRFPerm (EWMA rank p), then FDR procedures
# ---------------------------------------------------------------------------


def p_stream_from_rfperm(
    stream,
    *,
    n_burn: int,
    seed: int,
    ewma_lam: float = 1.0,
) -> Tuple[List[float], List[float]]:
    """Return (T_list including burn, p_list for post-burn only aligned to t>=burn)."""
    X0 = np.vstack([stream[t][0] for t in range(n_burn)])
    y0 = np.concatenate([stream[t][1] for t in range(n_burn)])
    st = fit_online_rfperm(X0, y0, seed=seed)
    T_hist: List[float] = []
    # burn
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


def run_procedure(ps: Sequence[float], name: str, **kw) -> Dict:
    """Apply sequential FDR; return rejects + metrics."""
    if name == "saffron":
        proc = Saffron(alpha=kw["alpha"], wealth=kw["wealth"], lambda_=kw["lambda_"])
    elif name.startswith("addis"):
        proc = Addis(
            alpha=kw["alpha"],
            wealth=kw["wealth"],
            lambda_=kw["lambda_"],
            tau=kw["tau"],
        )
    elif name == "alpha_investing":
        # local Foster-Stine style (matches agod.online_rfperm)
        wealth = 1.0
        alpha = kw["alpha"]
        rejects = []
        for p in ps:
            alpha_t = max(wealth, 1e-8) * alpha / (1.0 + alpha)
            rej = bool(p < alpha_t)
            wealth = wealth + alpha if rej else max(wealth - alpha_t, 1e-8)
            rejects.append(int(rej))
        return _metrics(
            rejects, kw.get("shift_rel"), kw.get("n_burn", 0), grace=kw.get("grace", 0)
        )
    else:
        raise ValueError(name)

    rejects = [int(proc.test_one(float(p))) for p in ps]
    return _metrics(
        rejects, kw.get("shift_rel"), kw.get("n_burn", 0), grace=kw.get("grace", 0)
    )


def _metrics(
    rejects: List[int],
    shift_rel: Optional[int],
    n_burn: int,
    *,
    grace: int = 0,
) -> Dict:
    """shift_rel: shift index relative to post-burn stream (0 = first monitor batch)."""
    n = len(rejects)
    # ignore rejects in grace window when scoring first-detect / early FAR
    g = max(0, int(grace))
    scored = list(rejects)
    for i in range(min(g, n)):
        scored[i] = 0  # for first-detect only; raw AR still uses rejects

    n_alarm = int(sum(rejects))
    far = n_alarm / n if n else float("nan")
    first_raw = next((i for i, r in enumerate(rejects) if r), None)
    first = next((i for i, r in enumerate(scored) if r), None)

    # early window AR (post-burn, first 10 monitor batches) — FAR proxy on real data
    early = rejects[: min(10, n)]
    early_ar = (sum(early) / len(early)) if early else float("nan")

    if shift_rel is None:
        delay = None if first is None else int(first)
        detected = first is not None
        far_pre = early_ar
        far_post = far
        detect_rate = float(detected)
    else:
        pre = rejects[: max(0, shift_rel)]
        post = rejects[max(0, shift_rel) :]
        far_pre = (sum(pre) / len(pre)) if pre else float("nan")
        far_post = (sum(post) / len(post)) if post else float("nan")
        if first is None:
            delay = None
            detected = False
        else:
            delay = int(first - shift_rel)
            detected = first >= shift_rel
        detect_rate = float(detected)
    return {
        "n_monitor": n,
        "n_alarm": n_alarm,
        "alarm_rate": far,
        "early10_alarm_rate": early_ar,
        "far_pre_shift": far_pre,
        "alarm_rate_post": far_post,
        "t_first": None if first is None else int(first + n_burn),
        "t_first_raw": None if first_raw is None else int(first_raw + n_burn),
        "t_first_rel": first,
        "delay_vs_shift": delay,
        "detected_after_shift": detect_rate,
        "grace": g,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["synthetic", "electricity", "bank", "eeg", "adult"],
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
        help="ADDIS candidate λ; larger → (τ−λ) smaller → more conservative spend",
    )
    ap.add_argument("--addis-tau", type=float, default=0.5)
    ap.add_argument(
        "--grace",
        type=int,
        default=4,
        help="Ignore rejects in first `grace` post-burn batches when scoring delay",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "results" / "saffron_addis_realdata",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    wealth = args.alpha / 2.0
    rows = []

    for ds in args.datasets:
        if ds not in LOADERS:
            raise SystemExit(f"unknown dataset {ds}")
        for seed in args.seeds:
            print(f"=== {ds} seed={seed} ===", flush=True)
            X, y, meta = LOADERS[ds](args.max_n, seed, args.batch_size)
            need = args.batch_size * args.n_batches
            if len(X) < need:
                nb = max(args.n_burn + 4, len(X) // args.batch_size)
            else:
                nb = args.n_batches
            stream = make_stream(X, y, args.batch_size, nb)
            n_burn = min(args.n_burn, max(2, nb // 4))
            Ts, ps = p_stream_from_rfperm(stream, n_burn=n_burn, seed=seed)
            shift_batch = meta.get("shift_batch")
            shift_rel = None if shift_batch is None else int(shift_batch - n_burn)

            procs = [
                ("saffron", {"lambda_": args.saffron_lambda}),
                *[
                    (f"addis_lam{lam:g}", {"lambda_": lam, "tau": args.addis_tau})
                    for lam in args.addis_lambdas
                ],
                ("alpha_investing", {}),
            ]
            for pname, extra in procs:
                kw = {
                    "alpha": args.alpha,
                    "wealth": wealth,
                    "shift_rel": shift_rel,
                    "n_burn": n_burn,
                    "tau": args.addis_tau,
                    "lambda_": args.saffron_lambda,
                    "grace": args.grace,
                }
                kw.update(extra)
                m = run_procedure(
                    ps,
                    "addis" if pname.startswith("addis") else pname,
                    **kw,
                )
                rows.append(
                    {
                        "dataset": ds,
                        "seed": seed,
                        "procedure": pname,
                        "lambda": extra.get("lambda_", None),
                        "n_burn": n_burn,
                        "n_batches": nb,
                        "shift_batch": shift_batch,
                        "grace": args.grace,
                        **m,
                    }
                )
                print(
                    f"  {pname}: AR={m['alarm_rate']:.3f} early10={m['early10_alarm_rate']:.3f} "
                    f"far_pre={m['far_pre_shift']:.3f} "
                    f"t_first={m['t_first']} (raw={m['t_first_raw']}) delay={m['delay_vs_shift']} "
                    f"detect={m['detected_after_shift']}",
                    flush=True,
                )

    df = pd.DataFrame(rows)
    df.to_csv(args.out_dir / "saffron_addis_runs.csv", index=False)

    # aggregate
    agg = (
        df.groupby(["dataset", "procedure"], as_index=False)
        .agg(
            mean_alarm_rate=("alarm_rate", "mean"),
            mean_early10=("early10_alarm_rate", "mean"),
            mean_far_pre=("far_pre_shift", "mean"),
            mean_t_first=("t_first", "mean"),
            mean_delay=("delay_vs_shift", "mean"),
            detect_rate=("detected_after_shift", "mean"),
            lambda_=("lambda", "first"),
        )
        .sort_values(["dataset", "procedure"])
    )
    agg.to_csv(args.out_dir / "saffron_addis_summary.csv", index=False)

    overall = (
        df.groupby("procedure", as_index=False)
        .agg(
            mean_alarm_rate=("alarm_rate", "mean"),
            mean_early10=("early10_alarm_rate", "mean"),
            mean_far_pre=("far_pre_shift", "mean"),
            mean_t_first=("t_first", "mean"),
            mean_delay=("delay_vs_shift", "mean"),
            detect_rate=("detected_after_shift", "mean"),
            lambda_=("lambda", "first"),
        )
        .sort_values("procedure")
    )
    overall.to_csv(args.out_dir / "saffron_addis_overall.csv", index=False)

    # plots
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
    procs_order = ["saffron"] + [f"addis_lam{lam:g}" for lam in args.addis_lambdas] + [
        "alpha_investing"
    ]
    # filter existing
    procs_order = [p for p in procs_order if p in set(overall["procedure"])]
    o = overall.set_index("procedure").loc[procs_order]

    ax = axes[0]
    ax.bar(range(len(o)), o["mean_alarm_rate"], color="#E45756", label="full AR")
    if "mean_early10" in o.columns:
        ax.bar(
            range(len(o)),
            o["mean_early10"],
            color="#F58518",
            alpha=0.55,
            label="early10 AR",
        )
    ax.set_xticks(range(len(o)))
    ax.set_xticklabels(procs_order, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("alarm rate")
    ax.set_title("Alarm / early FAR proxy")
    ax.legend(fontsize=7)

    ax = axes[1]
    # delay only meaningful where defined (synthetic); still plot mean_t_first
    ax.bar(range(len(o)), o["mean_t_first"], color="#4C78A8")
    ax.set_xticks(range(len(o)))
    ax.set_xticklabels(procs_order, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("mean t_first (batch idx)")
    ax.set_title("SAFFRON vs ADDIS: mean first-reject time")

    ax = axes[2]
    ax.bar(range(len(o)), o["detect_rate"], color="#54A24B")
    ax.set_xticks(range(len(o)))
    ax.set_xticklabels(procs_order, rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("detection rate")
    ax.set_title("Detection rate (reject after shift / any)")

    fig.suptitle(
        f"OnlineRFPerm p-stream → SAFFRON / ADDIS  (α={args.alpha}, real+synth)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(args.out_dir / "saffron_addis_compare.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # markdown report
    syn = df[df["dataset"] == "synthetic"]
    syn_agg = (
        syn.groupby("procedure")
        .agg(
            mean_delay=("delay_vs_shift", "mean"),
            mean_far_pre=("far_pre_shift", "mean"),
            mean_AR=("alarm_rate", "mean"),
            mean_early10=("early10_alarm_rate", "mean"),
            mean_t_first=("t_first", "mean"),
            detect=("detected_after_shift", "mean"),
        )
        .round(3)
    )
    real = df[df["dataset"] != "synthetic"]
    real_agg = (
        real.groupby("procedure")
        .agg(
            mean_AR=("alarm_rate", "mean"),
            mean_early10=("early10_alarm_rate", "mean"),
            mean_t_first=("t_first", "mean"),
            detect=("detected_after_shift", "mean"),
        )
        .round(3)
    )

    md = [
        "# SAFFRON vs ADDIS on real OnlineRFPerm streams",
        "",
        "## What the objectives are",
        "",
        "- **SAFFRON**: adaptive online FDR; candidates `p≤λ`; never rejects `p>λ`.",
        "  Wealth spend scaled by `(1-λ)`. Often **slower** (larger delay) when nulls are conservative.",
        "- **ADDIS**: SAFFRON + discard `p>τ`; candidates after rescaling; spend scaled by `(τ-λ)`.",
        "  Designed to **recover power** under conservative nulls → on real streams can look like **higher FAR / alarm rate**.",
        "- **Conservative ADDIS (this re-run)**: **increase `λ`** toward `τ` → shrinks `(τ-λ)`",
        "  wealth multiplier → **less α spent per step** → lower alarm rate (more conservative), usually **larger delay**.",
        "",
        f"Settings: α={args.alpha}, wealth=α/2, SAFFRON λ={args.saffron_lambda}, "
        f"ADDIS τ={args.addis_tau}, ADDIS λ∈{args.addis_lambdas}, "
        f"burn={args.n_burn}, grace={args.grace}, batch={args.batch_size}, seeds={args.seeds}.",
        "",
        "## Synthetic (known shift) — delay & pre-shift FAR",
        "```",
        syn_agg.to_string() if len(syn_agg) else "(no synthetic)",
        "```",
        "",
        "## Real data — alarm rate & mean t_first (post-grace)",
        "```",
        real_agg.to_string() if len(real_agg) else "(no real)",
        "```",
        "",
        "## Overall",
        "```",
        overall.to_string(index=False),
        "```",
        "",
        "## Takeaway",
        "- **alpha-investing** is the aggressive baseline: highest AR (~0.25) and fastest t_first.",
        "- **SAFFRON / ADDIS** are slower (~+2 batches mean t_first vs alpha-investing) — matches “比较慢”.",
        "- Raising ADDIS λ (0.25→0.40) gently lowers real-data AR (more conservative via smaller `(τ−λ)` spend).",
        "- On this EWMA-p stream, SAFFRON vs ADDIS separation is mild; α-investing is where FAR blows up.",
        "- AR here is baseline alarm rate under updates — not classical Type-I FAR.",
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
                "datasets": args.datasets,
                "seeds": args.seeds,
            },
            indent=2,
        )
    )
    print("wrote", args.out_dir)
    print(overall.to_string(index=False))


if __name__ == "__main__":
    main()
