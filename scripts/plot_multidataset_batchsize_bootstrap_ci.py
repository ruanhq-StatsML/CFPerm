#!/usr/bin/env python3
"""Multi-dataset dashboard · null-adjusted online-bootstrap CI vs batch size.

Board only (not a methodological claim). Same protocol on:
  stationary · concept · covariate-shift · DiabetesReadmission (source→target)

  python3 scripts/plot_multidataset_batchsize_bootstrap_ci.py
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "streaming_porisk_online_bootstrap"
DOCS = ROOT / "docs" / "method"
DATA = ROOT / "datasets" / "extracted"

SEED = 2026
N_REF = 1000
N_TRAIL = 5000
T_CHANGE_OBS = N_REF + 1500
N_BOOT = 200
ALPHA = 0.05
BETA = float(2.0 ** 0.5 - 1.0)
BATCHES = (1, 10, 20, 50)
BURN = 5


class OnlineARBootstrap:
    def __init__(self, n_boot: int = 200, seed: int = 0):
        self.n_boot = int(n_boot)
        self._rng = np.random.default_rng(seed)
        self.t = 0
        self.V = np.zeros(self.n_boot)
        self.Vbar = np.zeros(self.n_boot)
        self.Xstar = np.zeros(self.n_boot)
        self.mean = np.nan

    def update(self, x: float) -> None:
        t = self.t + 1
        rho = float(np.clip(1.0 - t ** (-BETA), 0.0, 1.0 - 1e-12))
        zeta = self._rng.normal(size=self.n_boot)
        self.V = 1.0 + rho * (self.V - 1.0) + np.sqrt(max(0.0, 1.0 - rho * rho)) * zeta
        if t == 1:
            self.Xstar = np.full(self.n_boot, float(x))
            self.Vbar = self.V.copy()
            self.mean = float(x)
        else:
            num = (t - 1) * self.Vbar * self.Xstar + float(x) * self.V
            den = (t - 1) * self.Vbar + self.V
            ok = np.abs(den) > 1e-15
            nxt = np.full(self.n_boot, float(x))
            nxt[ok] = num[ok] / den[ok]
            self.Xstar = nxt
            self.Vbar = (1.0 - 1.0 / t) * self.Vbar + self.V / t
            self.mean = ((t - 1) * float(self.mean) + float(x)) / t
        self.t = t

    def ci(self, alpha: float = 0.05) -> tuple[float, float]:
        if self.t < 2:
            return float("nan"), float("nan")
        lo, hi = np.quantile(self.Xstar, [alpha / 2.0, 1.0 - alpha / 2.0])
        return float(lo), float(hi)


def load_csv_xy(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open() as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [[float(v) for v in row] for row in reader]
    arr = np.asarray(rows, dtype=float)
    return arr[:, :-1], arr[:, -1]


def stationary_xy(*, n: int, seed: int, p: int = 8):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    beta = rng.normal(size=p)
    Y = X @ beta + rng.normal(size=n)
    return X, Y, None


def concept_xy(*, n: int, seed: int, p: int = 8, t_change: int = T_CHANGE_OBS):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    beta = rng.normal(size=p)
    signal = X @ beta
    Y = signal + rng.normal(size=n)
    Y[t_change:] = signal[t_change:] + 2.5 * X[t_change:, 0] + rng.normal(size=n - t_change)
    return X, Y, t_change


def covariate_xy(*, n: int, seed: int, p: int = 8, t_change: int = T_CHANGE_OBS):
    """P(X) mean-shift after t_change; P(Y|X) unchanged."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    X[t_change:, :3] += 2.0
    beta = rng.normal(size=p)
    Y = X @ beta + rng.normal(size=n)
    return X, Y, t_change


def diabetes_xy(*, n_ref: int = N_REF, n_pre: int = 500, n_post: int = 4500, seed: int = SEED):
    """Source → target domain stream. CP at end of source pre-trail."""
    Xs, Ys = load_csv_xy(DATA / "source_DiabetesReadmission.csv")
    Xt, Yt = load_csv_xy(DATA / "target_DiabetesReadmission.csv")
    rng = np.random.default_rng(seed)
    # shuffle within domains so the board isn't order artifacts
    is_ = rng.permutation(len(Ys))
    it_ = rng.permutation(len(Yt))
    Xs, Ys = Xs[is_], Ys[is_]
    Xt, Yt = Xt[it_], Yt[it_]
    n_need_src = n_ref + n_pre
    if len(Ys) < n_need_src:
        raise RuntimeError(f"source too short: {len(Ys)} < {n_need_src}")
    n_post = min(n_post, len(Yt))
    X = np.vstack([Xs[:n_need_src], Xt[:n_post]])
    Y = np.concatenate([Ys[:n_need_src], Yt[:n_post]])
    t_change = n_need_src
    return X, Y, t_change


DATASETS = {
    "stationary": dict(title="Stationary", maker=stationary_xy, has_cp=False),
    "concept": dict(title="Concept drift", maker=concept_xy, has_cp=True),
    "covariate": dict(title="Covariate shift", maker=covariate_xy, has_cp=True),
    "diabetes": dict(title="Diabetes (src→tgt)", maker=diabetes_xy, has_cp=True),
}


def freeze_and_score(X, Y, *, n_ref: int, batch: int):
    X0, Y0 = X[:n_ref], Y[:n_ref]
    beta_hat = np.linalg.lstsq(X0, Y0, rcond=None)[0]

    def one_score(s: int, e: int) -> float:
        resid = Y[s:e] - X[s:e] @ beta_hat
        n1 = e - s
        ee = float(np.clip(n1 / max(n_ref + n1, 1), 0.02, 0.98))
        return float(np.mean((resid * (1.0 - ee)) ** 2))

    ref_scores = []
    for s in range(0, n_ref, batch):
        e = min(n_ref, s + batch)
        if e > s:
            ref_scores.append(one_score(s, e))
    if not ref_scores:
        ref_scores = [one_score(0, n_ref)]
    mu_ref = float(np.mean(ref_scores))

    scores, times = [], []
    t = 0
    for s in range(n_ref, len(Y), batch):
        e = min(len(Y), s + batch)
        if e <= s:
            continue
        t += 1
        scores.append(one_score(s, e))
        times.append(t)
    return np.asarray(times, int), np.asarray(scores, float), mu_ref


def run_one(ds_key: str, batch: int):
    meta = DATASETS[ds_key]
    maker = meta["maker"]
    if ds_key == "diabetes":
        X, Y, t_change = maker(seed=SEED + 17)
        n_ref = N_REF
    else:
        X, Y, t_change = maker(n=N_REF + N_TRAIL, seed=SEED + hash(ds_key) % 10007)
        n_ref = N_REF
    times, scores, mu_ref = freeze_and_score(X, Y, n_ref=n_ref, batch=batch)
    delta = scores - mu_ref
    boot = OnlineARBootstrap(n_boot=N_BOOT, seed=SEED + 100 * batch + (hash(ds_key) % 97))
    rows = []
    for t, d in zip(times, delta):
        boot.update(float(d))
        lo, hi = boot.ci(ALPHA)
        rows.append(
            dict(
                t=int(t),
                obs=int(n_ref + t * batch),
                s=float(d),
                mean=float(boot.mean),
                lo=float(lo),
                hi=float(hi),
            )
        )
    post = [r for r in rows if r["t"] > BURN and np.isfinite(r["lo"])]
    los = np.asarray([r["lo"] for r in post], float)
    his = np.asarray([r["hi"] for r in post], float)
    ss = np.asarray([r["s"] for r in post], float)
    dlo = np.diff(los) if len(los) > 1 else np.asarray([np.nan])
    n_obs = float(max(len(post) * batch, 1))
    tv_lo_per_1k = float(np.nansum(np.abs(dlo)) / n_obs * 1000.0)
    cover0 = float(np.mean((los <= 0.0) & (his >= 0.0))) if post else float("nan")

    # pre / post relative to labeled CP (if any)
    if t_change is None:
        pre = post
        post_cp = []
        delay = float("nan")
        det = float("nan")
        cover0_pre = cover0
    else:
        pre = [r for r in post if r["obs"] < t_change]
        post_cp = [r for r in post if r["obs"] >= t_change]
        cover0_pre = (
            float(np.mean([(r["lo"] <= 0.0) and (r["hi"] >= 0.0) for r in pre]))
            if pre
            else float("nan")
        )
        det = float(np.mean([r["lo"] > 0.0 for r in post_cp])) if post_cp else float("nan")
        delay = float("nan")
        for r in post_cp:
            if r["lo"] > 0.0:
                # batches after CP
                delay = float((r["obs"] - t_change) / batch)
                break

    summary = dict(
        dataset=ds_key,
        batch=batch,
        mu_ref=float(mu_ref),
        tv_lo_per_1k=tv_lo_per_1k,
        std_s=float(np.nanstd(ss)) if len(ss) else float("nan"),
        cover0=cover0,
        cover0_pre=cover0_pre,
        post_det=det,
        delay=delay,
        t_change=t_change,
        n_ref=n_ref,
        n_total=int(len(Y)),
    )
    return rows, summary


def plot_dataset_paths(ds_key: str, results: dict):
    meta = DATASETS[ds_key]
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.0), sharey=False)
    for ax, batch in zip(axes.ravel(), BATCHES):
        rows, summary = results[batch]
        x = np.asarray([r["obs"] for r in rows], float)
        s = np.asarray([r["s"] for r in rows], float)
        m = np.asarray([r["mean"] for r in rows], float)
        lo = np.asarray([r["lo"] for r in rows], float)
        hi = np.asarray([r["hi"] for r in rows], float)
        ax.fill_between(x, lo, hi, color="#9ecae1", alpha=0.55, label="online AR-bootstrap CI on Δ")
        ax.plot(x, m, color="#2f6f4e", lw=1.8, label="running mean of Δ")
        ax.plot(
            x,
            s,
            color="#1f4e79",
            lw=0.6 if batch == 1 else 1.2,
            alpha=0.25 if batch == 1 else 0.55,
            label=r"batch excess $\Delta$",
        )
        ax.axhline(0.0, color="#a33b3b", lw=1.0, alpha=0.8, label="null 0")
        tc = summary["t_change"]
        if tc is not None:
            ax.axvline(tc, color="#d97706", ls=":", lw=1.4, label="labeled CP")
            # first significant post-CP
            for r in rows:
                if r["obs"] >= tc and r["lo"] > 0.0 and r["t"] > BURN:
                    ax.axvline(r["obs"], color="#2f6f4e", ls="-.", lw=1.2, label="first CI_lo>0")
                    break
        title = (
            f"{meta['title']} · B={batch} · μ_ref={summary['mu_ref']:.3f}\n"
            f"TV(lo)/1k={summary['tv_lo_per_1k']:.3f}"
        )
        if meta["has_cp"]:
            title += f" · cover0_pre={summary['cover0_pre']:.2f} · det={summary['post_det']:.2f}"
            if np.isfinite(summary["delay"]):
                title += f" · delay={summary['delay']:.0f}"
        else:
            title += f" · cover0={summary['cover0']:.2f}"
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("observation index")
        ax.grid(True, alpha=0.3)
        y_core = np.concatenate([m, lo, hi, np.array([0.0])])
        y_core = y_core[np.isfinite(y_core)]
        if y_core.size:
            lo_y, hi_y = np.nanpercentile(y_core, [1, 99])
            pad = 0.2 * max(hi_y - lo_y, 1e-3)
            ax.set_ylim(lo_y - pad, hi_y + pad)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    # de-dup legend
    seen = set()
    h2, l2 = [], []
    for h, lab in zip(handles, labels):
        if lab in seen:
            continue
        seen.add(lab)
        h2.append(h)
        l2.append(lab)
    fig.legend(h2, l2, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle(
        rf"{meta['title']} · null-adjusted online-bootstrap CI · $\Delta_t=s_t-\mu_{{\mathrm{{ref}}}}$"
        f"\nPalm & Nagler AR-bootstrap · n_boot={N_BOOT} · dashboard only",
        fontsize=12,
        y=1.08,
    )
    axes[0, 0].set_ylabel("excess Δ")
    axes[1, 0].set_ylabel("excess Δ")
    fig.tight_layout()
    dest = OUT / f"{ds_key}_batchsize_bootstrap_ci.png"
    fig.savefig(dest, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}", flush=True)


def plot_overview(all_summaries: list[dict]):
    ds_order = list(DATASETS.keys())
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))
    x = np.arange(len(BATCHES))
    width = 0.18
    colors = ["#1f4e79", "#c47b17", "#2f6f4e", "#7a3e9d"]

    for i, ds in enumerate(ds_order):
        rows = [s for s in all_summaries if s["dataset"] == ds]
        rows = sorted(rows, key=lambda r: r["batch"])
        axes[0].bar(x + (i - 1.5) * width, [r["tv_lo_per_1k"] for r in rows], width, color=colors[i], label=DATASETS[ds]["title"])
        cover = [r["cover0_pre"] if DATASETS[ds]["has_cp"] else r["cover0"] for r in rows]
        axes[1].bar(x + (i - 1.5) * width, cover, width, color=colors[i])
        det = [r["post_det"] if np.isfinite(r["post_det"]) else 0.0 for r in rows]
        axes[2].bar(x + (i - 1.5) * width, det, width, color=colors[i])

    axes[0].set_title(r"CI twitch · TV(lo)/1k")
    axes[1].set_title(r"cover0 (pre-CP / stationary)")
    axes[2].set_title(r"post-CP det · frac CI$_{lo}>0$")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels([str(b) for b in BATCHES])
        ax.set_xlabel("batch size")
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle(
        r"Multi-dataset board · null-adjusted $\Delta$ · batch-size trend (dashboard)",
        fontsize=12,
    )
    fig.tight_layout()
    dest = OUT / "multidataset_batchsize_bootstrap_ci_overview.png"
    fig.savefig(dest, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}", flush=True)


def write_tex(all_summaries: list[dict]):
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Multi-dataset dashboard (not a methodological claim). "
        r"Null-adjusted $\Delta_t=s_t-\mu_{\mathrm{ref}}$ with online AR-bootstrap CI. "
        r"cover0 / cover0 pre $=$ fraction of CIs containing $0$ (stationary / pre-CP). "
        r"post det $=$ post-CP fraction with CI$_{\mathrm{lo}}>0$. "
        r"delay $=$ batches from labeled CP to first significant CI.}",
        r"\label{tab:multidataset-batchsize-boot-ci}",
        r"\small",
        r"\begin{tabular}{@{}l r ccccc@{}}",
        r"\toprule",
        r"dataset & $B$ & $\mu_{\mathrm{ref}}$ & TV$(\mathrm{lo})/1\mathrm{k}$ & cover0/pre & post det & delay \\",
        r"\midrule",
    ]
    for ds in DATASETS:
        for b in BATCHES:
            s = next(x for x in all_summaries if x["dataset"] == ds and x["batch"] == b)
            cover = s["cover0_pre"] if DATASETS[ds]["has_cp"] else s["cover0"]
            det = s["post_det"] if np.isfinite(s["post_det"]) else float("nan")
            delay = s["delay"] if np.isfinite(s["delay"]) else float("nan")
            det_s = f"{det:.3f}" if np.isfinite(det) else "--"
            delay_s = f"{delay:.0f}" if np.isfinite(delay) else "--"
            lines.append(
                f"{DATASETS[ds]['title']} & {b} & {s['mu_ref']:.3f} & {s['tv_lo_per_1k']:.3f} & "
                f"{cover:.3f} & {det_s} & {delay_s} \\\\"
            )
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    tex = DOCS / "Streaming_PORisk_MultiDataset_BatchSize_BootstrapCI_tables_only.tex"
    tex.write_text("\n".join(lines))
    print(f"wrote {tex}", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    if not (DATA / "source_DiabetesReadmission.csv").exists():
        raise SystemExit(f"missing diabetes CSVs under {DATA}")

    all_summaries: list[dict] = []
    for ds in DATASETS:
        print(f"running {ds} batch-size CI grid", flush=True)
        results = {b: run_one(ds, b) for b in BATCHES}
        plot_dataset_paths(ds, results)
        for b in BATCHES:
            s = results[b][1]
            all_summaries.append(s)
            cover = s["cover0_pre"] if DATASETS[ds]["has_cp"] else s["cover0"]
            print(
                f"  {ds} B={b}: TV/1k={s['tv_lo_per_1k']:.3f} "
                f"cover0/pre={cover:.3f} det={s['post_det']} delay={s['delay']}",
                flush=True,
            )

    plot_overview(all_summaries)
    write_tex(all_summaries)
    print("done", flush=True)


if __name__ == "__main__":
    main()
