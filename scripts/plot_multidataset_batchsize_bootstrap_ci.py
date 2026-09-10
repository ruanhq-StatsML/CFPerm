#!/usr/bin/env python3
"""Per-dataset time-index dashboard · null-adjusted online-bootstrap CI.

Primary output: one figure per dataset with observation/time index on x,
batch sizes stacked as rows. Dashboard only (not a methodological claim).

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
        next(reader)
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
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    X[t_change:, :3] += 2.0
    beta = rng.normal(size=p)
    Y = X @ beta + rng.normal(size=n)
    return X, Y, t_change


def diabetes_xy(*, n_ref: int = N_REF, n_pre: int = 500, n_post: int = 4500, seed: int = SEED):
    Xs, Ys = load_csv_xy(DATA / "source_DiabetesReadmission.csv")
    Xt, Yt = load_csv_xy(DATA / "target_DiabetesReadmission.csv")
    rng = np.random.default_rng(seed)
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
    return X, Y, n_need_src


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
    else:
        X, Y, t_change = maker(n=N_REF + N_TRAIL, seed=SEED + (hash(ds_key) % 10007))
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

    if t_change is None:
        cover0_pre = cover0
        det = float("nan")
        delay = float("nan")
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


def _legend_unique(ax):
    handles, labels = ax.get_legend_handles_labels()
    seen, h2, l2 = set(), [], []
    for h, lab in zip(handles, labels):
        if lab in seen:
            continue
        seen.add(lab)
        h2.append(h)
        l2.append(lab)
    return h2, l2


def _draw_time_panel(ax, rows, summary, *, batch: int, show_xlabel: bool = False, ylabel: str | None = None):
    x = np.asarray([r["obs"] for r in rows], float)
    s = np.asarray([r["s"] for r in rows], float)
    m = np.asarray([r["mean"] for r in rows], float)
    lo = np.asarray([r["lo"] for r in rows], float)
    hi = np.asarray([r["hi"] for r in rows], float)
    x0, x1 = float(np.nanmin(x)), float(np.nanmax(x))
    tc = summary["t_change"]

    if tc is not None and x0 < tc < x1:
        ax.axvspan(x0, tc, color="#e8eef5", alpha=0.95, zorder=0, label="pre-CP window")
        ax.axvspan(tc, x1, color="#f7efe6", alpha=0.95, zorder=0, label="post-CP window")
        ax.axvline(tc, color="#d97706", ls=":", lw=1.8, label="labeled CP", zorder=3)
    else:
        ax.axvspan(x0, x1, color="#e8eef5", alpha=0.8, zorder=0)

    ax.fill_between(x, lo, hi, color="#9ecae1", alpha=0.7, label="online AR-bootstrap CI", zorder=2)
    ax.plot(
        x,
        s,
        color="#1f4e79",
        lw=0.5 if batch == 1 else 1.0,
        alpha=0.22 if batch == 1 else 0.45,
        label=r"batch $\Delta_t$",
        zorder=3,
    )
    ax.plot(x, m, color="#2f6f4e", lw=2.15, label=r"running mean of $\Delta$", zorder=4)
    ax.axhline(0.0, color="#a33b3b", lw=1.15, alpha=0.9, label="null 0", zorder=3)

    if tc is not None:
        for r in rows:
            if r["obs"] >= tc and r["lo"] > 0.0 and r["t"] > BURN:
                ax.axvline(
                    r["obs"],
                    color="#2f6f4e",
                    ls="-.",
                    lw=1.4,
                    label=r"first CI$_{\mathrm{lo}}>0$",
                    zorder=4,
                )
                break

    bits = [
        f"B={batch}",
        rf"$\mu_{{\mathrm{{ref}}}}={summary['mu_ref']:.3f}$",
        f"TV(lo)/1k={summary['tv_lo_per_1k']:.2f}",
    ]
    if tc is None:
        bits.append(f"cover0={summary['cover0']:.2f}")
    else:
        bits.append(f"cover0_pre={summary['cover0_pre']:.2f}")
        if np.isfinite(summary["post_det"]):
            bits.append(f"post-det={summary['post_det']:.2f}")
        if np.isfinite(summary["delay"]):
            bits.append(f"delay={summary['delay']:.0f} batches")
    ax.text(
        0.01,
        0.96,
        "  ·  ".join(bits),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.6,
        color="#222",
    )
    ax.set_ylabel(ylabel if ylabel is not None else rf"$B={batch}$" + "\n" + r"$\Delta$", fontsize=10)
    if show_xlabel:
        ax.set_xlabel("observation / time index", fontsize=11)
    ax.grid(True, alpha=0.28)
    y_core = np.concatenate(
        [m[np.isfinite(m)], lo[np.isfinite(lo)], hi[np.isfinite(hi)], np.array([0.0])]
    )
    if y_core.size:
        lo_y, hi_y = np.nanpercentile(y_core, [2, 98])
        pad = 0.18 * max(hi_y - lo_y, 1e-3)
        ax.set_ylim(lo_y - pad, hi_y + pad)
    ax.set_xlim(x0, x1)


def plot_dataset_timeindex(ds_key: str, results: dict):
    meta = DATASETS[ds_key]
    fig, axes = plt.subplots(len(BATCHES), 1, figsize=(12.2, 9.4), sharex=True)
    for i, (ax, batch) in enumerate(zip(axes, BATCHES)):
        rows, summary = results[batch]
        _draw_time_panel(ax, rows, summary, batch=batch, show_xlabel=(i == len(BATCHES) - 1))
    h2, l2 = _legend_unique(axes[-1])
    fig.legend(h2, l2, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.018))
    fig.suptitle(
        rf"{meta['title']} · time-index path of $\Delta_t=s_t-\mu_{{\mathrm{{ref}}}}$"
        f"\nOnline AR-bootstrap CI · n_boot={N_BOOT} · rows = batch size",
        fontsize=12.5,
        y=1.05,
    )
    fig.tight_layout()
    dest = OUT / f"{ds_key}_timeindex_bootstrap_ci.png"
    fig.savefig(dest, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}", flush=True)


def plot_gallery_b20(all_results: dict):
    focus_b = 20
    keys = list(DATASETS.keys())
    fig, axes = plt.subplots(len(keys), 1, figsize=(12.4, 10.2), sharex=False)
    for i, (ax, ds) in enumerate(zip(axes, keys)):
        rows, summary = all_results[ds][focus_b]
        _draw_time_panel(
            ax,
            rows,
            summary,
            batch=focus_b,
            show_xlabel=(i == len(keys) - 1),
            ylabel=DATASETS[ds]["title"] + "\n" + r"$\Delta$",
        )
    h2, l2 = _legend_unique(axes[0])
    fig.legend(h2, l2, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.015))
    fig.suptitle(
        rf"Per-dataset time-index gallery · $B={focus_b}$ · $\Delta_t=s_t-\mu_{{\mathrm{{ref}}}}$",
        fontsize=12.5,
        y=1.04,
    )
    fig.tight_layout()
    dest = OUT / "multidataset_timeindex_gallery_b20.png"
    fig.savefig(dest, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}", flush=True)


def write_tex(all_summaries: list[dict]):
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Per-dataset time-index board (dashboard). "
        r"Null-adjusted $\Delta_t=s_t-\mu_{\mathrm{ref}}$ with online AR-bootstrap CI along observation index. "
        r"cover0 / cover0 pre $=$ fraction of CIs containing $0$. "
        r"post det $=$ post-CP fraction with CI$_{\mathrm{lo}}>0$. "
        r"delay $=$ batches from labeled CP to first significant CI.}",
        r"\label{tab:multidataset-timeindex-boot-ci}",
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
            det_s = f"{s['post_det']:.3f}" if np.isfinite(s["post_det"]) else "--"
            delay_s = f"{s['delay']:.0f}" if np.isfinite(s["delay"]) else "--"
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
    all_results: dict = {}
    for ds in DATASETS:
        print(f"running {ds} time-index CI grid", flush=True)
        results = {b: run_one(ds, b) for b in BATCHES}
        all_results[ds] = results
        plot_dataset_timeindex(ds, results)
        for b in BATCHES:
            s = results[b][1]
            all_summaries.append(s)
            cover = s["cover0_pre"] if DATASETS[ds]["has_cp"] else s["cover0"]
            print(
                f"  {ds} B={b}: TV/1k={s['tv_lo_per_1k']:.3f} "
                f"cover0/pre={cover:.3f} det={s['post_det']} delay={s['delay']}",
                flush=True,
            )

    plot_gallery_b20(all_results)
    write_tex(all_summaries)
    print("done", flush=True)


if __name__ == "__main__":
    main()
