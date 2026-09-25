#!/usr/bin/env python3
"""Stationary DGP · null-adjusted online-bootstrap CI vs batch size.

Null adjustment: μ_ref = mean frozen PO-risk on D_ref mini-batches.
Trail excess Δ_t = s_t − μ_ref. Bootstrap the centered stream; stationary
CIs for the running mean of Δ should cover 0.

Claim: too-small batches still make the hanging CI twitchy after centering.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "streaming_porisk_online_bootstrap"
DOCS = ROOT / "docs" / "method"

SEED = 2026
N_REF = 1000
N_TRAIL = 5000
N_BOOT = 200
ALPHA = 0.05
BETA = float(2.0 ** 0.5 - 1.0)
BATCHES = (1, 10, 20, 50)
BURN = 5


class OnlineARBootstrap:
    """Palm & Nagler Algorithm 1, vectorized B chains."""

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


def stationary_stream(*, n: int, seed: int, p: int = 8):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    beta = rng.normal(size=p)
    Y = X @ beta + rng.normal(size=n)
    return X, Y


def freeze_and_score(X, Y, *, n_ref: int, batch: int):
    """Freeze OLS on original batch; score ref and trail mini-batches the same way.

    Returns trail (times, scores), plus mu_ref = mean PO-risk on D_ref mini-batches
    (null adjustment floor).
    """
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
    # If batch > n_ref somehow, fall back to whole-ref score
    if not ref_scores:
        ref_scores = [one_score(0, n_ref)]
    mu_ref = float(np.mean(ref_scores))

    scores, times, n_pts = [], [], []
    t = 0
    for s in range(n_ref, len(Y), batch):
        e = min(len(Y), s + batch)
        if e <= s:
            continue
        t += 1
        scores.append(one_score(s, e))
        times.append(t)
        n_pts.append(e - s)
    return (
        np.asarray(times, int),
        np.asarray(scores, float),
        np.asarray(n_pts, int),
        mu_ref,
        np.asarray(ref_scores, float),
    )


def run_batch(batch: int):
    X, Y = stationary_stream(n=N_REF + N_TRAIL, seed=SEED)
    times, scores, n_pts, mu_ref, ref_scores = freeze_and_score(
        X, Y, n_ref=N_REF, batch=batch
    )
    # Null adjustment: excess vs D_ref average PO-risk. Bootstrap the centered stream.
    delta = scores - mu_ref
    boot = OnlineARBootstrap(n_boot=N_BOOT, seed=SEED + batch)
    rows = []
    for t, d in zip(times, delta):
        boot.update(float(d))
        lo, hi = boot.ci(ALPHA)
        rows.append(
            dict(
                t=int(t),
                s=float(d),  # centered excess
                mean=float(boot.mean),
                lo=float(lo),
                hi=float(hi),
                width=float(hi - lo),
                n_pts=int(batch),
                mu_ref=float(mu_ref),
            )
        )
    post = [r for r in rows if r["t"] > BURN and np.isfinite(r["width"])]
    widths = np.asarray([r["width"] for r in post], float)
    los = np.asarray([r["lo"] for r in post], float)
    his = np.asarray([r["hi"] for r in post], float)
    means = np.asarray([r["mean"] for r in post], float)
    ss = np.asarray([r["s"] for r in post], float)
    dlo = np.diff(los) if len(los) > 1 else np.asarray([np.nan])
    n_obs = float(max(len(post) * batch, 1))
    tv_lo_per_1k = float(np.nansum(np.abs(dlo)) / n_obs * 1000.0)
    # Stationary: after null adjustment, CI for running mean should cover 0
    cover0 = float(np.mean((los <= 0.0) & (his >= 0.0))) if post else float("nan")
    summary = dict(
        batch=batch,
        n_batches=len(rows),
        mu_ref=float(mu_ref),
        mean_width=float(np.nanmean(widths)),
        std_s=float(np.nanstd(ss)),
        tv_lo_per_1k=tv_lo_per_1k,
        mean_abs_dlo_batch=float(np.nanmean(np.abs(dlo))),
        cover0=cover0,
        cover_mean=float(np.mean((means >= los) & (means <= his))),
    )
    return rows, summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    print("running stationary batch-size CI grid", flush=True)
    results = {b: run_batch(b) for b in BATCHES}

    # Path board: CI + running mean (the estimand) + light s_b
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.0), sharey=False)
    for ax, batch in zip(axes.ravel(), BATCHES):
        rows, summary = results[batch]
        t = np.asarray([r["t"] for r in rows], float)
        s = np.asarray([r["s"] for r in rows], float)
        m = np.asarray([r["mean"] for r in rows], float)
        lo = np.asarray([r["lo"] for r in rows], float)
        hi = np.asarray([r["hi"] for r in rows], float)
        x_obs = N_REF + t * batch
        ax.fill_between(x_obs, lo, hi, color="#9ecae1", alpha=0.55, label="online AR-bootstrap CI on Δ")
        ax.plot(x_obs, m, color="#2f6f4e", lw=1.8, label="running mean of Δ")
        # B=1 score is too noisy to overlay at full opacity
        ax.plot(
            x_obs,
            s,
            color="#1f4e79",
            lw=0.6 if batch == 1 else 1.2,
            alpha=0.25 if batch == 1 else 0.55,
            label=r"batch excess $\Delta=s_b-\mu_{\mathrm{ref}}$",
        )
        ax.axhline(0.0, color="#a33b3b", lw=1.0, alpha=0.8, label="null 0 after D_ref adjust")
        ax.axvline(N_REF + BURN * batch, color="#888", ls="--", lw=1)
        ax.set_title(
            f"Stationary · B={batch}  ·  μ_ref={summary['mu_ref']:.3f}\n"
            f"width={summary['mean_width']:.3f}  ·  "
            f"TV(lo)/1k={summary['tv_lo_per_1k']:.3f}  ·  "
            f"cover0={summary['cover0']:.2f}",
            fontsize=10,
        )
        ax.set_xlabel("observation index")
        ax.grid(True, alpha=0.3)
        # zoom y around the CI / mean, not the B=1 spikes
        y_core = np.concatenate([m, lo, hi, np.array([0.0])])
        y_core = y_core[np.isfinite(y_core)]
        if y_core.size:
            lo_y, hi_y = np.nanpercentile(y_core, [1, 99])
            pad = 0.2 * max(hi_y - lo_y, 1e-3)
            ax.set_ylim(lo_y - pad, hi_y + pad)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle(
        r"Stationary · null-adjusted online-bootstrap CI  ·  $\Delta_t=s_t-\mu_{\mathrm{ref}}$"
        f"\nPalm & Nagler AR-bootstrap · n_boot={N_BOOT} · "
        r"$\mu_{\mathrm{ref}}=$ mean PO-risk on $D_{\mathrm{ref}}$ mini-batches",
        fontsize=12,
        y=1.08,
    )
    axes[0, 0].set_ylabel("score")
    axes[1, 0].set_ylabel("score")
    fig.tight_layout()
    dest = OUT / "stationary_batchsize_bootstrap_ci.png"
    fig.savefig(dest, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}", flush=True)

    # summary: score noise + CI total variation per 1k obs
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))
    bs = list(BATCHES)
    axes[0].bar([str(b) for b in bs], [results[b][1]["std_s"] for b in bs], color="#1f4e79")
    axes[0].set_title(r"Excess noise  ·  std$(\Delta)$")
    axes[0].set_xlabel("batch size")
    axes[0].grid(True, axis="y", alpha=0.3)
    axes[1].bar([str(b) for b in bs], [results[b][1]["tv_lo_per_1k"] for b in bs], color="#c47b17")
    axes[1].set_title(r"CI twitch  ·  TV(lo) per 1000 observations")
    axes[1].set_xlabel("batch size")
    axes[1].grid(True, axis="y", alpha=0.3)
    fig.suptitle(
        r"Stationary · null-adjusted $\Delta=s-\mu_{\mathrm{ref}}$ · smaller B ⇒ twitchier CI",
        fontsize=11,
    )
    fig.tight_layout()
    dest2 = OUT / "stationary_batchsize_bootstrap_ci_summary.png"
    fig.savefig(dest2, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest2}", flush=True)

    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Stationary DGP with null adjustment. "
        r"$\Delta_t=s_t-\mu_{\mathrm{ref}}$, where $\mu_{\mathrm{ref}}$ is the mean frozen PO-risk "
        r"on $D_{\mathrm{ref}}$ mini-batches. Online AR-bootstrap CI ($n_{\mathrm{boot}}=200$) targets the "
        r"running mean of $\Delta$. cover0 $=$ fraction of post-burn batches whose CI contains $0$. "
        r"TV$(\mathrm{lo})/1\mathrm{k}$ is total variation of the lower endpoint per $1000$ observations.}",
        r"\label{tab:stationary-batchsize-boot-ci}",
        r"\small",
        r"\begin{tabular}{@{}r ccccc@{}}",
        r"\toprule",
        r"$B$ & $\mu_{\mathrm{ref}}$ & mean width & std$(\Delta)$ & TV$(\mathrm{lo})/1\mathrm{k}$ & cover0 \\",
        r"\midrule",
    ]
    for b in BATCHES:
        s = results[b][1]
        lines.append(
            f"{b} & {s['mu_ref']:.3f} & {s['mean_width']:.3f} & {s['std_s']:.3f} & "
            f"{s['tv_lo_per_1k']:.3f} & {s['cover0']:.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    tex = DOCS / "Streaming_PORisk_Stationary_BatchSize_BootstrapCI_tables_only.tex"
    tex.write_text("\n".join(lines))
    print(f"wrote {tex}", flush=True)

    for b in BATCHES:
        s = results[b][1]
        print(
            f"B={b}: mu_ref={s['mu_ref']:.3f} width={s['mean_width']:.3f} "
            f"std_d={s['std_s']:.3f} TV_lo/1k={s['tv_lo_per_1k']:.3f} cover0={s['cover0']:.3f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
