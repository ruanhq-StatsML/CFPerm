#!/usr/bin/env python3
"""Concept-drift DGP · null-adjusted online-bootstrap CI vs batch size.

Same protocol as the stationary board:
  μ_ref = mean frozen PO-risk on D_ref mini-batches
  Δ_t = s_t − μ_ref
  Palm & Nagler AR-bootstrap on Δ

On concept drift, larger B makes the post-change excess readable (CI leaves 0).
Too-small B: twitchy stream UQ, hard to call a significant shift.

  python3 scripts/plot_concept_batchsize_bootstrap_ci.py
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
T_CHANGE_OBS = N_REF + 1500  # absolute observation index of concept jump
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


def concept_stream(*, n: int, seed: int, p: int = 8, t_change: int = T_CHANGE_OBS):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    beta = rng.normal(size=p)
    signal = X @ beta
    Y = signal + rng.normal(size=n)
    # Concept drift: P(Y|X) jumps after t_change
    Y[t_change:] = signal[t_change:] + 2.5 * X[t_change:, 0] + rng.normal(size=n - t_change)
    return X, Y, t_change


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


def run_batch(batch: int):
    X, Y, t_change = concept_stream(n=N_REF + N_TRAIL, seed=SEED)
    times, scores, mu_ref = freeze_and_score(X, Y, n_ref=N_REF, batch=batch)
    delta = scores - mu_ref
    boot = OnlineARBootstrap(n_boot=N_BOOT, seed=SEED + batch)
    rows = []
    for t, d in zip(times, delta):
        boot.update(float(d))
        lo, hi = boot.ci(ALPHA)
        # observation index of this batch end
        obs_end = N_REF + int(t) * batch
        rows.append(
            dict(
                t=int(t),
                obs=int(obs_end),
                s=float(d),
                mean=float(boot.mean),
                lo=float(lo),
                hi=float(hi),
                width=float(hi - lo),
                post=bool(obs_end > t_change),
            )
        )

    post_burn = [r for r in rows if r["t"] > BURN and np.isfinite(r["width"])]
    pre = [r for r in post_burn if not r["post"]]
    post = [r for r in post_burn if r["post"]]

    def _tv_lo(block, batch_size):
        if len(block) < 2:
            return float("nan")
        los = np.asarray([r["lo"] for r in block], float)
        n_obs = float(max(len(block) * batch_size, 1))
        return float(np.nansum(np.abs(np.diff(los))) / n_obs * 1000.0)

    # first post-change batch where CI entirely above 0 (significant excess)
    first_sig = None
    for r in post:
        if np.isfinite(r["lo"]) and r["lo"] > 0:
            first_sig = r
            break
    delay_batches = (
        float(first_sig["t"] - max(1, int(np.ceil((t_change - N_REF) / batch))))
        if first_sig is not None
        else float("nan")
    )
    # DetRate: fraction of post-change batches with lo>0
    det_rate = float(np.mean([r["lo"] > 0 for r in post])) if post else float("nan")
    cover0_pre = (
        float(np.mean([(r["lo"] <= 0) and (r["hi"] >= 0) for r in pre])) if pre else float("nan")
    )

    summary = dict(
        batch=batch,
        mu_ref=float(mu_ref),
        t_change=int(t_change),
        t_change_batch=int(np.ceil((t_change - N_REF) / batch)),
        n_batches=len(rows),
        mean_width=float(np.nanmean([r["width"] for r in post_burn])),
        std_pre=float(np.nanstd([r["s"] for r in pre])) if pre else float("nan"),
        std_post=float(np.nanstd([r["s"] for r in post])) if post else float("nan"),
        tv_lo_per_1k=_tv_lo(post_burn, batch),
        cover0_pre=cover0_pre,
        det_rate=det_rate,
        delay_batches=delay_batches,
        first_sig_obs=float(first_sig["obs"]) if first_sig is not None else float("nan"),
        mean_post=float(np.nanmean([r["mean"] for r in post])) if post else float("nan"),
    )
    return rows, summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    print("running concept-drift batch-size CI grid", flush=True)
    results = {b: run_batch(b) for b in BATCHES}

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.0), sharey=False)
    for ax, batch in zip(axes.ravel(), BATCHES):
        rows, summary = results[batch]
        t = np.asarray([r["t"] for r in rows], float)
        s = np.asarray([r["s"] for r in rows], float)
        m = np.asarray([r["mean"] for r in rows], float)
        lo = np.asarray([r["lo"] for r in rows], float)
        hi = np.asarray([r["hi"] for r in rows], float)
        x_obs = np.asarray([r["obs"] for r in rows], float)
        ax.fill_between(x_obs, lo, hi, color="#9ecae1", alpha=0.55, label="online AR-bootstrap CI on Δ")
        ax.plot(x_obs, m, color="#2f6f4e", lw=1.8, label="running mean of Δ")
        ax.plot(
            x_obs,
            s,
            color="#1f4e79",
            lw=0.6 if batch == 1 else 1.2,
            alpha=0.22 if batch == 1 else 0.50,
            label=r"batch excess $\Delta=s_b-\mu_{\mathrm{ref}}$",
        )
        ax.axhline(0.0, color="#a33b3b", lw=1.0, alpha=0.85, label="null 0")
        ax.axvline(summary["t_change"], color="#c47b17", ls=":", lw=1.4, label="concept CP")
        ax.axvline(N_REF + BURN * batch, color="#888", ls="--", lw=1)
        if np.isfinite(summary["first_sig_obs"]):
            ax.axvline(summary["first_sig_obs"], color="#2f6f4e", ls="-.", lw=1.2, alpha=0.9)
        det = summary["det_rate"]
        delay = summary["delay_batches"]
        delay_s = f"{delay:.0f}" if np.isfinite(delay) else "---"
        ax.set_title(
            f"Concept · B={batch}  ·  μ_ref={summary['mu_ref']:.3f}\n"
            f"TV(lo)/1k={summary['tv_lo_per_1k']:.3f}  ·  "
            f"post det(lo>0)={det:.2f}  ·  delay={delay_s} batches",
            fontsize=10,
        )
        ax.set_xlabel("observation index")
        ax.grid(True, alpha=0.3)
        y_core = np.concatenate([m, lo, hi, np.array([0.0])])
        y_core = y_core[np.isfinite(y_core)]
        if y_core.size:
            lo_y, hi_y = np.nanpercentile(y_core, [1, 99])
            pad = 0.2 * max(hi_y - lo_y, 1e-3)
            ax.set_ylim(lo_y - pad, hi_y + pad)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle(
        r"Concept drift · null-adjusted online-bootstrap CI  ·  $\Delta_t=s_t-\mu_{\mathrm{ref}}$"
        f"\nPalm & Nagler · n_boot={N_BOOT} · smaller B ⇒ harder to read a significant shift",
        fontsize=12,
        y=1.08,
    )
    axes[0, 0].set_ylabel("excess")
    axes[1, 0].set_ylabel("excess")
    fig.tight_layout()
    dest = OUT / "concept_batchsize_bootstrap_ci.png"
    fig.savefig(dest, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest}", flush=True)

    # summary trends
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5))
    bs = list(BATCHES)
    axes[0].bar([str(b) for b in bs], [results[b][1]["tv_lo_per_1k"] for b in bs], color="#c47b17")
    axes[0].set_title(r"CI twitch  ·  TV(lo)/1k obs")
    axes[0].set_xlabel("batch size")
    axes[0].grid(True, axis="y", alpha=0.3)

    axes[1].bar([str(b) for b in bs], [results[b][1]["det_rate"] for b in bs], color="#1f4e79")
    axes[1].set_title(r"Post-CP significance  ·  fraction with CI$_{\mathrm{lo}}>0$")
    axes[1].set_xlabel("batch size")
    axes[1].set_ylim(0, 1.05)
    axes[1].grid(True, axis="y", alpha=0.3)

    delays = [results[b][1]["delay_batches"] for b in bs]
    axes[2].bar(
        [str(b) for b in bs],
        [d if np.isfinite(d) else 0 for d in delays],
        color="#2f6f4e",
    )
    axes[2].set_title("Delay to first significant CI (batches)")
    axes[2].set_xlabel("batch size")
    axes[2].grid(True, axis="y", alpha=0.3)
    fig.suptitle(
        r"Concept drift · null-adjusted $\Delta$ · batch-size trend",
        fontsize=11,
    )
    fig.tight_layout()
    dest2 = OUT / "concept_batchsize_bootstrap_ci_summary.png"
    fig.savefig(dest2, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {dest2}", flush=True)

    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Concept-drift DGP with null adjustment $\Delta_t=s_t-\mu_{\mathrm{ref}}$. "
        r"Online AR-bootstrap CI on $\Delta$. Post-CP det.\ $=$ fraction of post-change batches with "
        r"CI$_{\mathrm{lo}}>0$ (readable excess vs $D_{\mathrm{ref}}$). Delay $=$ batches from the labeled CP "
        r"to the first such significant CI. Small $B$ is twitchy and hard to call a shift.}",
        r"\label{tab:concept-batchsize-boot-ci}",
        r"\small",
        r"\begin{tabular}{@{}r ccccc@{}}",
        r"\toprule",
        r"$B$ & $\mu_{\mathrm{ref}}$ & TV$(\mathrm{lo})/1\mathrm{k}$ & cover0 pre & post det & delay \\",
        r"\midrule",
    ]
    for b in BATCHES:
        s = results[b][1]
        delay_s = f"{s['delay_batches']:.0f}" if np.isfinite(s["delay_batches"]) else "---"
        lines.append(
            f"{b} & {s['mu_ref']:.3f} & {s['tv_lo_per_1k']:.3f} & {s['cover0_pre']:.3f} & "
            f"{s['det_rate']:.3f} & {delay_s} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    tex = DOCS / "Streaming_PORisk_Concept_BatchSize_BootstrapCI_tables_only.tex"
    tex.write_text("\n".join(lines))
    print(f"wrote {tex}", flush=True)

    for b in BATCHES:
        s = results[b][1]
        delay_s = f"{s['delay_batches']:.0f}" if np.isfinite(s["delay_batches"]) else "---"
        print(
            f"B={b}: mu_ref={s['mu_ref']:.3f} TV/1k={s['tv_lo_per_1k']:.3f} "
            f"cover0_pre={s['cover0_pre']:.3f} det={s['det_rate']:.3f} delay={delay_s}",
            flush=True,
        )


if __name__ == "__main__":
    main()
