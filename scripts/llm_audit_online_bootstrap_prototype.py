#!/usr/bin/env python3
"""Palm–Nagler online-bootstrap + last-two hop overlay on LLM-audit labels.

Aligns the concept-board protocol with the audit tables:

    μ_ref = mean frozen-probe Brier on D_ref mini-batches
    Δ_t   = s_t − μ_ref
    fire  iff online AR-bootstrap CI_lo(Δ) > 0

Last-two OnlineRFPerm ``hop_fires`` is reported on the same stream as a
related but distinct gate (refit every hop, ratio ≥ γ).

Streams
-------
HH helpful/harmless × consistent/hop (already on disk).
BeaverTails / WildGuard / ToxicChat native labels, plus an in-memory
hop overlay (flip Y after ``cut_batch`` with p=0.92). HH chosen is not Y.

Usage::

    PYTHONPATH=. python3 scripts/llm_audit_online_bootstrap_prototype.py
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path = [p for p in sys.path if "/workspace/datasets" not in p]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from agod.online_ar_bootstrap import (  # noqa: E402
    BETA,
    first_significant,
    run_delta_bootstrap,
)
from agod.rf_probe import (  # noqa: E402
    brier_score,
    error_floor,
    fit_online_rf,
    hop_fires,
    probe_err,
    shift_ratio,
)

XY = ROOT / "results" / "manuscript" / "llm_audit"
OUT = ROOT / "results" / "manuscript" / "llm_audit_online_bootstrap"
DOCS = ROOT / "docs" / "manuscript"

X_COLS = [
    "x_n_toks",
    "x_n_chars",
    "x_avg_word",
    "x_qmark",
    "x_bang",
    "x_hedge",
    "x_formal",
    "x_i_count",
    "x_newlines",
    "x_upper",
    "x_refuse",
    "x_please",
    "x_thank",
]

N_BOOT = 200
ALPHA = 0.05
GATE = 1.5
N_REF_BATCHES = 4
CUT_BATCH = 4
FLIP_RATE = 0.92
BURN = 2
SEED = 2026


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if obj is None:
        return None
    return obj


def load_xy(path: Path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"empty table {path}")
    keys = {k.lower() for k in rows[0]}
    if "chosen" in keys or "rejected" in keys:
        raise RuntimeError(f"{path} still has chosen/rejected — that is not Y")
    y = np.asarray([int(r["y"]) for r in rows], dtype=int)
    batch = np.asarray([int(r["batch"]) for r in rows], dtype=int)
    X = np.asarray([[float(r[c]) for c in X_COLS] for r in rows], dtype=float)
    return X, y, batch


def apply_preference_hop(y, batch, *, cut_batch, flip_rate=0.92, seed=0):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=int).copy()
    batch = np.asarray(batch, dtype=int)
    after = batch >= int(cut_batch)
    flip = after & (rng.random(len(y)) < float(flip_rate))
    y[flip] = 1 - y[flip]
    return y, float(flip.mean()) if after.any() else 0.0


def freeze_and_score(X, y, batch, *, n_ref_batches: int, seed: int):
    """Frozen-ref Brier, concept-board protocol, classification probe."""
    batch = np.asarray(batch, dtype=int)
    ref = batch < int(n_ref_batches)
    if not np.any(ref) or np.all(ref):
        raise RuntimeError("need both D_ref and a trail")
    probe = fit_online_rf(X[ref], y[ref], seed=seed, task="acc")

    def one_score(b: int) -> float:
        m = batch == int(b)
        return brier_score(probe, X[m], y[m])

    ref_batches = sorted(int(b) for b in np.unique(batch[ref]))
    ref_scores = [one_score(b) for b in ref_batches]
    mu_ref = float(np.mean(ref_scores))
    trail_batches = sorted(int(b) for b in np.unique(batch[~ref]))
    scores = np.asarray([one_score(b) for b in trail_batches], dtype=float)
    err01 = np.asarray(
        [probe_err(probe, X[batch == b], y[batch == b], task="acc") for b in trail_batches],
        dtype=float,
    )
    return {
        "probe": probe,
        "mu_ref": mu_ref,
        "ref_scores": [float(s) for s in ref_scores],
        "ref_batches": ref_batches,
        "trail_batches": trail_batches,
        "scores": scores,
        "err01": err01,
        "n_ref": int(ref.sum()),
        "n_trail": int((~ref).sum()),
    }


def last_two_hops(X, y, batch, *, gate: float, seed: int):
    """Consecutive OOS hop_fires over the full stream (not frozen-ref)."""
    batch = np.asarray(batch, dtype=int)
    batches = sorted(int(b) for b in np.unique(batch))
    e_prev = None
    hist = []
    for t in batches[1:]:
        prev = batch == (t - 1)
        cur = batch == t
        if not np.any(prev) or not np.any(cur):
            continue
        probe = fit_online_rf(X[prev], y[prev], seed=int(seed) + t, task="acc")
        e_now = probe_err(probe, X[cur], y[cur], task="acc")
        e_fl = error_floor("acc", int(cur.sum()))
        fired = hop_fires(e_now, e_prev, gate=gate, e_floor=e_fl)
        ratio = 1.0 if e_prev is None else shift_ratio(e_now, e_prev, e_floor=e_fl)
        hist.append(
            {
                "abs_batch": int(t),
                "fired": bool(fired),
                "e_now": float(e_now),
                "e_prev": None if e_prev is None else float(e_prev),
                "ratio": float(ratio),
            }
        )
        e_prev = e_now
    return hist


def summarize_boot(rows, *, mu_ref: float, trail_batches, cut_batch: int | None, burn: int):
    post_burn = [r for r in rows if r["t"] > int(burn) and np.isfinite(r.get("width", np.nan))]
    labeled = cut_batch is not None
    pre, post = [], []
    for r, abs_b in zip(rows, trail_batches):
        rec = {**r, "abs_batch": int(abs_b)}
        if labeled and abs_b < int(cut_batch):
            pre.append(rec)
        else:
            post.append(rec)
    # trail starts at cut when D_ref = batches < cut, so pre is empty by design
    first = first_significant(rows, after_t=0)
    first_post = None
    if labeled:
        for r, abs_b in zip(rows, trail_batches):
            if abs_b >= int(cut_batch) and r.get("fire"):
                first_post = {**r, "abs_batch": int(abs_b)}
                break
    else:
        first_post = {**first, "abs_batch": int(trail_batches[first["t"] - 1])} if first else None

    def _cover0(block):
        if not block:
            return float("nan")
        return float(np.mean([(r["lo"] <= 0.0) and (r["hi"] >= 0.0) for r in block if np.isfinite(r["lo"])]))

    det = (
        float(np.mean([r["fire"] for r in post if np.isfinite(r["lo"])]))
        if any(np.isfinite(r["lo"]) for r in post)
        else float("nan")
    )
    delay = float("nan")
    if labeled and first_post is not None:
        delay = float(first_post["abs_batch"] - int(cut_batch))
    return {
        "mu_ref": float(mu_ref),
        "n_trail_batches": len(rows),
        "mean_width": float(np.nanmean([r["width"] for r in post_burn])) if post_burn else float("nan"),
        "cover0_all": _cover0(post_burn),
        "det_rate": det,
        "delay_batches": delay,
        "first_sig_t": None if first is None else int(first["t"]),
        "first_sig_abs_batch": None if first_post is None else int(first_post["abs_batch"]),
        "n_fires": int(sum(1 for r in rows if r.get("fire"))),
        "mean_delta": float(np.mean([r["s"] for r in rows])) if rows else float("nan"),
        "mean_post_delta": float(np.mean([r["s"] for r in post])) if post else float("nan"),
    }


def run_stream(
    name: str,
    title: str,
    X,
    y,
    batch,
    *,
    regime: str,
    n_ref_batches: int,
    cut_batch: int,
    gate: float,
    seed: int,
    labeled_cut: bool,
):
    frozen = freeze_and_score(X, y, batch, n_ref_batches=n_ref_batches, seed=seed)
    delta = frozen["scores"] - frozen["mu_ref"]
    rows, _boot = run_delta_bootstrap(delta, n_boot=N_BOOT, seed=seed, alpha=ALPHA)
    for r, abs_b, err01 in zip(rows, frozen["trail_batches"], frozen["err01"]):
        r["abs_batch"] = int(abs_b)
        r["err01"] = float(err01)
        r["s_raw"] = float(frozen["scores"][r["t"] - 1])
    hops = last_two_hops(X, y, batch, gate=gate, seed=seed)
    hop_at_cut = next((h for h in hops if h["abs_batch"] == int(cut_batch)), None)
    first_hop = next((h for h in hops if h["fired"]), None)
    cut = int(cut_batch) if labeled_cut else None
    boot_sum = summarize_boot(
        rows,
        mu_ref=frozen["mu_ref"],
        trail_batches=frozen["trail_batches"],
        cut_batch=cut,
        burn=BURN,
    )
    rec = {
        "name": name,
        "title": title,
        "regime": regime,
        "n": int(len(y)),
        "n_batches": int(batch.max() + 1),
        "y_pass_rate": float(np.mean(y)),
        "n_ref_batches": int(n_ref_batches),
        "cut_batch": int(cut_batch) if labeled_cut else None,
        "gate": float(gate),
        "beta": float(BETA),
        "n_boot": N_BOOT,
        "alpha": ALPHA,
        "score": "frozen_rf_brier",
        "mu_ref": frozen["mu_ref"],
        "ref_scores": frozen["ref_scores"],
        "bootstrap": boot_sum,
        "rows": rows,
        "hops": hops,
        "hop_fire_at_cut": None if hop_at_cut is None else bool(hop_at_cut["fired"]),
        "n_hop_fires": int(sum(1 for h in hops if h["fired"])),
        "first_hop_abs_batch": None if first_hop is None else int(first_hop["abs_batch"]),
        "cut_window_hops": [h for h in hops if abs(h["abs_batch"] - int(cut_batch)) <= 2],
    }
    return rec


def _fmt(x):
    if x is None:
        return "—"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float):
        if not np.isfinite(x):
            return "—"
        return f"{x:.3f}"
    return str(x)


def render_md(runs, *, n_ref_batches, cut_batch, gate, flip_rate) -> str:
    lines = [
        "# LLM audit — Palm–Nagler online bootstrap prototype",
        "",
        "Two gates on the same `(Y, X, batch)` tables. HH chosen is not Y.",
        "",
        "## Alignment",
        "",
        "| | Frozen-ref online AR-bootstrap | Last-two `hop_fires` |",
        "|---|---|---|",
        "| Probe | Fit once on `D_ref` (batches `< n_ref`) | Refit on `B_{t-1}` every hop |",
        "| Score `s_t` | Brier of the frozen RF on `B_t` | 0-1 OOS error `e_now` |",
        r"| Null | \(\mu_{\mathrm{ref}}\) = mean Brier on `D_ref` mini-batches | previous hop `e_prev` |",
        r"| Statistic | \(\Delta_t = s_t - \mu_{\mathrm{ref}}\) | `e_now / e_prev` |",
        "| UQ | Palm & Nagler AR-bootstrap CI | none (hard gate) |",
        r"| Fire | \(\mathrm{CI}_{\mathrm{lo}} > 0\) | ratio \(\ge \gamma\) and `e_prev ≥ e_floor` |",
        "| Needs | ≥2 trail updates before a CI exists | first hop always quiet (`e_prev is None`) |",
        "",
        f"Defaults: `n_ref_batches={n_ref_batches}`, labeled cut at batch `{cut_batch}`, "
        f"`γ={gate}`, `n_boot={N_BOOT}`, `β=√2−1`, hop overlay `p={flip_rate}`.",
        "",
        "Frozen-ref answers: is current auditor error *significantly above* the reference window?",
        "Last-two answers: did `P(Y|X)` hop between adjacent windows?",
        "",
        "## Results",
        "",
        "| Stream | Regime | pass | μ_ref | mean Δ | boot fires | first CI_lo>0 | delay | hop@cut | n hop fires |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|---:|",
    ]
    for r in runs:
        b = r["bootstrap"]
        lines.append(
            "| {title} | {regime} | {pass_} | {mu} | {md} | {nf} | {fs} | {delay} | {hop} | {nh} |".format(
                title=r["title"],
                regime=r["regime"],
                pass_=f"{r['y_pass_rate']:.2f}",
                mu=f"{r['mu_ref']:.3f}",
                md=_fmt(b["mean_delta"]),
                nf=b["n_fires"],
                fs="—" if b["first_sig_abs_batch"] is None else b["first_sig_abs_batch"],
                delay=_fmt(b["delay_batches"]),
                hop=_fmt(r["hop_fire_at_cut"]) if r["regime"] != "native" else "—",
                nh=r["n_hop_fires"],
            )
        )
    lines += [
        "",
        "Hop regime: expect CI_lo>0 after the cut (delay can be 1+ because the CI needs two trail points) "
        "and last-two fire at the cut. Consistent / native: CI should cover 0; last-two should stay mostly quiet.",
        "",
        "## Cut-window last-two log",
        "",
    ]
    for r in runs:
        if r["regime"] == "native":
            continue
        lines.append(f"### {r['title']} — `{r['regime']}`")
        lines.append("")
        lines.append("| abs batch | fire | e_prev | e_now | ratio |")
        lines.append("|---:|---|---:|---:|---:|")
        for h in r.get("cut_window_hops") or []:
            lines.append(
                "| {b} | {f} | {p} | {n} | {ratio} |".format(
                    b=h["abs_batch"],
                    f="yes" if h["fired"] else "no",
                    p=_fmt(h["e_prev"]),
                    n=_fmt(h["e_now"]),
                    ratio=_fmt(h["ratio"]),
                )
            )
        lines.append("")
    lines += [
        "## Serving",
        "",
        "| Gate | Quiet | Fire |",
        "|---|---|---|",
        "| AR-bootstrap `CI_lo>0` | current error is compatible with `D_ref` | frozen auditor is significantly worse; do not treat current judge as gold |",
        "| last-two hop | adjacent windows still the same map | policy-pack / judge swap between the last two batches; Top-k re-review |",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/llm_audit_online_bootstrap_prototype.py",
        "```",
        "",
    ]
    return "\n".join(lines) + "\n"


def plot_grid(runs, dest: Path, *, title: str, ncols: int = 2):
    n = len(runs)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.6 * ncols, 3.4 * nrows), squeeze=False)
    for ax, rec in zip(axes.ravel(), runs):
        rows = rec["rows"]
        x = np.asarray([r["abs_batch"] for r in rows], float)
        m = np.asarray([r["mean"] for r in rows], float)
        lo = np.asarray([r["lo"] for r in rows], float)
        hi = np.asarray([r["hi"] for r in rows], float)
        s = np.asarray([r["s"] for r in rows], float)
        ax.fill_between(x, lo, hi, color="#9ecae1", alpha=0.55, label="AR-bootstrap CI on Δ")
        ax.plot(x, m, color="#2f6f4e", lw=1.8, label="running mean of Δ")
        ax.plot(x, s, color="#1f4e79", lw=1.1, alpha=0.55, label=r"$\Delta_t=s_t-\mu_{\mathrm{ref}}$")
        ax.axhline(0.0, color="#a33b3b", lw=1.0, alpha=0.85, label="null 0")
        if rec.get("cut_batch") is not None:
            ax.axvline(rec["cut_batch"], color="#c47b17", ls=":", lw=1.4, label="cut")
        fs = rec["bootstrap"].get("first_sig_abs_batch")
        if fs is not None:
            ax.axvline(fs, color="#2f6f4e", ls="-.", lw=1.2, alpha=0.9)
        b = rec["bootstrap"]
        ax.set_title(
            f"{rec['title']} · {rec['regime']}\n"
            f"μ_ref={rec['mu_ref']:.3f}  fires={b['n_fires']}  "
            f"delay={_fmt(b['delay_batches'])}",
            fontsize=10,
        )
        ax.set_xlabel("batch")
        ax.grid(True, alpha=0.3)
        y_core = np.concatenate([m, lo, hi, np.array([0.0])])
        y_core = y_core[np.isfinite(y_core)]
        if y_core.size:
            lo_y, hi_y = np.nanpercentile(y_core, [1, 99])
            pad = 0.2 * max(hi_y - lo_y, 1e-3)
            ax.set_ylim(lo_y - pad, hi_y + pad)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(title, fontsize=12, y=1.08)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_tex_table(runs, dest: Path):
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{LLM-audit streams. Frozen-ref Palm--Nagler AR-bootstrap on "
        r"$\Delta_t=s_t-\mu_{\mathrm{ref}}$ (Brier of a shallow RF fitted on $D_{\mathrm{ref}}$). "
        r"Fire $=$ $\mathrm{CI}_{\mathrm{lo}}>0$. Delay $=$ batches from the labeled cut to the first fire. "
        r"Last-two hop@cut is the OnlineRFPerm consecutive-OOS gate, not the bootstrap.}",
        r"\label{tab:llm-audit-online-bootstrap}",
        r"\small",
        r"\begin{tabular}{@{}l l r r r r r@{}}",
        r"\toprule",
        r"Stream & regime & $\mu_{\mathrm{ref}}$ & mean $\Delta$ & boot fires & delay & hop@cut \\",
        r"\midrule",
    ]
    for r in runs:
        b = r["bootstrap"]
        delay_s = _fmt(b["delay_batches"])
        hop_s = _fmt(r["hop_fire_at_cut"]) if r["regime"] != "native" else "---"
        lines.append(
            f"{r['title']} & {r['regime']} & {r['mu_ref']:.3f} & {b['mean_delta']:.3f} & "
            f"{b['n_fires']} & {delay_s} & {hop_s} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    dest.write_text("\n".join(lines))


SPECS = [
    {
        "name": "hh_helpful_consistent",
        "title": "HH helpful",
        "file": "xy_hh_helpful_consistent.csv",
        "regime": "consistent",
        "overlay_hop": False,
        "labeled_cut": True,
        "seed": SEED,
    },
    {
        "name": "hh_helpful_hop",
        "title": "HH helpful",
        "file": "xy_hh_helpful_hop.csv",
        "regime": "hop",
        "overlay_hop": False,
        "labeled_cut": True,
        "seed": SEED,
    },
    {
        "name": "hh_harmless_consistent",
        "title": "HH harmless",
        "file": "xy_hh_harmless_consistent.csv",
        "regime": "consistent",
        "overlay_hop": False,
        "labeled_cut": True,
        "seed": SEED + 10,
    },
    {
        "name": "hh_harmless_hop",
        "title": "HH harmless",
        "file": "xy_hh_harmless_hop.csv",
        "regime": "hop",
        "overlay_hop": False,
        "labeled_cut": True,
        "seed": SEED + 10,
    },
    {
        "name": "beavertails_native",
        "title": "BeaverTails",
        "file": "xy_beavertails.csv",
        "regime": "native",
        "overlay_hop": False,
        "labeled_cut": False,
        "seed": SEED + 20,
    },
    {
        "name": "beavertails_hop",
        "title": "BeaverTails",
        "file": "xy_beavertails.csv",
        "regime": "hop",
        "overlay_hop": True,
        "labeled_cut": True,
        "seed": SEED + 20,
    },
    {
        "name": "wildguard_native",
        "title": "WildGuard",
        "file": "xy_wildguard.csv",
        "regime": "native",
        "overlay_hop": False,
        "labeled_cut": False,
        "seed": SEED + 30,
    },
    {
        "name": "wildguard_hop",
        "title": "WildGuard",
        "file": "xy_wildguard.csv",
        "regime": "hop",
        "overlay_hop": True,
        "labeled_cut": True,
        "seed": SEED + 30,
    },
    {
        "name": "toxicchat_native",
        "title": "ToxicChat",
        "file": "xy_toxicchat.csv",
        "regime": "native",
        "overlay_hop": False,
        "labeled_cut": False,
        "seed": SEED + 40,
    },
    {
        "name": "toxicchat_hop",
        "title": "ToxicChat",
        "file": "xy_toxicchat.csv",
        "regime": "hop",
        "overlay_hop": True,
        "labeled_cut": True,
        "seed": SEED + 40,
    },
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-ref-batches", type=int, default=N_REF_BATCHES)
    ap.add_argument("--cut-batch", type=int, default=CUT_BATCH)
    ap.add_argument("--gate", type=float, default=GATE)
    ap.add_argument("--flip-rate", type=float, default=FLIP_RATE)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    runs = []
    for spec in SPECS:
        path = XY / spec["file"]
        print(f"run {spec['name']} from {path.name}", flush=True)
        X, y, batch = load_xy(path)
        if spec["overlay_hop"]:
            y, _ = apply_preference_hop(
                y,
                batch,
                cut_batch=args.cut_batch,
                flip_rate=args.flip_rate,
                seed=spec["seed"] + 17,
            )
        rec = run_stream(
            spec["name"],
            spec["title"],
            X,
            y,
            batch,
            regime=spec["regime"],
            n_ref_batches=args.n_ref_batches,
            cut_batch=args.cut_batch,
            gate=args.gate,
            seed=spec["seed"],
            labeled_cut=spec["labeled_cut"],
        )
        rec["file"] = spec["file"]
        rec["y_is_chosen"] = False
        runs.append(rec)
        b = rec["bootstrap"]
        print(
            f"  μ_ref={rec['mu_ref']:.3f} meanΔ={b['mean_delta']:.3f} "
            f"boot_fires={b['n_fires']} first_sig={b['first_sig_abs_batch']} "
            f"hop@cut={rec['hop_fire_at_cut']} n_hop={rec['n_hop_fires']}",
            flush=True,
        )
        slim = {k: v for k, v in rec.items() if k not in {"rows", "hops"}}
        slim["rows"] = rec["rows"]
        slim["hops"] = rec["hops"]
        (OUT / f"{spec['name']}.json").write_text(json.dumps(jsonable(slim), indent=2) + "\n")

    hh = [r for r in runs if r["name"].startswith("hh_")]
    real_native = [r for r in runs if r["name"].endswith("_native")]
    real_hop = [r for r in runs if r["name"] in {"beavertails_hop", "wildguard_hop", "toxicchat_hop"}]
    plot_grid(
        hh,
        OUT / "hh_online_bootstrap_ci.png",
        title=r"HH auditor · frozen-ref $\Delta_t=s_t-\mu_{\mathrm{ref}}$ · Palm & Nagler",
        ncols=2,
    )
    plot_grid(
        real_native,
        OUT / "real_native_online_bootstrap_ci.png",
        title=r"Real labels (native) · frozen-ref $\Delta$ · no synthetic hop",
        ncols=3,
    )
    plot_grid(
        real_hop,
        OUT / "real_hop_online_bootstrap_ci.png",
        title=r"Real labels + Y-flip hop overlay · frozen-ref $\Delta$",
        ncols=3,
    )

    summary = {
        "protocol": "frozen_ref_delta_ar_bootstrap",
        "last_two_gate": "hop_fires",
        "n_ref_batches": args.n_ref_batches,
        "cut_batch": args.cut_batch,
        "gate": args.gate,
        "flip_rate": args.flip_rate,
        "beta": float(BETA),
        "n_boot": N_BOOT,
        "alpha": ALPHA,
        "y_is_chosen": False,
        "runs": jsonable(
            [
                {
                    k: r[k]
                    for k in (
                        "name",
                        "title",
                        "regime",
                        "file",
                        "n",
                        "n_batches",
                        "y_pass_rate",
                        "mu_ref",
                        "bootstrap",
                        "hop_fire_at_cut",
                        "n_hop_fires",
                        "first_hop_abs_batch",
                        "y_is_chosen",
                    )
                }
                for r in runs
            ]
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    md = render_md(
        runs,
        n_ref_batches=args.n_ref_batches,
        cut_batch=args.cut_batch,
        gate=args.gate,
        flip_rate=args.flip_rate,
    )
    (OUT / "REPORT.md").write_text(md)
    write_tex_table(runs, OUT / "online_bootstrap_tables_only.tex")
    (DOCS / "use_case_04_online_bootstrap.md").write_text(md)
    print(md)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
