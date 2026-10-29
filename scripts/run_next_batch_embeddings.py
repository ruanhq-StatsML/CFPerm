#!/usr/bin/env python3
"""Run Wu InstDisc, Yu SDC, and typed GPM identification on typed streams.

  python3 scripts/run_next_batch_embeddings.py
  python3 scripts/run_next_batch_embeddings.py --quick
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from gpm_fsds import gpm_identification, run_gpm_method
from instance_discrimination import run_instdisc_stream
from msrvtt_attribution_plots import COLORS, GRID, INK, MUTED, _save, _style
from msrvtt_multimodal_attribution import GROUP_NAMES, write_json
from prototype_drift import run_prototype_stream
from typed_shift_stepsize import make_typed_stream

OUT = ROOT / "results" / "next_batch_embeddings"
DOCS = ROOT / "docs" / "method"


def _make(regime, n_batches, n_per, seed):
    if regime == "cov_only":
        return make_typed_stream(
            n_batches=n_batches,
            n_per=n_per,
            seed=seed,
            cov={"video": 0.18, "audio": 0.05, "text": 0.0},
        )
    return make_typed_stream(
        n_batches=n_batches,
        n_per=n_per,
        seed=seed,
        concept={"video": 1.0},
        concept_at=max(n_batches // 2, 3),
    )


def _mean_sd(rows, path):
    v = np.array([path(r) for r in rows], dtype=float)
    return float(np.nanmean(v)), float(np.nanstd(v, ddof=1) if len(v) > 1 else 0.0)


def run_suite(seeds, n_batches, n_per, dim=24, epochs=3, steps=4):
    table = {}
    traces = {}
    for regime in ("cov_only", "concept_only"):
        inst, proto, gpm_id, gpm_run = [], [], [], []
        traces[regime] = {"inst": [], "proto": [], "gpm": []}
        for seed in seeds:
            stream = _make(regime, n_batches, n_per, seed)
            rec_i, _ = run_instdisc_stream(stream, dim=dim, epochs=epochs, seed=seed)
            rec_p = run_prototype_stream(stream, seed=seed)
            rec_g = gpm_identification(stream)
            modes = {}
            for mode in ("none", "always", "typed"):
                modes[mode] = run_gpm_method(
                    stream, mode=mode, steps_per_batch=steps, warmup_steps=2, seed=seed
                )
            inst.append(rec_i)
            proto.append(rec_p)
            gpm_id.append(rec_g)
            gpm_run.append(modes)
            traces[regime]["inst"].append(rec_i)
            traces[regime]["proto"].append(rec_p)
            traces[regime]["gpm"].append(modes)
        cell = {"instdisc": {}, "prototype": {}, "gpm": {}}
        for g in GROUP_NAMES:
            cell["instdisc"][g] = {
                "acc": _mean_sd(inst, lambda r, h=g: r["mean_acc"][h]),
                "stale": _mean_sd(inst, lambda r, h=g: r["mean_stale"][h]),
                "collision": _mean_sd(inst, lambda r, h=g: r["mean_collision"][h]),
            }
            cell["prototype"][g] = {
                "ncm_self": _mean_sd(proto, lambda r, h=g: r["mean_ncm_self"][h]),
                "ncm_stale": _mean_sd(proto, lambda r, h=g: r["mean_ncm_stale"][h]),
                "ncm_sdc": _mean_sd(proto, lambda r, h=g: r["mean_ncm_sdc"][h]),
                "proto_cos_sdc": _mean_sd(proto, lambda r, h=g: r["mean_proto_cos_sdc"][h]),
                "proto_cos_stale": _mean_sd(proto, lambda r, h=g: r["mean_proto_cos_stale"][h]),
            }
            cell["gpm"][g] = {
                "residual": _mean_sd(gpm_id, lambda r, h=g: r["mean_residual"][h]),
                "c": _mean_sd(gpm_id, lambda r, h=g: r["mean_c"][h]),
            }
        for mode in ("none", "always", "typed"):
            cell["gpm"][mode] = {
                "bwt": _mean_sd(gpm_run, lambda r, m=mode: r[m]["bwt"]),
                "post_acc": _mean_sd(gpm_run, lambda r, m=mode: r[m]["post_acc"]),
                "online_acc": _mean_sd(gpm_run, lambda r, m=mode: r[m]["online_acc"]),
            }
        table[regime] = cell
    return {"table": table, "traces": traces, "seeds": [int(s) for s in seeds]}


def plot_suite(suite, path):
    _style()
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(11.4, 6.2))
    regimes = (("cov_only", "Covariate only"), ("concept_only", "Concept only"))
    for i, (regime, title) in enumerate(regimes):
        inst = suite["traces"][regime]["inst"]
        proto = suite["traces"][regime]["proto"]
        gpm = suite["traces"][regime]["gpm"]
        t = [h["round"] for h in inst[0]["history"][1:]]
        ax = axes[i, 0]
        for g in GROUP_NAMES:
            y = np.array([[h[g]["collision"] for h in r["history"][1:]] for r in inst])
            ax.plot(t, y.mean(0), color=COLORS[g], lw=2.0, label=g)
            ax.fill_between(t, y.mean(0) - y.std(0), y.mean(0) + y.std(0), color=COLORS[g], alpha=0.12, lw=0)
        ax.set_ylabel(title + "\ninstance collision")
        ax.set_xlabel("batch")
        ax.grid(True, color=GRID)
        if i == 0:
            ax.set_title("InstDisc collision", loc="left", fontsize=12, fontweight="bold")
            ax.legend(frameon=False, fontsize=8)
        ax = axes[i, 1]
        for g in GROUP_NAMES:
            y = np.array([[h[g]["proto_cos_sdc"] for h in r["history"][1:]] for r in proto])
            ax.plot(t, np.nanmean(y, 0), color=COLORS[g], lw=2.0, label=g)
        ax.set_ylabel("SDC proto cosine")
        ax.set_xlabel("batch")
        ax.grid(True, color=GRID)
        if i == 0:
            ax.set_title("Yu SDC prototypes", loc="left", fontsize=12, fontweight="bold")
        ax = axes[i, 2]
        for mode, ls, col in (
            ("none", "-", "#9AA3AE"),
            ("always", "--", "#A33B24"),
            ("typed", "-", "#C45C26"),
        ):
            y = np.array([[h["acc"] for h in r[mode]["history"]] for r in gpm], dtype=float)
            tt = [h["round"] for h in gpm[0][mode]["history"]]
            ax.plot(tt, y.mean(0), ls=ls, color=col, lw=2.0, label=mode)
        ax.set_ylabel("probe acc")
        ax.set_xlabel("batch")
        ax.grid(True, color=GRID)
        if i == 0:
            ax.set_title("GPM on the linear head", loc="left", fontsize=12, fontweight="bold")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Next-batch embeddings: Wu bank, SDC prototypes, typed GPM",
        fontsize=13.2,
        fontweight="bold",
        color=INK,
        x=0.04,
        ha="left",
    )
    fig.text(
        0.04,
        0.01,
        "InstDisc collision is off-diagonal cosine (location hop collapses instances). "
        "SDC moves old class means with the current-batch field. "
        "Typed GPM projects only when ĉ is loud and δ̂ is quiet.",
        fontsize=8.2,
        color=MUTED,
    )
    return _save(fig, path)


def write_tex_table(suite, path):
    table = suite["table"]
    lines = [
        r"% Next-batch embeddings. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{Wu instance-discrimination memory bank, Yu SDC prototypes, and GPM gated by FSDS/TSS $(\hat c,\hat\delta)$.",
        r"Entries are mean (s.d.) over seeds. Collision is off-diagonal cosine of the current batch (instances collapse under a location hop).",
        r"SDC prototype cosine is compensated old class means vs the true new means. GPM modes: no projection, always project, typed gate.}",
        r"\label{tab:next-batch-embeddings}",
        r"\small",
        r"\begin{tabular}{@{}ll ccc ccc@{}}\toprule",
        r"Regime & Head / mode & InstDisc coll. & SDC proto-cos & GPM residual & BWT & post-acc \\",
        r"\midrule",
    ]

    def fmt_pair(mu_sd):
        return r"$%.3f$ ($%.3f$)" % (mu_sd[0], mu_sd[1])

    names = {"cov_only": "covariate only", "concept_only": "concept only"}
    for regime in ("cov_only", "concept_only"):
        cell = table[regime]
        first = True
        for g in GROUP_NAMES:
            lines.append(
                r"%s & %s & %s & %s & %s & --- & --- \\"
                % (
                    names[regime] if first else "",
                    g,
                    fmt_pair(cell["instdisc"][g]["collision"]),
                    fmt_pair(cell["prototype"][g]["proto_cos_sdc"]),
                    fmt_pair(cell["gpm"][g]["residual"]),
                )
            )
            first = False
        for mode in ("none", "always", "typed"):
            lines.append(
                r" & GPM %s & --- & --- & --- & %s & %s \\"
                % (
                    mode,
                    fmt_pair(cell["gpm"][mode]["bwt"]),
                    fmt_pair(cell["gpm"][mode]["post_acc"]),
                )
            )
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    lines.extend(
        [
            r"\end{tabular}",
            r"\end{table}",
            "",
        ]
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    n_batches = 6 if args.quick else 8
    n_per = 36 if args.quick else 48
    seeds = list(range(2026, 2026 + (2 if args.quick else 4)))
    print("suite", seeds, n_batches, n_per, flush=True)
    suite = run_suite(seeds, n_batches, n_per, dim=16 if args.quick else 24, epochs=2 if args.quick else 3)
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {"table": suite["table"], "seeds": suite["seeds"]}
    write_json(OUT / "next_batch_embeddings.json", payload)
    plot_suite(suite, OUT / "next_batch_embeddings.png")
    write_tex_table(suite, OUT / "Next_batch_embeddings.tex")
    write_tex_table(suite, DOCS / "Next_batch_embeddings.tex")
    for regime, cell in suite["table"].items():
        print("==", regime, flush=True)
        for g in GROUP_NAMES:
            print(
                "  %s  coll=%.3f  sdc_cos=%.3f  gpm_res=%.3f"
                % (
                    g,
                    cell["instdisc"][g]["collision"][0],
                    cell["prototype"][g]["proto_cos_sdc"][0],
                    cell["gpm"][g]["residual"][0],
                ),
                flush=True,
            )
        for mode in ("none", "always", "typed"):
            print(
                "  GPM %s  bwt=%.3f  post=%.3f"
                % (mode, cell["gpm"][mode]["bwt"][0], cell["gpm"][mode]["post_acc"][0]),
                flush=True,
            )


if __name__ == "__main__":
    main()
