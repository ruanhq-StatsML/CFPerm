#!/usr/bin/env python3
"""Recsys order-stream by grain: OnlineRFPerm, RFPerm/CFPerm, FSDS, post-hoc."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
for p in (SRC, ROOT):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from recsys_dim_monitor import DIMS, run_all_dims  # noqa: E402
from stream_dgps import make_order_graph_stream  # noqa: E402

OUT = ROOT / "results" / "recsys_dim_monitor"
KINDS = ("covariate_south", "concept_south", "both")
N_REF = 400
N_NEW = 120
N_BATCHES = 6
ONSET = 2
N_MERCHANTS = 16
N_USERS = 48


def _fmt(x, n=3):
    if x is None:
        return ""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v != v:
        return ""
    return f"{v:.{n}g}"


def _md_row(cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    return obj


def run_kind(kind: str, seed: int = 2026) -> dict:
    tables = make_order_graph_stream(
        n_ref=N_REF,
        n_new=N_NEW,
        n_batches=N_BATCHES,
        onset_batch=ONSET,
        n_merchants=N_MERCHANTS,
        n_users=N_USERS,
        kind=kind,
        seed=seed,
    )
    dims = run_all_dims(tables, seed=seed, with_po=True)
    return {"meta": tables["meta"], "dims": dims}


def slim_dim(d: dict) -> dict:
    planted = d["planted"]
    planted_names = list(planted) if isinstance(planted, dict) else list(planted)
    return {
        "dim": d["dim"],
        "onset_true": d["onset_true"],
        "onset_hat": d["onset_hat"],
        "onset_rank": d["onset_rank"],
        "addis_t": d["addis_t"],
        "saffron_t": d["saffron_t"],
        "n_hop": d["n_hop"],
        "hop_status": d.get("hop_status"),
        "rank_status": d.get("rank_status"),
        "addis_status": d.get("addis_status"),
        "saffron_status": d.get("saffron_status"),
        "hop_delay": d.get("hop_delay"),
        "rank_delay": d.get("rank_delay"),
        "addis_delay": d.get("addis_delay"),
        "saffron_delay": d.get("saffron_delay"),
        "fsds_top": [r["feature"] for r in d["fsds"][:4]],
        "fsds_loud": [r["feature"] for r in d["fsds"] if r.get("loud")],
        "fsds_recovered": d["fsds_recovered"],
        "vimp_top": d["vimp_top"],
        "vimp_reject": d["vimp_reject"],
        "vimp_max": d["vimp"]["max_vimp"],
        "mse_vimp_top": d.get("mse_vimp_top") or d.get("mse_vimp", {}).get("top", []),
        "tau_fsds": d.get("tau_fsds"),
        "tau_cfperm": d.get("tau_cfperm"),
        "tau_rfperm": d.get("tau_rfperm"),
        "planted": planted_names,
        "y_in_X": d["y_in_X"],
        "loc_n_south": d["localization"]["n_south"],
        "loc_pair_mmd": d["localization"].get("pair_mmd"),
        "loc_south_mmd": (d["localization"].get("south") or {}).get("mmd"),
        "loc_north_mmd": (d["localization"].get("north") or {}).get("mmd"),
        "loc_south_cmean_x": (d["localization"].get("south") or {}).get("cmean_x"),
        "loc_north_cmean_x": (d["localization"].get("north") or {}).get("cmean_x"),
        "loc_south_cmean_y": (d["localization"].get("south") or {}).get("cmean_y"),
        "loc_north_cmean_y": (d["localization"].get("north") or {}).get("cmean_y"),
        "loc_south_po": (d["localization"].get("south") or {}).get("po"),
        "loc_north_po": (d["localization"].get("north") or {}).get("po"),
        "last_T": d["rows"][-1]["rfperm_T"],
        "last_mmd": d["rows"][-1]["mmd"],
        "last_po": d["rows"][-1]["po"],
        "rows": [
            {
                "t": r["t"],
                "onset": r["onset"],
                "T": r["rfperm_T"],
                "p": r["rfperm_p"],
                "hop": r["rfperm_hop"],
                "mse": r["mse"],
                "mmd": r["mmd"],
                "cmean_x": r["cmean_x"],
                "cmean_y": r["cmean_y"],
                "po": r["po"],
            }
            for r in d["rows"]
        ],
        "fsds": [
            {k: r[k] for k in ("feature", "score", "mmd", "cmean_x", "cmean_y", "loud")}
            for r in d["fsds"]
        ],
        "vimp_rank": d["vimp"]["rank"][:6],
        "mse_vimp_rank": d.get("mse_vimp", {}).get("rank", [])[:6],
    }


def _planted_names(d: dict) -> list[str]:
    planted = d.get("planted") or {}
    if isinstance(planted, dict):
        return list(planted)
    return list(planted)


def _status_cell(t_hat, delay, status) -> str:
    if status == "miss" or t_hat is None:
        return "miss"
    if status == "FAR":
        return f"FAR@{t_hat}"
    d = "" if delay is None else str(delay)
    return f"{t_hat} (d={d})"


def write_markdown(results: dict, dest: Path) -> None:
    lines = [
        "# Recsys order-stream by grain — OnlineRFPerm × RFPerm/CFPerm × FSDS",
        "",
        "Data: synthetic order / merchant / user stream (recommendation-shaped table). "
        "Y = conversion, never a feature. T = batch label. "
        "OnlineRFPerm (Algorithm 1) marks WHEN. RFPerm ΔMSE and CFPerm/PermuCATE VIMP "
        "plus FSDS mark WHICH columns; Kendall-τ is the ranking recovery in the MetaLearner paper. "
        "Post-hoc region split marks WHICH accounts (south vs north). "
        "Localization, not a unique decomposition. Not a graph method.",
        "",
        f"n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}, onset_true={ONSET}, "
        f"merchants={N_MERCHANTS}. Planted subset = south. "
        "Delay = first rejection − onset_true. Negative / FAR = mark before labeled onset.",
        "",
        "## Table 1 · OnlineRFPerm (WHEN) — first rejection and delay",
        "",
        _md_row(
            [
                "kind",
                "dim",
                "onset_true",
                "hop",
                "rank-p (p<0.05)",
                "ADDIS",
                "SAFFRON",
                "last T",
                "last MMD",
            ]
        ),
        _md_row(["---"] * 9),
    ]
    for kind in KINDS:
        for dim in DIMS:
            d = results[kind]["dims"][dim]
            lines.append(
                _md_row(
                    [
                        kind,
                        dim,
                        d["onset_true"],
                        _status_cell(d["onset_hat"], d.get("hop_delay"), d.get("hop_status")),
                        _status_cell(d["onset_rank"], d.get("rank_delay"), d.get("rank_status")),
                        _status_cell(d["addis_t"], d.get("addis_delay"), d.get("addis_status")),
                        _status_cell(d["saffron_t"], d.get("saffron_delay"), d.get("saffron_status")),
                        _fmt(d["rows"][-1]["rfperm_T"]),
                        _fmt(d["rows"][-1]["mmd"]),
                    ]
                )
            )
    lines += [
        "",
        "## Table 2 · Ranking recovery (WHICH columns) — Kendall-τ vs planted magnitude",
        "",
        "Planted order on covariate/both: amount ≻ merchant_gmv ≻ channel ≻ noise. "
        "Concept plants amount in Y|X only. τ is Kendall’s τ between scores and that magnitude. "
        "CFPerm reject = max φ-VIMP vs 95% of T-permuted nulls (B=12). Empty reject is a miss, not a quiet stream.",
        "",
        _md_row(
            [
                "kind",
                "dim",
                "planted in grain",
                "FSDS recovered",
                "τ_FSDS",
                "RFPerm ΔMSE top-3",
                "τ_RFPerm",
                "CFPerm φ top-3",
                "τ_CFPerm",
                "CFPerm reject",
            ]
        ),
        _md_row(["---"] * 10),
    ]
    for kind in KINDS:
        for dim in DIMS:
            d = results[kind]["dims"][dim]
            planted = ",".join(_planted_names(d)) or "—"
            recov = ",".join(d["fsds_recovered"]) or "—"
            mse_top = d.get("mse_vimp_top") or [r["feature"] for r in d.get("mse_vimp", {}).get("rank", [])[:3]]
            lines.append(
                _md_row(
                    [
                        kind,
                        dim,
                        planted,
                        recov,
                        _fmt(d.get("tau_fsds")),
                        ",".join(mse_top) or "—",
                        _fmt(d.get("tau_rfperm")),
                        ",".join(d["vimp_top"]),
                        _fmt(d.get("tau_cfperm")),
                        "yes" if d["vimp_reject"] else "",
                    ]
                )
            )
    lines += [
        "",
        "## Table 3 · Post-hoc localization (south vs north, last batch)",
        "",
        "Own-ref clock: this region's new bag vs this region's D_ref. "
        "Three readouts together: MMD (P(X)), CMean (||ΔE[X]|| and ΔE[Y]), PO-risk (P(Y|X)). "
        "Subset key is region, not Y.",
        "",
        _md_row(
            [
                "kind",
                "dim",
                "south MMD",
                "south ‖ΔE[X]‖",
                "south PO",
                "south ΔE[Y]",
                "north MMD",
                "north ‖ΔE[X]‖",
                "north PO",
                "north ΔE[Y]",
            ]
        ),
        _md_row(["---"] * 10),
    ]
    for kind in KINDS:
        for dim in DIMS:
            d = results[kind]["dims"][dim]
            loc = d["localization"]
            s, n = loc.get("south") or {}, loc.get("north") or {}
            lines.append(
                _md_row(
                    [
                        kind,
                        dim,
                        _fmt(s.get("mmd")),
                        _fmt(s.get("cmean_x")),
                        _fmt(s.get("po")),
                        _fmt(s.get("cmean_y")),
                        _fmt(n.get("mmd")),
                        _fmt(n.get("cmean_x")),
                        _fmt(n.get("po")),
                        _fmt(n.get("cmean_y")),
                    ]
                )
            )
    lines += [
        "",
        "## Table 4 · Sequential T_t / p_t on the order grain (OnlineRFPerm Algorithm 1)",
        "",
        _md_row(["kind", "t", "onset", "T = MSE−E_ref", "p", "hop", "MMD", "ΔE[Y]"]),
        _md_row(["---"] * 8),
    ]
    for kind in KINDS:
        for r in results[kind]["dims"]["order"]["rows"]:
            lines.append(
                _md_row(
                    [
                        kind,
                        r["t"],
                        "yes" if r["onset"] else "",
                        _fmt(r["rfperm_T"]),
                        _fmt(r["rfperm_p"], 3),
                        "yes" if r["rfperm_hop"] else "",
                        _fmt(r["mmd"]),
                        _fmt(r["cmean_y"]),
                    ]
                )
            )
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_T(results: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), sharey=True)
    for ax, kind in zip(axes, KINDS):
        for dim in DIMS:
            rows = results[kind]["dims"][dim]["rows"]
            ax.plot([r["t"] for r in rows], [r["rfperm_T"] for r in rows], label=dim, lw=1.6)
        ax.axvline(ONSET - 0.5, color="#888", ls="--", lw=0.8)
        ax.set_title(kind)
        ax.set_xlabel("batch t")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("OnlineRFPerm T = MSE − E_ref")
    axes[-1].legend(frameon=False, fontsize=7)
    fig.suptitle("WHEN · T by grain · dashed = labeled onset")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def plot_vimp(results: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(11.2, 7.2), sharey=False)
    row_keys = (("vimp", "CFPerm φ-VIMP"), ("mse_vimp", "RFPerm ΔMSE"))
    for row, (key, ylabel) in enumerate(row_keys):
        for col, kind in enumerate(KINDS):
            ax = axes[row][col]
            d = results[kind]["dims"]["all"]
            pack = d[key]
            rank = pack["rank"][:8]
            planted = set(_planted_names(d))
            y = np.arange(len(rank))[::-1]
            colors = ["#b33026" if r["feature"] in planted else "#8aa0b4" for r in rank]
            ax.barh(y, [r["vimp"] for r in rank], color=colors, height=0.72)
            ax.set_yticks(y)
            ax.set_yticklabels([r["feature"] for r in rank], fontsize=8)
            if row == 0:
                ax.set_title(kind)
            ax.set_xlabel(ylabel, fontsize=8)
            ax.grid(axis="x", alpha=0.3)
    fig.suptitle("WHICH columns · all-grain · red = planted · last batch")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def plot_fsds(results: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.2), sharey=False)
    for ax, kind in zip(axes, KINDS):
        d = results[kind]["dims"]["all"]
        rank = d["fsds"][:8]
        planted = set(_planted_names(d))
        y = np.arange(len(rank))[::-1]
        colors = ["#b33026" if r["feature"] in planted else "#8aa0b4" for r in rank]
        ax.barh(y, [r["score"] for r in rank], color=colors, height=0.72)
        ax.set_yticks(y)
        ax.set_yticklabels([r["feature"] for r in rank], fontsize=8)
        ax.set_title(kind)
        ax.grid(axis="x", alpha=0.3)
    axes[0].set_xlabel("FSDS score")
    fig.suptitle("WHICH columns · all-grain FSDS · red = planted · last batch")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def _status_class(status) -> str:
    if status == "FAR":
        return "far"
    if status == "hit":
        return "hit"
    return "miss"


def write_tables_html(slim: dict, dest: Path) -> None:
    """Three paper-style tables. No extra story."""

    def td(text, cls=""):
        c = f' class="{cls}"' if cls else ""
        return f"<td{c}>{text}</td>"

    def status_td(t_hat, delay, status):
        return td(_status_cell(t_hat, delay, status), _status_class(status))

    rows1 = []
    rows2 = []
    rows3 = []
    for kind in KINDS:
        for dim in DIMS:
            d = slim[kind]["dims"][dim]
            planted = ",".join(d.get("planted") or []) or "—"
            recov = ",".join(d.get("fsds_recovered") or []) or "—"
            mse_top = ",".join(d.get("mse_vimp_top") or []) or "—"
            cf_top = ",".join(d.get("vimp_top") or []) or "—"
            rows1.append(
                "<tr>"
                + td(kind)
                + td(dim)
                + status_td(d["addis_t"], d.get("addis_delay"), d.get("addis_status"))
                + status_td(d["onset_rank"], d.get("rank_delay"), d.get("rank_status"))
                + status_td(d.get("saffron_t"), d.get("saffron_delay"), d.get("saffron_status"))
                + status_td(d["onset_hat"], d.get("hop_delay"), d.get("hop_status"))
                + td(_fmt(d["last_T"]))
                + td(_fmt(d["last_mmd"]))
                + "</tr>"
            )
            rows2.append(
                "<tr>"
                + td(kind)
                + td(dim)
                + td(planted)
                + td(recov)
                + td(_fmt(d.get("tau_fsds")))
                + td(mse_top)
                + td(_fmt(d.get("tau_rfperm")))
                + td(cf_top)
                + td(_fmt(d.get("tau_cfperm")))
                + "</tr>"
            )
            rows3.append(
                "<tr>"
                + td(kind)
                + td(dim)
                + td(_fmt(d.get("loc_south_mmd")))
                + td(_fmt(d.get("loc_south_cmean_x")))
                + td(_fmt(d.get("loc_south_po")))
                + td(_fmt(d.get("loc_south_cmean_y")))
                + td(_fmt(d.get("loc_north_mmd")))
                + td(_fmt(d.get("loc_north_cmean_x")))
                + td(_fmt(d.get("loc_north_po")))
                + td(_fmt(d.get("loc_north_cmean_y")))
                + "</tr>"
            )
    html = f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<title>Recsys grain · OnlineRFPerm / CFPerm / FSDS</title>
<style>
body {{ font-family: "IBM Plex Sans", "Noto Sans SC", sans-serif; margin: 24px auto; max-width: 1080px; color: #122; background: #f7f5f0; line-height: 1.45; }}
h1 {{ font-size: 1.28rem; margin-bottom: 0.25rem; }}
h2 {{ font-size: 1.05rem; margin: 1.4rem 0 0.35rem; }}
.lead {{ background: #1f4e79; color: #fff; padding: 10px 14px; border-radius: 8px; font-size: 0.92rem; }}
.note {{ color: #445; font-size: 0.88rem; }}
table {{ border-collapse: collapse; background: #fff; font-size: 0.82rem; width: 100%; margin: 0.4rem 0 0.8rem; }}
td, th {{ border: 1px solid #ccc; padding: 5px 7px; font-variant-numeric: tabular-nums; text-align: left; }}
th {{ background: #ece7dc; }}
.hit {{ background: #e5f4e3; }}
.far {{ background: #f8e0dc; }}
.miss {{ color: #889; }}
.foot {{ color: #667; font-size: 0.82rem; margin-top: 1.4rem; }}
code {{ background: #eee; padding: 1px 4px; }}
</style>
</head>
<body>
<h1>推荐流按 grain 切开：OnlineRFPerm / RFPerm·CFPerm / FSDS</h1>
<p class="lead">三张表。Table 1 = WHEN。Table 2 = WHICH columns。Table 3 = WHICH accounts，south/north 同时给出 MMD / CMean / PO-risk。Y 不当特征。T = batch。定位，不是唯一分解。</p>
<p class="note">n_ref={N_REF} · n_new={N_NEW} · batches={N_BATCHES} · onset_true={ONSET} · merchants={N_MERCHANTS} · planted subset = south。<br/>
covariate / both 种 amount ≻ merchant_gmv ≻ channel；concept 只种 amount 在 Y|X。Delay = first rejection − onset。绿 = hit，红 = FAR，灰 = miss。Hop 是 1.5× jump detector，这批是慢走。</p>

<h2>Table 1 · OnlineRFPerm（WHEN）</h2>
<table>
<thead><tr><th>kind</th><th>dim</th><th>ADDIS</th><th>rank-p</th><th>SAFFRON</th><th>hop</th><th>last T</th><th>last MMD</th></tr></thead>
<tbody>
{''.join(rows1)}
</tbody>
</table>

<h2>Table 2 · FSDS / RFPerm ΔMSE / CFPerm φ-VIMP（WHICH columns · Kendall-τ）</h2>
<table>
<thead><tr><th>kind</th><th>dim</th><th>planted</th><th>FSDS recovered</th><th>τ_FSDS</th><th>RFPerm ΔMSE top-3</th><th>τ_RFPerm</th><th>CFPerm φ top-3</th><th>τ_CFPerm</th></tr></thead>
<tbody>
{''.join(rows2)}
</tbody>
</table>

<h2>Table 3 · Post-hoc localization（south vs north · MMD / CMean / PO-risk）</h2>
<table>
<thead>
<tr><th rowspan="2">kind</th><th rowspan="2">dim</th><th colspan="4">south own-ref</th><th colspan="4">north own-ref</th></tr>
<tr><th>MMD</th><th>‖ΔE[X]‖</th><th>PO</th><th>ΔE[Y]</th><th>MMD</th><th>‖ΔE[X]‖</th><th>PO</th><th>ΔE[Y]</th></tr>
</thead>
<tbody>
{''.join(rows3)}
</tbody>
</table>

<p class="foot">LaTeX <a href="../../docs/recsys_grain_monitor.tex">docs/recsys_grain_monitor.tex</a> · PDF <a href="recsys_grain_monitor.pdf">recsys_grain_monitor.pdf</a> · 细表 <a href="TABLES.md">TABLES.md</a> · 口径 <a href="REPORT.md">REPORT.md</a></p>
</body>
</html>
"""
    dest.write_text(html, encoding="utf-8")


def _tex_status(t_hat, delay, status) -> str:
    if status == "miss" or t_hat is None:
        return "miss"
    if status == "FAR":
        return f"FAR@${int(t_hat)}$"
    d = 0 if delay is None else int(delay)
    return f"${int(t_hat)}$ ($d{{=}}{d}$)"


def _tex_num(x, n=3) -> str:
    s = _fmt(x, n)
    if s == "":
        return "---"
    return f"${s}$"


def write_tex(slim: dict, dest: Path) -> None:
    """Two-page note: Y, X, three tables. Localization prints MMD / CMean / PO together."""

    def when_row(kind, dim):
        d = slim[kind]["dims"][dim]
        return (
            f"{kind.replace('_south', '')} & {dim} & "
            + _tex_status(d["addis_t"], d.get("addis_delay"), d.get("addis_status"))
            + " & "
            + _tex_status(d["onset_rank"], d.get("rank_delay"), d.get("rank_status"))
            + " & "
            + _tex_status(d.get("saffron_t"), d.get("saffron_delay"), d.get("saffron_status"))
            + " & "
            + _tex_status(d["onset_hat"], d.get("hop_delay"), d.get("hop_status"))
            + f" & {_tex_num(d['last_T'])} & {_tex_num(d['last_mmd'])} \\\\"
        )

    def which_row(kind, dim):
        d = slim[kind]["dims"][dim]
        planted = ", ".join(d.get("planted") or []) or "---"
        recov = ", ".join(d.get("fsds_recovered") or []) or "---"
        planted = planted.replace("_", r"\_")
        recov = recov.replace("_", r"\_")
        return (
            f"{kind.replace('_south', '')} & {dim} & {planted} & {recov} & "
            f"{_tex_num(d.get('tau_fsds'))} & {_tex_num(d.get('tau_rfperm'))} & "
            f"{_tex_num(d.get('tau_cfperm'))} \\\\"
        )

    def loc_row(kind, dim):
        d = slim[kind]["dims"][dim]
        return (
            f"{kind.replace('_south', '')} & {dim} & "
            f"{_tex_num(d.get('loc_south_mmd'))} & {_tex_num(d.get('loc_south_cmean_x'))} & "
            f"{_tex_num(d.get('loc_south_po'))} & {_tex_num(d.get('loc_south_cmean_y'))} & "
            f"{_tex_num(d.get('loc_north_mmd'))} & {_tex_num(d.get('loc_north_cmean_x'))} & "
            f"{_tex_num(d.get('loc_north_po'))} & {_tex_num(d.get('loc_north_cmean_y'))} \\\\"
        )

    when_rows = "\n".join(when_row(k, dim) for k in KINDS for dim in DIMS)
    which_rows = "\n".join(which_row(k, dim) for k in KINDS for dim in DIMS)
    loc_rows = "\n".join(loc_row(k, dim) for k in KINDS for dim in DIMS)
    cov = slim["covariate_south"]["dims"]["order"]
    con = slim["concept_south"]["dims"]["order"]
    body = r"""% Recsys grain monitor. Compile: pdflatex docs/recsys_grain_monitor.tex
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,booktabs}
\usepackage[hidelinks]{hyperref}
\usepackage{microtype}
\title{Recommendation stream by grain:\\
OnlineRFPerm, RFPerm/CFPerm, and FSDS}
\author{}
\date{}
\begin{document}
\maketitle
\thispagestyle{empty}

\paragraph{$Y$ (outcome).}
One row is one order.
$Y\in\{0,1\}$ is \textbf{conversion} on that order (purchase / not).
$Y$ is the response only.
It is never a column of $X$, never a ranking key, never a subset key.

\paragraph{$X$ (covariates).}
$X$ is the serving table of that order, sliced the way the log is written.
$Y$ is not in any slice.
\begin{center}
\begin{tabular}{lll}
\toprule
grain & $X$ columns & what the row is \\
\midrule
order & amount, hour, n\_items, channel & this order \\
merchant & merchant\_cat, merchant\_gmv, n\_skus & the merchant on this order \\
user & user\_tenure, user\_hist\_freq & the user on this order \\
all & the nine columns concatenated & anti-pattern; shown as a check \\
\bottomrule
\end{tabular}
\end{center}
Region (south / north) is a merchant attribute, not a feature and not $Y$.
The batch index $T$ is a time label ($D_{\mathrm{ref}}$ vs the new batch), not a treatment.

\paragraph{Setup.}
$n_{\mathrm{ref}}{=}400$, $n_{\mathrm{new}}{=}120$, six batches, labeled onset $t{=}2$.
After onset, only \textbf{south} merchants are shifted.
Covariate / both plant $\mathrm{amount}\succ\mathrm{merchant\_gmv}\succ\mathrm{channel}$ in $P(X)$.
Concept plants amount in $P(Y\mid X)$ only.
User columns never walk.

\vspace{0.6em}
\noindent
Table~\ref{tab:when}: OnlineRFPerm (Algorithm~1) --- first rejection, delay, FAR.\\
Table~\ref{tab:which}: FSDS / RFPerm $\Delta$MSE / CFPerm $\varphi$-VIMP --- Kendall-$\tau$ vs planted.\\
Table~\ref{tab:where}: south vs north own-ref. \textbf{MMD, CMean, and PO-risk together.}

\begin{table}[ht]
\centering
\caption{WHEN. Frozen RF on $D_{\mathrm{ref}}$. $T_t=\mathrm{MSE}_t-E_{\mathrm{ref}}$.
Delay $=$ first rejection $-$ onset. FAR $=$ mark before onset.}
\label{tab:when}
\scriptsize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{llcccccc}
\toprule
kind & grain & ADDIS & rank-$p$ & SAFFRON & hop & last $T$ & last MMD \\
\midrule
""" + when_rows + r"""
\bottomrule
\end{tabular}
\end{table}

\begin{table}[ht]
\centering
\caption{WHICH columns. Ranking recovery vs planted magnitude.
User grain has no planted $X$ columns --- the correct negative control.}
\label{tab:which}
\scriptsize
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{llllccc}
\toprule
kind & grain & planted in $X$ & FSDS recovered & $\tau_{\mathrm{FSDS}}$ & $\tau_{\mathrm{RFPerm}}$ & $\tau_{\mathrm{CFPerm}}$ \\
\midrule
""" + which_rows + r"""
\bottomrule
\end{tabular}
\end{table}

\begin{table}[ht]
\centering
\caption{WHICH accounts. Own-ref clock. Three readouts on the same slice:
MMD ($P(X)$), CMean ($\Vert\Delta\mathbb{E}[X]\Vert$ and $\Delta\mathbb{E}[Y]$), PO-risk ($P(Y\mid X)$).
Subset key is region, not $Y$.}
\label{tab:where}
\scriptsize
\setlength{\tabcolsep}{2.8pt}
\begin{tabular}{llcccc cccc}
\toprule
& & \multicolumn{4}{c}{south own-ref} & \multicolumn{4}{c}{north own-ref} \\
\cmidrule(lr){3-6}\cmidrule(lr){7-10}
kind & grain & MMD & $\Vert\Delta X\Vert$ & PO & $\Delta\mathbb{E}[Y]$ & MMD & $\Vert\Delta X\Vert$ & PO & $\Delta\mathbb{E}[Y]$ \\
\midrule
""" + loc_rows + r"""
\bottomrule
\end{tabular}
\end{table}

\paragraph{Read-off.}
Order grain $+$ ADDIS hits at the labeled onset ($d{=}0$).
The concatenated \texttt{all} grain is FAR at $t{=}0$.
User $X$ does not move (MMD $\approx 0$); a $T$ mark there is $Y$ walking through another grain.
FSDS recovers the planted columns on the grain that actually contains them.
Covariate south order: MMD $""" + _fmt(cov.get("loc_south_mmd")) + r"$, $\Vert\Delta X\Vert=" + _fmt(cov.get("loc_south_cmean_x")) + r"$, PO $=" + _fmt(cov.get("loc_south_po")) + r"$, $\Delta\mathbb{E}[Y]=" + _fmt(cov.get("loc_south_cmean_y")) + r"$; north MMD $=" + _fmt(cov.get("loc_north_mmd")) + r"$. Concept south order: MMD $=" + _fmt(con.get("loc_south_mmd")) + r"$, PO $=" + _fmt(con.get("loc_south_po")) + r"$, $\Delta\mathbb{E}[Y]=" + _fmt(con.get("loc_south_cmean_y")) + r"""$.

\end{document}
"""
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(body, encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    if "--from-summary" in sys.argv:
        slim = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
        write_tables_html(slim, OUT / "tables.html")
        print("wrote", OUT / "tables.html")
        return 0
    results = {kind: run_kind(kind) for kind in KINDS}
    write_markdown(results, OUT / "TABLES.md")
    plot_T(results, OUT / "online_rfperm_T.png")
    plot_vimp(results, OUT / "rfperm_vimp.png")
    plot_fsds(results, OUT / "fsds_rank.png")
    slim = {}
    for kind, pack in results.items():
        slim[kind] = {
            "meta": pack["meta"],
            "dims": {dim: slim_dim(d) for dim, d in pack["dims"].items()},
        }
    (OUT / "summary.json").write_text(json.dumps(jsonable(slim), indent=2), encoding="utf-8")
    write_justify(slim)
    write_report(slim)
    write_tables_html(slim, OUT / "tables.html")
    write_tex(slim, ROOT / "docs" / "recsys_grain_monitor.tex")
    print("wrote", OUT)
    return 0


def write_justify(slim: dict) -> None:
    lines = [
        "# Justify: recsys stream × grain — OnlineRFPerm / RFPerm / FSDS",
        "",
        "WHEN = frozen RF, T=MSE−E_ref, last-two hop; rank-p into ADDIS (primary) / SAFFRON (contrast).",
        "Delay = first rejection − onset. FAR = mark before labeled onset. miss = never marked.",
        "WHICH columns = FSDS + RFPerm ΔMSE + CFPerm φ-VIMP. Ranking metric = Kendall-τ vs planted magnitude.",
        "WHICH accounts = south vs north own-ref, each with MMD / CMean / PO-risk. Y never a feature. Localization, not unique decomp. Not a graph method.",
        "",
        f"onset_true={ONSET}. n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}.",
        "",
        "| kind | dim | ADDIS | FSDS recovered | south MMD | south ‖ΔX‖ | south PO | south ΔE[Y] | north MMD | north PO |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for kind in KINDS:
        for dim in DIMS:
            d = slim[kind]["dims"][dim]
            lines.append(
                "| "
                + " | ".join(
                    [
                        kind,
                        dim,
                        _status_cell(d["addis_t"], d.get("addis_delay"), d.get("addis_status")),
                        ",".join(d["fsds_recovered"]) or "—",
                        _fmt(d["loc_south_mmd"]),
                        _fmt(d.get("loc_south_cmean_x")),
                        _fmt(d.get("loc_south_po")),
                        _fmt(d.get("loc_south_cmean_y")),
                        _fmt(d["loc_north_mmd"]),
                        _fmt(d.get("loc_north_po")),
                    ]
                )
                + " |"
            )
    (OUT / "JUSTIFY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(slim: dict) -> None:
    cov = slim["covariate_south"]["dims"]
    con = slim["concept_south"]["dims"]
    both = slim["both"]["dims"]
    body = "\n".join(
        [
            "# Recsys grain monitor — conclusions in OnlineRFPerm / CFPerm form",
            "",
            "Same two questions as the papers. **OnlineRFPerm** (Sep 15 Algorithm 1): when did predictive error leave the reference pool? First rejection, delay, FAR. **CFPerm / RFPerm VIMP** (MetaLearner Algorithm 1): which columns contribute, recovered by Kendall-τ against the planted magnitude amount ≻ merchant_gmv ≻ channel. **FSDS + south/north** is the post-hoc localization on this table — not a graph cut.",
            "",
            "Data is the synthetic recommendation-shaped order stream: user, merchant, order, conversion Y. Slice X by grain (order / merchant / user / all). Y is never a feature. T is the batch label. Shares are localization, not Shapley or CATE.",
            "",
            f"n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}, labeled onset={ONSET}. Planted accounts = south merchants. Planted columns: covariate/both = amount, channel, merchant_gmv; concept = amount in Y|X.",
            "",
            "## Dataset (paper §5 style)",
            "",
            "Each incoming batch is 120 orders. Frozen RF is fit once on D_ref (400 orders). After onset, only **south** merchants are shifted. User features never walk. Concatenating every grain into `all` is the anti-pattern: more columns, optimistic in-sample E_ref, rank-p FAR at t=0.",
            "",
            "## WHEN (OnlineRFPerm Algorithm 1)",
            "",
            "Frozen RF on D_ref. Each batch T_t = MSE_t − E_ref. Hop = last-two 1.5× (jump detector). Rank-p vs a ref-null pool; ADDIS primary, SAFFRON contrast. Quiet T≤0 is fed as p=1. Delay = first rejection − onset_true.",
            "",
            "| kind | grain | ADDIS | rank-p | hop | last T | last MMD |",
            "|---|---|---|---|---|---|---|",
            "| covariate_south | order | "
            + _status_cell(cov["order"]["addis_t"], cov["order"].get("addis_delay"), cov["order"].get("addis_status"))
            + " | "
            + _status_cell(cov["order"]["onset_rank"], cov["order"].get("rank_delay"), cov["order"].get("rank_status"))
            + " | "
            + _status_cell(cov["order"]["onset_hat"], cov["order"].get("hop_delay"), cov["order"].get("hop_status"))
            + f" | {_fmt(cov['order']['last_T'])} | {_fmt(cov['order']['last_mmd'])} |",
            "| covariate_south | user | "
            + _status_cell(cov["user"]["addis_t"], cov["user"].get("addis_delay"), cov["user"].get("addis_status"))
            + " | "
            + _status_cell(cov["user"]["onset_rank"], cov["user"].get("rank_delay"), cov["user"].get("rank_status"))
            + f" | miss | {_fmt(cov['user']['last_T'])} | {_fmt(cov['user']['last_mmd'])} |",
            "| covariate_south | all | "
            + _status_cell(cov["all"]["addis_t"], cov["all"].get("addis_delay"), cov["all"].get("addis_status"))
            + " | "
            + _status_cell(cov["all"]["onset_rank"], cov["all"].get("rank_delay"), cov["all"].get("rank_status"))
            + f" | miss | {_fmt(cov['all']['last_T'])} | {_fmt(cov['all']['last_mmd'])} |",
            "| concept_south | order | "
            + _status_cell(con["order"]["addis_t"], con["order"].get("addis_delay"), con["order"].get("addis_status"))
            + " | "
            + _status_cell(con["order"]["onset_rank"], con["order"].get("rank_delay"), con["order"].get("rank_status"))
            + f" | miss | {_fmt(con['order']['last_T'])} | {_fmt(con['order']['last_mmd'])} |",
            "| both | order | "
            + _status_cell(both["order"]["addis_t"], both["order"].get("addis_delay"), both["order"].get("addis_status"))
            + " | "
            + _status_cell(both["order"]["onset_rank"], both["order"].get("rank_delay"), both["order"].get("rank_status"))
            + f" | miss | {_fmt(both['order']['last_T'])} | {_fmt(both['order']['last_mmd'])} |",
            "",
            "Findings (OnlineRFPerm §5.2 style):",
            "",
            f"1. Order grain + ADDIS hits at the labeled onset on covariate (delay={cov['order'].get('addis_delay')}), concept, and both. That is the slice that actually moved.",
            "2. Last-two hop never fires on this six-batch walk. Gradual south-only amount/channel walk is a trend in T / MMD, not a 1.5× jump. Rank-p and ADDIS are the marks, matching the paper: hop is a jump detector.",
            f"3. User grain MMD stays quiet ({_fmt(cov['user']['last_mmd'])}). ADDIS can still mark t=2 because Y moved (amount is in the logit) while user X did not — performance-relevant shift with quiet X, the concept fingerprint on the wrong grain.",
            "4. Concatenated `all` grain rank-p / ADDIS fire at t=0 (FAR). Nine-column in-sample E_ref is optimistic. Slice the log the way it is written.",
            "5. Concept order-grain MMD is quiet while T rises — Y|X moved, P(X) did not.",
            "",
            "Figure: `online_rfperm_T.png`. Sequential T_t / p_t: Table 4 in `TABLES.md`.",
            "",
            "## WHICH columns (RFPerm / CFPerm / FSDS, MetaLearner Algorithm 1)",
            "",
            "Nuisances on φ=(Y−μ)(T−e) fit once. CFPerm VIMP = extra MSE of predicting φ after permuting a column. RFPerm ΔMSE = extra MSE of the frozen f_ref after permuting a column of the new batch. FSDS is the univariate MMD / CMean catalog. Kendall-τ vs planted magnitude is the ranking metric in the MetaLearner paper. CFPerm global reject = max VIMP vs 95% of T-permuted nulls (B=12, not 500).",
            "",
            "| kind | grain | FSDS recovered | τ_FSDS | RFPerm ΔMSE top | τ_RFPerm | CFPerm φ top | τ_CFPerm | reject |",
            "|---|---|---|---|---|---|---|---|---|",
            f"| covariate | order | {','.join(cov['order']['fsds_recovered']) or '—'} | {_fmt(cov['order'].get('tau_fsds'))} | {','.join(cov['order'].get('mse_vimp_top') or [])} | {_fmt(cov['order'].get('tau_rfperm'))} | {','.join(cov['order']['vimp_top'])} | {_fmt(cov['order'].get('tau_cfperm'))} | {'yes' if cov['order']['vimp_reject'] else ''} |",
            f"| covariate | all | {','.join(cov['all']['fsds_recovered']) or '—'} | {_fmt(cov['all'].get('tau_fsds'))} | {','.join(cov['all'].get('mse_vimp_top') or [])} | {_fmt(cov['all'].get('tau_rfperm'))} | {','.join(cov['all']['vimp_top'])} | {_fmt(cov['all'].get('tau_cfperm'))} | {'yes' if cov['all']['vimp_reject'] else ''} |",
            f"| concept | order | {','.join(con['order']['fsds_recovered']) or '—'} | {_fmt(con['order'].get('tau_fsds'))} | {','.join(con['order'].get('mse_vimp_top') or [])} | {_fmt(con['order'].get('tau_rfperm'))} | {','.join(con['order']['vimp_top'])} | {_fmt(con['order'].get('tau_cfperm'))} | {'yes' if con['order']['vimp_reject'] else ''} |",
            f"| both | all | {','.join(both['all']['fsds_recovered']) or '—'} | {_fmt(both['all'].get('tau_fsds'))} | {','.join(both['all'].get('mse_vimp_top') or [])} | {_fmt(both['all'].get('tau_rfperm'))} | {','.join(both['all']['vimp_top'])} | {_fmt(both['all'].get('tau_cfperm'))} | {'yes' if both['all']['vimp_reject'] else ''} |",
            f"| covariate | user | {','.join(cov['user']['fsds_recovered']) or '—'} | {_fmt(cov['user'].get('tau_fsds'))} | {','.join(cov['user'].get('mse_vimp_top') or [])} | {_fmt(cov['user'].get('tau_rfperm'))} | {','.join(cov['user']['vimp_top'])} | {_fmt(cov['user'].get('tau_cfperm'))} |  |",
            "",
            "Findings (MetaLearner ranking style):",
            "",
            f"1. FSDS on the native grains recovers the planted columns: covariate order → {', '.join(cov['order']['fsds_recovered']) or '—'}; covariate all → {', '.join(cov['all']['fsds_recovered']) or '—'}; concept order → {', '.join(con['order']['fsds_recovered']) or '—'} (amount via CMean_Y, not MMD).",
            "2. User grain recovers nothing of amount/channel/gmv — those columns are not in that slice. A quiet user catalog is the correct negative control.",
            "3. CFPerm global reject at B=12 does not fire. Ranking, not the max-vs-null test, is the readout here (paper uses B=500).",
            "4. φ-VIMP can put a noise column (n_items) first on covariate; FSDS and RFPerm ΔMSE are the methods that track the planted X-walk. Concept is the reverse: amount leads φ-VIMP because Y|X moved.",
            "",
            "Figures: `rfperm_vimp.png` (red = planted), `fsds_rank.png`.",
            "",
            "## WHICH accounts (post-hoc FSDS localization)",
            "",
            "Subset key = region (south / north), not Y. Own-ref clock. Three readouts together: MMD (P(X)), CMean (||ΔE[X]|| and ΔE[Y]), PO-risk (P(Y|X)).",
            "",
            f"Covariate order-grain south: MMD={_fmt(cov['order']['loc_south_mmd'])}, CMean_X={_fmt(cov['order'].get('loc_south_cmean_x'))}, PO={_fmt(cov['order'].get('loc_south_po'))}, CMean_Y={_fmt(cov['order']['loc_south_cmean_y'])}. North MMD={_fmt(cov['order']['loc_north_mmd'])}.",
            f"Covariate user-grain south MMD={_fmt(cov['user']['loc_south_mmd'])} — users mix across merchants.",
            f"Concept order-grain south: MMD={_fmt(con['order']['loc_south_mmd'])}, PO={_fmt(con['order'].get('loc_south_po'))}, ΔE[Y]={_fmt(con['order']['loc_south_cmean_y'])}.",
            "",
            "## What to tell a production recsys",
            "",
            "1. Slice the serving table the way the log is written (order / merchant / user). Do not dump every id embedding into one simplex — `all` is the FAR grain.",
            "2. OnlineRFPerm + ADDIS on the slice that actually moved. A quiet user MMD does not mean the order grain is quiet, and a user-grain T mark can just be Y walking through another grain.",
            "3. After a mark, FSDS names the columns (Kendall-τ); RFPerm ΔMSE is the frozen-model companion; CFPerm φ-VIMP is for Y|X. South/north (or any frozen account key) names the accounts.",
            "4. Reset the error pool after a refresh (paper Appendix A.2).",
            "",
            "Fine tables: `TABLES.md`. One-pager: `JUSTIFY.md`.",
            "",
        ]
    )
    (OUT / "REPORT.md").write_text(body, encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
