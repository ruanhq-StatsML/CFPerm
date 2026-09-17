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
    return {
        "dim": d["dim"],
        "onset_true": d["onset_true"],
        "onset_hat": d["onset_hat"],
        "onset_rank": d["onset_rank"],
        "addis_t": d["addis_t"],
        "saffron_t": d["saffron_t"],
        "n_hop": d["n_hop"],
        "fsds_top": [r["feature"] for r in d["fsds"][:4]],
        "fsds_loud": [r["feature"] for r in d["fsds"] if r.get("loud")],
        "fsds_recovered": d["fsds_recovered"],
        "vimp_top": d["vimp_top"],
        "vimp_reject": d["vimp_reject"],
        "vimp_max": d["vimp"]["max_vimp"],
        "planted": d["planted"],
        "y_in_X": d["y_in_X"],
        "loc_n_south": d["localization"]["n_south"],
        "loc_pair_mmd": d["localization"].get("pair_mmd"),
        "loc_south_mmd": (d["localization"].get("south") or {}).get("mmd"),
        "loc_north_mmd": (d["localization"].get("north") or {}).get("mmd"),
        "loc_south_cmean_y": (d["localization"].get("south") or {}).get("cmean_y"),
        "loc_north_cmean_y": (d["localization"].get("north") or {}).get("cmean_y"),
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
    }


def write_markdown(results: dict, dest: Path) -> None:
    lines = [
        "# Recsys order-stream by grain — OnlineRFPerm × RFPerm/CFPerm × FSDS",
        "",
        "Data: synthetic order / merchant / user stream (recommendation-shaped table). "
        "Y = conversion, never a feature. T = batch label. "
        "OnlineRFPerm marks WHEN (Algorithm 1). RFPerm/CFPerm VIMP and FSDS mark WHICH columns. "
        "Post-hoc region split marks WHICH accounts (south vs north). "
        "Localization, not a unique decomposition. Not a graph method.",
        "",
        f"n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}, onset_true={ONSET}, "
        f"merchants={N_MERCHANTS}. Planted subset = south.",
        "",
        "## Table 1 · OnlineRFPerm (WHEN)",
        "",
        _md_row(
            [
                "kind",
                "dim",
                "onset_true",
                "onset_hat (hop)",
                "onset_rank (p<0.05)",
                "ADDIS first rej",
                "SAFFRON first rej",
                "n_hop",
                "last T",
                "last MMD",
            ]
        ),
        _md_row(["---"] * 10),
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
                        "" if d["onset_hat"] is None else d["onset_hat"],
                        "" if d["onset_rank"] is None else d["onset_rank"],
                        "" if d["addis_t"] is None else d["addis_t"],
                        "" if d["saffron_t"] is None else d["saffron_t"],
                        d["n_hop"],
                        _fmt(d["rows"][-1]["rfperm_T"]),
                        _fmt(d["rows"][-1]["mmd"]),
                    ]
                )
            )
    lines += [
        "",
        "## Table 2 · RFPerm/CFPerm VIMP and FSDS (WHICH columns)",
        "",
        _md_row(
            [
                "kind",
                "dim",
                "planted",
                "FSDS loud",
                "FSDS recovered",
                "VIMP top-3",
                "CFPerm reject",
            ]
        ),
        _md_row(["---"] * 7),
    ]
    for kind in KINDS:
        for dim in DIMS:
            d = results[kind]["dims"][dim]
            planted = ",".join(d["planted"]) or "—"
            loud = ",".join(r["feature"] for r in d["fsds"] if r.get("loud")) or "—"
            recov = ",".join(d["fsds_recovered"]) or "—"
            lines.append(
                _md_row(
                    [
                        kind,
                        dim,
                        planted,
                        loud,
                        recov,
                        ",".join(d["vimp_top"]),
                        "yes" if d["vimp_reject"] else "",
                    ]
                )
            )
    lines += [
        "",
        "## Table 3 · Post-hoc localization (south vs north, last batch)",
        "",
        _md_row(
            [
                "kind",
                "dim",
                "n_south",
                "south MMD vs own-ref",
                "north MMD vs own-ref",
                "pair MMD",
                "south ΔE[Y]",
                "north ΔE[Y]",
            ]
        ),
        _md_row(["---"] * 8),
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
                        loc.get("n_south"),
                        _fmt(s.get("mmd")),
                        _fmt(n.get("mmd")),
                        _fmt(loc.get("pair_mmd")),
                        _fmt(s.get("cmean_y")),
                        _fmt(n.get("cmean_y")),
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

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.2), sharey=False)
    for ax, kind in zip(axes, KINDS):
        d = results[kind]["dims"]["all"]
        rank = d["vimp"]["rank"][:8]
        planted = set(d["planted"])
        y = np.arange(len(rank))[::-1]
        colors = ["#b33026" if r["feature"] in planted else "#8aa0b4" for r in rank]
        ax.barh(y, [r["vimp"] for r in rank], color=colors, height=0.72)
        ax.set_yticks(y)
        ax.set_yticklabels([r["feature"] for r in rank], fontsize=8)
        ax.set_title(kind)
        ax.grid(axis="x", alpha=0.3)
    axes[0].set_xlabel("RFPerm VIMP on φ")
    fig.suptitle("WHICH columns · all-grain VIMP · red = planted · last batch")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def write_reports(results: dict) -> None:
    # filled after slim exists; JUSTIFY/REPORT written in main from slim numbers
    pass


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = {kind: run_kind(kind) for kind in KINDS}
    write_markdown(results, OUT / "TABLES.md")
    plot_T(results, OUT / "online_rfperm_T.png")
    plot_vimp(results, OUT / "rfperm_vimp.png")
    slim = {}
    for kind, pack in results.items():
        slim[kind] = {
            "meta": pack["meta"],
            "dims": {dim: slim_dim(d) for dim, d in pack["dims"].items()},
        }
    (OUT / "summary.json").write_text(json.dumps(jsonable(slim), indent=2), encoding="utf-8")
    write_justify(slim)
    write_report(slim)
    print("wrote", OUT)
    return 0


def write_justify(slim: dict) -> None:
    lines = [
        "# Justify: recsys stream × grain — OnlineRFPerm / RFPerm / FSDS",
        "",
        "WHEN = frozen RF, T=MSE−E_ref, last-two hop; rank-p into ADDIS (primary) / SAFFRON (contrast).",
        "WHICH columns = FSDS + RFPerm VIMP on frozen φ (T=batch). WHICH accounts = south vs north.",
        "Y never a feature. Localization, not unique decomp. Not a graph method.",
        "",
        f"onset_true={ONSET}. n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}.",
        "",
        "| kind | dim | hop | rank-p | ADDIS | FSDS recovered | VIMP top | CFPerm | south MMD | north MMD |",
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
                        "" if d["onset_hat"] is None else str(d["onset_hat"]),
                        "" if d["onset_rank"] is None else str(d["onset_rank"]),
                        "" if d["addis_t"] is None else str(d["addis_t"]),
                        ",".join(d["fsds_recovered"]) or "—",
                        ",".join(d["vimp_top"][:3]),
                        "yes" if d["vimp_reject"] else "",
                        _fmt(d["loc_south_mmd"]),
                        _fmt(d["loc_north_mmd"]),
                    ]
                )
                + " |"
            )
    (OUT / "JUSTIFY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(slim: dict) -> None:
    cov = slim["covariate_south"]["dims"]
    con = slim["concept_south"]["dims"]
    body = f"""# Recsys grain monitor — conclusions in OnlineRFPerm / CFPerm form

Same two questions as the papers. **OnlineRFPerm** (Sep 15 Algorithm 1): when did predictive error leave the reference pool? **CFPerm / RFPerm VIMP** (MetaLearner Algorithm 1): which columns contribute to the PO-risk between \(D_{{\\mathrm{{ref}}}}\) and the new batch? **FSDS + south/north** is the post-hoc localization on this table — not a graph cut.

Data is the synthetic recommendation-shaped order stream: user, merchant, order, conversion Y. Slice X by grain. Y is never a feature. T is the batch label. Shares are localization, not Shapley or CATE.

n_ref={N_REF}, n_new={N_NEW}, batches={N_BATCHES}, labeled onset={ONSET}. Planted accounts = south merchants.

## WHEN (OnlineRFPerm)

Frozen RF on \(D_{{\\mathrm{{ref}}}}\). Each batch \(T_t=\\mathrm{{MSE}}_t-E_{{\\mathrm{{ref}}}}\). Hop = last-two 1.5×. Rank-p vs a ref-null pool; ADDIS primary, SAFFRON contrast. Quiet \(T\\le 0\) is fed as \(p=1\).

| kind | grain that should fire | what the stream did |
|---|---|---|
| covariate_south | order / all (amount, channel in X) | order last T={_fmt(cov["order"]["last_T"])}, MMD={_fmt(cov["order"]["last_mmd"])}; user last T={_fmt(cov["user"]["last_T"])} (user X is not planted) |
| concept_south | order via \(Y\\mid X\), MMD quiet | order MMD={_fmt(con["order"]["last_mmd"])}, ΔE[Y] south={_fmt(con["order"]["loc_south_cmean_y"])} |
| both | order / all | same X-shift as covariate, PO still to be read |

Hop / ADDIS can stay empty on a short six-batch walk — that matches the paper: last-two is a jump detector; a gradual walk shows up as T / MMD trends. Figure: `online_rfperm_T.png`.

## WHICH columns (RFPerm / FSDS)

Nuisances on φ=(Y−μ)(T−e) fit once. VIMP = extra MSE of predicting φ after permuting a column. CFPerm reject = max VIMP above the 95% quantile of T-permuted nulls. FSDS is the univariate MMD / CMean catalog on the same grain.

Covariate, all-grain FSDS recovered: {", ".join(cov["all"]["fsds_recovered"]) or "—"}. VIMP top: {", ".join(cov["all"]["vimp_top"])}.
Concept FSDS recovered: {", ".join(con["all"]["fsds_recovered"]) or "—"} (amount should lead on CMean_Y, not MMD).
User grain should not recover amount/channel — those columns are not in that slice.

Figure: `rfperm_vimp.png` (red = planted).

## WHICH accounts (post-hoc)

Subset key = region (south / north), not Y. Own-ref MMD: this region's new bag vs this region's \(D_{{\\mathrm{{ref}}}}\). Pair MMD compares the two regions inside the new batch (heterogeneity, not drift).

Covariate order-grain: south MMD={_fmt(cov["order"]["loc_south_mmd"])}, north MMD={_fmt(cov["order"]["loc_north_mmd"])}.
User grain south vs north should stay closer — users mix across merchants.

## What to tell a production recsys

1. Slice the serving table the way the log is written (order / merchant / user). Do not dump every id embedding into one simplex.
2. OnlineRFPerm on the slice that actually moved. A quiet user grain does not mean the order grain is quiet.
3. After a mark, FSDS + RFPerm VIMP name the columns; south/north (or any frozen account key) names the accounts.
4. Reset the error pool after a refresh (paper Appendix A.2).

Fine tables: `TABLES.md`. One-pager: `JUSTIFY.md`.
"""
    (OUT / "REPORT.md").write_text(body, encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
