#!/usr/bin/env python3
"""Order-grain FSDS + two-layer subset attribution (MMD, PO-risk, Conditional Mean)."""
from __future__ import annotations

import html as html_lib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
for p in (SRC, ROOT):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from graph_fsds_localize import run_pipeline  # noqa: E402
from stream_dgps import make_order_graph_stream  # noqa: E402

OUT = ROOT / "results" / "graph_fsds_localize"
KINDS = ("covariate_south", "concept_south", "both")
N_REF = 480
N_NEW = 160
N_BATCHES = 3
ONSET = 1
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


def slim_rank(rank, k=4):
    rows = []
    for r in rank[:k]:
        rows.append(
            {
                "feature": r["feature"],
                "score": r["score"],
                "mmd": r["mmd"],
                "cmean_x": r["cmean_x"],
                "cmean_y": r["cmean_y"],
                "loud": r["loud"],
            }
        )
    return rows


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
    rows = []
    for t in range(N_BATCHES):
        use_logo = t == N_BATCHES - 1
        rec = run_pipeline(
            tables,
            t=t,
            grain="order",
            mode="localize",
            subset_by="level_set",
            top_k=4,
            seed=seed,
            min_n=16,
            with_logo=use_logo,
            with_po=True,
        )
        rec["onset"] = t >= ONSET
        rows.append(rec)
    return {"meta": tables["meta"], "rows": rows}


def portrait_map(portraits):
    return {p["subset"]: p for p in portraits}


def write_markdown(results: dict, dest: Path) -> None:
    lines = [
        "# Feature library → FSDS unify → subset scan (order grain)",
        "",
        "Default cut = **own-ref subset scan** (level set `{φ≥τ}` / coverage prefix). "
        "`merchant_id` / `user_id` lift Ŝ onto orders. No graph is built. "
        "Scan cost is a sort, O(N log N), N = n_merchants. "
        "Shares are localization proxies, not a unique decomposition. Y is never a feature. "
        "Board: `library.html`.",
        "",
        f"n_ref={N_REF}, n_new={N_NEW}, merchants={N_MERCHANTS}, batches={N_BATCHES}, onset={ONSET}. "
        "Planted region is **south** (second half of merchant ids).",
        "",
        "## Subset scan by grain (lift to orders)",
        "",
    ]
    gheader = [
        "kind",
        "t",
        "onset",
        "cut",
        "loud south_frac",
        "J(scan mer,south)",
        "J(coverage mer,south)",
        "J(mass mer,south)",
        "J(user scan,south)",
        "n_loud mer",
        "FSDS selected",
        "fingerprint",
    ]
    lines.append(_md_row(gheader))
    lines.append(_md_row(["---"] * len(gheader)))
    for kind in KINDS:
        for rec in results[kind]["rows"]:
            g = rec.get("graph") or {}
            ls = g.get("level_set") or {}
            layer2 = rec.get("loud_subset")
            fp = ""
            south_frac = ""
            for p in rec.get("portraits") or []:
                if str(p.get("subset")) == str(layer2):
                    fp = p.get("fingerprint") or ""
                    south_frac = p.get("south_frac")
                    break
            mer = ls.get("merchant") or {}
            lines.append(
                _md_row(
                    [
                        kind,
                        rec["t"],
                        "yes" if rec["onset"] else "",
                        g.get("cut") or "",
                        _fmt(south_frac),
                        _fmt(ls.get("jaccard_merchant_vs_south")),
                        _fmt(ls.get("jaccard_coverage_merchant_vs_south")),
                        _fmt(ls.get("jaccard_mass_merchant_vs_south")),
                        _fmt(ls.get("jaccard_user_vs_south")),
                        mer.get("n_loud_nodes"),
                        ",".join(rec["fsds"]["selected_names"]),
                        fp,
                    ]
                )
            )
    lines += [
        "",
        "## Feature library (native-grain catalog · last batch)",
        "",
        _md_row(
            [
                "kind",
                "grain",
                "feature",
                "role",
                "planted",
                "selected",
                "loud",
                "score",
                "mmd",
                "cmean_x",
                "cmean_y",
            ]
        ),
        _md_row(["---"] * 11),
    ]
    for kind in KINDS:
        last = results[kind]["rows"][-1]
        for r in (last.get("library") or {}).get("rows") or []:
            lines.append(
                _md_row(
                    [
                        kind,
                        r.get("grain"),
                        r.get("feature"),
                        r.get("role") or "",
                        r.get("planted_how") or "",
                        "yes" if r.get("selected") else "",
                        "yes" if r.get("loud") else "",
                        _fmt(r.get("score")),
                        _fmt(r.get("mmd")),
                        _fmt(r.get("cmean_x")),
                        _fmt(r.get("cmean_y")),
                    ]
                )
            )
    lines += [
        "",
        "## Own-ref vs full-ref (node clock) and multi-layer lift-to-order",
        "",
        _md_row(
            [
                "kind",
                "t",
                "south own MMD",
                "north own MMD",
                "south own ‖ΔX‖",
                "north own ‖ΔX‖",
                "J(scan mer,south)",
                "J(coverage mer,south)",
                "J(user scan,south)",
                "MMD-slice south_frac",
                "ΔY-slice south_frac",
            ]
        ),
        _md_row(["---"] * 11),
    ]
    for kind in KINDS:
        for rec in results[kind]["rows"]:
            if not rec.get("onset"):
                continue
            g = rec.get("graph") or {}
            o = g.get("own_vs_full") or {}
            ls = g.get("level_set") or {}
            mer = ls.get("merchant") or {}
            slices = mer.get("slices") or {}
            lines.append(
                _md_row(
                    [
                        kind,
                        rec["t"],
                        _fmt(o.get("south_own_mmd")),
                        _fmt(o.get("north_own_mmd")),
                        _fmt(o.get("south_own_cmean_x")),
                        _fmt(o.get("north_own_cmean_x")),
                        _fmt(ls.get("jaccard_merchant_vs_south")),
                        _fmt(ls.get("jaccard_coverage_merchant_vs_south")),
                        _fmt(ls.get("jaccard_user_vs_south")),
                        _fmt((slices.get("mmd") or {}).get("south_frac")),
                        _fmt((slices.get("cmean_y") or {}).get("south_frac")),
                    ]
                )
            )
    lines += [
        "",
        "## Subset portraits (scan loud vs other · MMD + PO + CMean)",
        "",
    ]
    header = [
        "kind",
        "t",
        "subset",
        "n",
        "south_frac",
        "π_MMD",
        "π_PO",
        "π_CMean",
        "MMD",
        "PO",
        "ΔE[Y]",
        "‖ΔE[X]‖",
        "fingerprint",
    ]
    lines.append(_md_row(header))
    lines.append(_md_row(["---"] * len(header)))
    for kind in KINDS:
        for rec in results[kind]["rows"]:
            if not rec.get("onset"):
                continue
            for p in rec.get("portraits") or []:
                lines.append(
                    _md_row(
                        [
                            kind,
                            rec["t"],
                            p.get("subset"),
                            p.get("n"),
                            _fmt(p.get("south_frac")),
                            _fmt(p.get("pi_mmd")),
                            _fmt(p.get("pi_po")),
                            _fmt(p.get("pi_cmean")),
                            _fmt(p.get("mmd")),
                            _fmt(p.get("po")),
                            _fmt(p.get("cmean_y")),
                            _fmt(p.get("cmean_x")),
                            p.get("fingerprint") or "",
                        ]
                    )
                )
    lines += [
        "",
        "## FSDS ranking (last batch, top of each grain)",
        "",
        _md_row(["kind", "grain", "feature", "score", "mmd", "cmean_x", "cmean_y", "loud"]),
        _md_row(["---"] * 8),
    ]
    for kind in KINDS:
        last = results[kind]["rows"][-1]
        for grain, rank in last["fsds"]["rank"].items():
            for r in slim_rank(rank, k=3):
                lines.append(
                    _md_row(
                        [
                            kind,
                            grain,
                            r["feature"],
                            _fmt(r["score"]),
                            _fmt(r["mmd"]),
                            _fmt(r["cmean_x"]),
                            _fmt(r["cmean_y"]),
                            "yes" if r["loud"] else "",
                        ]
                    )
                )
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_portraits(results: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    metrics = [("pi_mmd", "π_MMD"), ("pi_po", "π_PO"), ("pi_cmean", "π_CMean")]
    for ax, kind in zip(axes, KINDS):
        rec = results[kind]["rows"][-1]
        portraits = rec.get("portraits") or []
        if not portraits:
            ax.set_title(kind)
            continue
        x = np.arange(len(metrics))
        w = 0.8 / max(len(portraits), 1)
        for i, p in enumerate(portraits[:4]):
            vals = [float(p.get(k) or 0.0) for k, _ in metrics]
            frac = float(p.get("south_frac") or 0.0)
            color = (0.75, 0.18, 0.15, 0.35 + 0.65 * frac)
            ax.bar(
                x + (i - 0.5 * (len(portraits[:4]) - 1)) * w,
                vals,
                w * 0.9,
                label=f"{p['subset']} s={frac:.2f}",
                color=color,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([t for _, t in metrics])
        ax.set_title(kind)
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(frameon=False, fontsize=7)
    axes[0].set_ylabel("subset share")
    fig.suptitle("Subset-scan loud vs other · last batch · redder = higher south fraction")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def _library_rows(rec) -> list:
    return list((rec.get("library") or {}).get("rows") or [])


def plot_library_heatmap(results: dict, dest: Path) -> None:
    """Last-batch FSDS catalog: features × (MMD, CMean_X, CMean_Y, score)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cols = [("mmd", "MMD"), ("cmean_x", "‖Δμ_X‖"), ("cmean_y", "|Δμ_Y|"), ("score", "score")]
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 5.2), sharey=True)
    for ax, kind in zip(axes, KINDS):
        rows = _library_rows(results[kind]["rows"][-1])
        labels = []
        raw = []
        planted = []
        selected = []
        for r in rows:
            star = "*" if r.get("planted") else " "
            labels.append(f"{star}{r['grain'][:3]}/{r['feature']}")
            raw.append([max(float(r.get(k) or 0.0), 0.0) for k, _ in cols])
            planted.append(bool(r.get("planted")))
            selected.append(bool(r.get("selected")))
        M = np.asarray(raw, dtype=float)
        floors = np.array([0.02, 0.15, 0.05, 0.8], dtype=float)
        N = M.copy()
        for j in range(N.shape[1]):
            mx = float(N[:, j].max()) if N.size else 0.0
            N[:, j] = N[:, j] / max(mx, float(floors[j]))
        im = ax.imshow(N, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([t for _, t in cols], fontsize=8)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8, fontfamily="monospace")
        for i in range(len(labels)):
            for j in range(len(cols)):
                val = raw[i][j]
                ax.text(
                    j,
                    i,
                    f"{val:.2g}",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="#111" if N[i, j] < 0.55 else "#fff",
                )
            if selected[i]:
                ax.add_patch(
                    plt.Rectangle(
                        (-0.5, i - 0.5),
                        len(cols),
                        1.0,
                        fill=False,
                        edgecolor="#1f4e79",
                        linewidth=1.4,
                    )
                )
        ax.set_title(kind, fontsize=10)
        ax.set_xlabel("relative to max(column, floor)")
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="relative")
    fig.suptitle("Feature library · last batch · * planted · box = FSDS selected · Y not in catalog")
    fig.savefig(dest, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_library_time(results: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0), sharey=True)
    for ax, kind in zip(axes, KINDS):
        recs = results[kind]["rows"]
        ts = [int(r["t"]) for r in recs]
        catalog = _library_rows(recs[-1])
        for item in catalog:
            ys = []
            for rec in recs:
                hit = next(
                    (r for r in _library_rows(rec) if r["feature"] == item["feature"]),
                    {},
                )
                ys.append(float(hit.get("score") or 0.0))
            kw = {"lw": 2.2 if item.get("planted") else 1.0, "alpha": 0.95 if item.get("planted") else 0.45}
            ax.plot(ts, ys, label=item["feature"], **kw)
        ax.set_title(kind)
        ax.set_xticks(ts)
        ax.set_xlabel("batch t")
        ax.grid(alpha=0.3)
        ax.axvline(ONSET - 0.5, color="#888", ls="--", lw=0.8)
    axes[0].set_ylabel("FSDS score")
    axes[-1].legend(frameon=False, fontsize=7, loc="upper left")
    fig.suptitle("Feature library over the stream · thick = planted · dashed = onset")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def plot_merchant_scan(results: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 5.4), sharex=False)
    for ax, kind in zip(axes, KINDS):
        rec = results[kind]["rows"][-1]
        pack = rec.get("scan_rank") or {}
        rows = list(pack.get("rows") or [])
        tau = float(pack.get("tau") or 0.0)
        if not rows:
            ax.set_title(kind)
            continue
        rows = list(reversed(rows))
        y = np.arange(len(rows))
        phis = [float(r["phi"]) for r in rows]
        colors = ["#b33026" if r.get("region") == "south" else "#8aa0b4" for r in rows]
        ax.barh(y, phis, color=colors, height=0.72, edgecolor="#222", linewidth=0.3)
        ax.axvline(tau, color="#1f4e79", ls="--", lw=1.1, label=f"τ={tau:.2g}")
        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{'●' if r.get('loud') else '○'} m{r['merchant_id']}" for r in rows],
            fontsize=7,
            fontfamily="monospace",
        )
        ax.set_title(kind)
        ax.grid(axis="x", alpha=0.3)
        ax.legend(frameon=False, fontsize=7, loc="lower right")
    axes[0].set_xlabel("own-ref φ")
    fig.suptitle("Merchant subset scan · last batch · red=south, blue=north · ● in {φ≥τ} · no network")
    fig.tight_layout()
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def _td(x) -> str:
    return html_lib.escape("" if x is None else str(x))


def write_library_html(results: dict, dest: Path) -> None:
    blocks = [
        "<!DOCTYPE html>",
        '<html lang="zh"><head><meta charset="utf-8"/>',
        "<title>特征库 · 订单流</title>",
        "<style>",
        'body { font-family: "IBM Plex Sans", "Noto Sans SC", sans-serif; margin: 24px; color: #122; background: #f7f5f0; }',
        "h1 { font-size: 1.45rem; }",
        ".note { max-width: 920px; line-height: 1.55; }",
        "figure { margin: 16px 0; }",
        "img { max-width: 100%; background: #fff; border: 1px solid #ddd; }",
        "table { border-collapse: collapse; background: #fff; font-size: 0.88rem; margin: 8px 0 22px; }",
        "td, th { border: 1px solid #ccc; padding: 5px 8px; font-variant-numeric: tabular-nums; }",
        "th { background: #ece7dc; }",
        "tr.planted { background: #f8e6d8; }",
        "tr.selected td.feat { font-weight: 650; }",
        "code { background: #eee; padding: 1px 4px; }",
        ".pill { display: inline-block; padding: 1px 7px; border-radius: 999px; background: #1f4e79; color: #fff; font-size: 0.75rem; }",
        ".quiet { color: #667; }",
        "</style></head><body>",
        "<h1>特征库 · 订单流 localization</h1>",
        '<p class="note">',
        "Subset scan 是正经的定位方法（Kulldorff / Neill LTSS）：每个实体打一口 own-ref 钟，",
        "<b>按 φ 排序取前缀</b>。扫描本身是一次排序，<code>O(N log N)</code>，",
        "N = 商户数（这里 16），不是边数，也不是 2<sup>N</sup>。",
        "贵的是打钟（每袋一次 MMD），那是任何 subset localization 都要付的。",
        "<b>这张看板是特征目录，不是网络图。</b> merchant_id / user_id 是外键，用来 lift，不建 graph。",
        "Y 不当特征。* / 橙色行 = 种下的列；蓝框 / 粗体 = FSDS 选中。",
        "</p>",
        "<figure><img src=\"library_heatmap.png\" alt=\"library heatmap\"><figcaption>最后一批：目录列 × 三支读数（列内 max-norm）</figcaption></figure>",
        "<figure><img src=\"library_time.png\" alt=\"library over time\"><figcaption>FSDS score 随 batch；粗线 = 种下的列</figcaption></figure>",
        "<figure><img src=\"merchant_scan.png\" alt=\"merchant scan\"><figcaption>商户粒 subset scan：红 = 南区。● 进 {φ≥τ}。没有画边。</figcaption></figure>",
    ]
    for kind in KINDS:
        last = results[kind]["rows"][-1]
        lib = last.get("library") or {}
        selected = ", ".join(lib.get("selected") or []) or "—"
        blocks += [
            f"<h2>{html_lib.escape(kind)}</h2>",
            f'<p>FSDS 选中 <span class="pill">{html_lib.escape(selected)}</span>',
            f' · loud subset = {html_lib.escape(str(last.get("loud_subset") or ""))}</p>',
            "<table><thead><tr>",
            "<th>grain</th><th>feature</th><th>role</th><th>planted</th><th>selected</th>",
            "<th>loud</th><th>score</th><th>MMD</th><th>‖Δμ_X‖</th><th>|Δμ_Y|</th>",
            "</tr></thead><tbody>",
        ]
        for r in lib.get("rows") or []:
            cls = []
            if r.get("planted"):
                cls.append("planted")
            if r.get("selected"):
                cls.append("selected")
            attr = f' class="{" ".join(cls)}"' if cls else ""
            blocks.append(
                "<tr{attr}><td>{grain}</td><td class=\"feat\">{feat}</td><td>{role}</td>"
                "<td>{how}</td><td>{sel}</td><td>{loud}</td>"
                "<td>{score}</td><td>{mmd}</td><td>{cx}</td><td>{cy}</td></tr>".format(
                    attr=attr,
                    grain=_td(r.get("grain")),
                    feat=_td(r.get("feature")),
                    role=_td(r.get("role")),
                    how=_td(r.get("planted_how") or ""),
                    sel="yes" if r.get("selected") else "",
                    loud="yes" if r.get("loud") else "",
                    score=_td(_fmt(r.get("score"))),
                    mmd=_td(_fmt(r.get("mmd"))),
                    cx=_td(_fmt(r.get("cmean_x"))),
                    cy=_td(_fmt(r.get("cmean_y"))),
                )
            )
        blocks.append("</tbody></table>")
        scan = last.get("scan_rank") or {}
        if scan.get("rows"):
            blocks.append(
                f'<p class="quiet">merchant scan τ={_fmt(scan.get("tau"))} · '
                f'n_loud={scan.get("n_loud")} / {scan.get("n")} · '
                f'{html_lib.escape(str(scan.get("read") or ""))}</p>'
            )
    blocks += ["</body></html>"]
    dest.write_text("\n".join(blocks) + "\n", encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results = {kind: run_kind(kind) for kind in KINDS}
    plot_portraits(results, OUT / "subset_shares.png")
    plot_library_heatmap(results, OUT / "library_heatmap.png")
    plot_library_time(results, OUT / "library_time.png")
    plot_merchant_scan(results, OUT / "merchant_scan.png")
    write_markdown(results, OUT / "TABLES.md")
    write_library_html(results, OUT / "library.html")
    slim = {}
    for kind, r in results.items():
        slim[kind] = {
            "meta": r["meta"],
            "rows": [
                {
                    "t": rec["t"],
                    "onset": rec["onset"],
                    "selected": rec["fsds"]["selected_names"],
                    "loud_subset": rec["loud_subset"],
                    "portraits": rec["portraits"],
                    "library": rec.get("library"),
                    "scan_rank": rec.get("scan_rank"),
                    "logo_full": rec.get("logo_full"),
                    "logo_subset": rec.get("logo_subset"),
                    "graph": rec.get("graph"),
                    "leakage": rec["leakage"],
                    "fsds_rank": {
                        g: slim_rank(rank, k=4)
                        for g, rank in rec["fsds"]["rank"].items()
                    },
                }
                for rec in r["rows"]
            ],
        }
    (OUT / "summary.json").write_text(json.dumps(jsonable(slim), indent=2), encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
