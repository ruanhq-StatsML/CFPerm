#!/usr/bin/env python3
"""Layer-freeze board: read PO-risk to see where to start updating.

k+1 models; model_i starts updating from the top i layers.
Incoming batch is T=1, n_ref=10000. Raw PO-risk jitters on a small
n_new; a causal moving average is the stability readout. No
online-bootstrap — repeated MLP inference cannot be afforded.
MA below 2× baseline → all-layer backprop.

    PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dl_model_registry import DLModelRegistry  # noqa: E402
from layer_freeze_cv import run_deviation_gate, run_layer_freeze_cv  # noqa: E402
from streaming_po_risk import (  # noqa: E402
    MIN_STREAM_N,
    REF_N,
    annotate_moving_average,
    ma_window,
    moving_average,
)

OUT = ROOT / "results" / "layer_freeze_online_cv"
DEFAULT_STREAM = {
    "electricity": MIN_STREAM_N,
    "covertype": REF_N,
    "airlines": MIN_STREAM_N,
}


def _standardize_from_ref(X, n_ref: int):
    mu = X[:n_ref].mean(axis=0)
    sd = X[:n_ref].std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - mu) / sd


def _frame_to_xy(df, y_raw, title):
    import pandas as pd

    Y = np.asarray(y_raw).ravel()
    if Y.dtype == object or str(Y.dtype).startswith("str") or str(Y.dtype) == "category":
        ys = pd.Series(Y.astype(str))
        # delay / UP / 1 as the positive class when present
        pos = {"1", "UP", "True", "true", "Y", "yes"}
        if set(ys.unique()) <= {"0", "1"} or any(v in pos for v in ys.unique()):
            Y = ys.isin(sorted(pos)).astype(int).to_numpy()
        else:
            top = ys.value_counts().index[0]
            Y = (ys == top).astype(int).to_numpy()
    else:
        u = np.unique(Y)
        Y = (Y == (2 if 2 in u else u[np.argmax([np.sum(Y == v) for v in u])])).astype(int)
    cols = []
    blocks = []
    for c in df.columns:
        s = df[c]
        if str(s.dtype) in {"object", "category"} or s.dtype == object:
            codes, _ = pd.factorize(s.astype(str), sort=False)
            blocks.append(codes.astype(float))
        else:
            blocks.append(np.asarray(s, dtype=float))
        cols.append(str(c))
    X = np.column_stack(blocks)
    return X, Y, title, cols


def load_electricity():
    from sklearn.datasets import fetch_openml

    bunch = fetch_openml("electricity", version=1, as_frame=True, parser="auto")
    df = bunch.data.copy()
    drop = [c for c in df.columns if str(c).lower() == "date"]
    return _frame_to_xy(df.drop(columns=drop), bunch.target, "electricity (NSW, ordered in time)")


def load_covertype():
    from sklearn.datasets import fetch_covtype

    bunch = fetch_covtype()
    X = np.asarray(bunch.data, dtype=float)
    Y = (np.asarray(bunch.target) == 2).astype(int)
    return X, Y, "covertype (geographic order, class 2 vs rest)", [f"x{j}" for j in range(X.shape[1])]


def load_airlines():
    from sklearn.datasets import fetch_openml

    bunch = fetch_openml("airlines", version=1, as_frame=True, parser="auto")
    return _frame_to_xy(bunch.data.copy(), bunch.target, "airlines (flight delay, time-ordered)")


LOADERS = {
    "electricity": load_electricity,
    "covertype": load_covertype,
    "airlines": load_airlines,
}


def plot_boards(result: dict, title: str, out_dir: Path) -> list[str]:
    rows = result["rows"]
    k = result["k"]
    names = result["layer_names"]
    ts = [r["t"] for r in rows]
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    fig, ax = plt.subplots(figsize=(8.6, 4.0))
    raw = np.asarray([r["po_stream"] for r in rows], dtype=float)
    n_new = int(rows[0]["n_new"]) if rows else 1
    w = ma_window(n_new)
    ma = moving_average(raw, w) if raw.size else raw
    ax.plot(ts, raw, color="#9bb4cc", lw=1.0, marker="o", ms=3.5, label="raw")
    ax.plot(ts, ma, color="#1f4e79", lw=2.2, label=f"MA({w})")
    ax.axhline(result["po_base"], color="#888", ls="--", lw=1.6, label="ref-split baseline")
    ax.axhline(2.0 * result["po_base"], color="#b33", ls=":", lw=1.0, label="2× baseline")
    for r in rows:
        if r["large_deviation"]:
            ax.scatter([r["t"]], [r["po_stream"]], s=80, color="#b33", zorder=4)
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_ylabel("PO-risk")
    ax.set_title(title + " — MA of PO-risk (no bootstrap)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = out_dir / "po_risk_gate.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p.name)

    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    freeze_y = []
    for r in rows:
        freeze_y.append(k if r["all_trainable"] else r["i_star"])
    ax.step(ts, freeze_y, where="mid", color="#1f4e79", lw=2.2)
    ax.scatter(ts, freeze_y, color="#1f4e79", zorder=3)
    ax.set_yticks(range(k + 1), labels=names)
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_title("all trainable unless large deviation → freeze from this layer")
    ax.set_ylim(-0.4, k + 0.4)
    fig.tight_layout()
    p = out_dir / "freeze_from.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p.name)

    large_rows = [r for r in rows if r["large_deviation"] and r["layers"]]
    if large_rows:
        colors = plt.cm.tab10(np.linspace(0, 0.8, k + 1))
        fig, ax = plt.subplots(figsize=(8.6, 4.2))
        ts_l = [r["t"] for r in large_rows]
        for i in range(k + 1):
            ys = []
            ok = True
            for r in large_rows:
                hit = next((x for x in r["layers"] if x["i"] == i), None)
                if hit is None:
                    ok = False
                    break
                ys.append(hit["po_fit"])
            if ok:
                ax.plot(ts_l, ys, color=colors[i], lw=2.0, marker="o", ms=3.5, label=names[i])
        ax.set_xlabel("incoming batch (T=1), large deviation only")
        ax.set_ylabel("PO-risk (conditional on freeze-depth)")
        ax.set_title("from which layer to freeze")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        p = out_dir / "po_risk_by_layer.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        written.append(p.name)
    return written


def render_html(spec: dict, result: dict, images: list[str], out_path: Path) -> None:
    rec_i = result["recommend_i"]
    k = result["k"]
    n_large = result.get("n_large", 0)
    if n_large == 0 or rec_i == k:
        rec = "没有大 deviation，或大偏差段仍全开 → 一直 trainable"
    else:
        rec = f"有大 deviation 的段：从 model_{rec_i} 开始冻（median）"
    n_new = result["rows"][0]["n_new"] if result["rows"] else ""
    cards = []
    for img in images:
        cards.append(f'<figure><img src="{img}" alt="{img}"><figcaption>{img}</figcaption></figure>')
    table = [
        "<table><thead><tr><th>t</th><th>PO-risk</th><th>baseline</th><th>large?</th><th>action</th></tr></thead><tbody>"
    ]
    for r in result["rows"]:
        action = "all trainable" if r["all_trainable"] else f"freeze from {r['freeze_from']}"
        flag = "yes" if r["large_deviation"] else ""
        table.append(
            f"<tr><td>{r['t']}</td><td>{r['po_stream']:.4g}</td><td>{r['po_base']:.4g}</td>"
            f"<td>{flag}</td><td>{action}</td></tr>"
        )
    table.append("</tbody></table>")
    layer_tab = ""
    last_large = next((r for r in reversed(result["rows"]) if r["large_deviation"] and r["layers"]), None)
    if last_large:
        layer_tab = "<h2>Last large-deviation batch — 从哪一层冻</h2><table><thead><tr><th>model</th><th>PO-risk</th></tr></thead><tbody>"
        for x in last_large["layers"]:
            mark = " ★" if x["i"] == last_large["i_star"] else ""
            layer_tab += f"<tr><td>{x['name']}{mark}</td><td>{x['po_fit']:.4g}</td></tr>"
        layer_tab += "</tbody></table>"
    html = f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<title>PO-risk 看板</title>
<style>
body {{ font-family: "IBM Plex Sans", "Noto Sans SC", sans-serif; margin: 24px; color: #122; background: #f7f5f0; }}
h1 {{ font-size: 1.5rem; }}
.note {{ max-width: 860px; line-height: 1.5; }}
.rec {{ background: #1f4e79; color: #fff; padding: 12px 16px; border-radius: 8px; display: inline-block; }}
figure {{ margin: 18px 0; }}
img {{ max-width: 100%; background: #fff; border: 1px solid #ddd; }}
table {{ border-collapse: collapse; background: #fff; }}
td, th {{ border: 1px solid #ccc; padding: 6px 10px; font-variant-numeric: tabular-nums; }}
code {{ background: #eee; padding: 1px 4px; }}
</style>
</head>
<body>
<h1>除非有大 deviation，否则一直 trainable</h1>
<p class="note">
{spec["title"]}. n_ref={result["n_ref"]}, n_new={n_new}, hidden={result["hidden_dims"]}.
表格 PO-risk 单独维护 <b>outcome model</b> μ(Y|X) 和 <b>propensity</b> e(T|X)。
新 batch 是 <b>T=1</b>。raw PO-risk 会抖就做 causal moving average；
<b>MA 稳（不过 2× baseline）→ 全量每一层 back-propagate</b>。
不做 online-bootstrap。何时 update 是业务逻辑，不是这个数。
</p>
<p class="rec">{rec}</p>
{"".join(cards)}
<h2>每一段 — 直接读 PO-risk</h2>
{"".join(table)}
{layer_tab}
</body>
</html>
"""
    out_path.write_text(html, encoding="utf-8")


def render_report(spec: dict, result: dict) -> str:
    rec_i = result["recommend_i"]
    n_new = result["rows"][0]["n_new"] if result["rows"] else ""
    n_large = result.get("n_large", 0)
    k = result["k"]
    rec = "all trainable"
    if n_large and rec_i < k:
        rec = f"on large-deviation batches, freeze from model_{rec_i}"
    elif n_large:
        rec = "large deviation present, freeze prototype still says all trainable"
    lines = [
        f"# PO-risk board — {spec['title']}",
        "",
        "No large deviation → all-layer backprop. Large deviation → from which layer to freeze.",
        "Raw PO-risk jitters; a causal moving average is the stability readout. No online-bootstrap.",
        "Tabular PO-risk keeps a separate outcome model and a separate propensity model.",
        f"T=1 on the incoming batch. n_ref={result['n_ref']}, n_new={n_new}. baseline={result['po_base']:.4g}.",
        "",
        f"- hidden_dims = `{result['hidden_dims']}`, k = {k}",
        f"- n_batches = {result['n_batches']}, n_large = {n_large}",
        f"- **{rec}**",
        "",
        "| t | PO-risk | baseline | large | action |",
        "|---:|---:|---:|---|---|",
    ]
    for r in result["rows"]:
        action = "all trainable" if r["all_trainable"] else f"freeze from {r['freeze_from']}"
        flag = "yes" if r["large_deviation"] else ""
        lines.append(f"| {r['t']} | {r['po_stream']:.3g} | {r['po_base']:.3g} | {flag} | {action} |")
    large_rows = [r for r in result["rows"] if r["large_deviation"] and r["layers"]]
    if large_rows:
        lines += [
            "",
            "Large-deviation batches, PO-risk conditional on freeze-depth:",
            "",
            "| t | " + " | ".join(result["layer_names"]) + " | freeze from |",
            "|---:|" + "|".join(["---:"] * (k + 1)) + "|---|",
        ]
        for r in large_rows:
            by = {x["i"]: x["po_fit"] for x in r["layers"]}
            cells = " | ".join(f"{by[i]:.3g}" for i in range(k + 1))
            action = "all trainable" if r["all_trainable"] else f"freeze from {r['freeze_from']}"
            lines.append(f"| {r['t']} | {cells} | {action} |")
    lines += ["", "Read the PO-risk. Nothing else.", ""]
    return "\n".join(lines)


def run_one(name: str, n_ref: int, batch: int, max_batches: int, hidden, online_epochs: int) -> dict:
    X, Y, title, cols = LOADERS[name]()
    X = _standardize_from_ref(X, n_ref)
    registry = DLModelRegistry(
        epochs=6,
        patience=4,
        scheduler="cosine",
        optimizer="adamw",
        lr=1e-3,
        verbose=False,
        amp=False,
    )
    result = run_layer_freeze_cv(
        X,
        Y,
        n_ref=n_ref,
        batch_size_stream=batch,
        hidden_dims=tuple(hidden),
        online_epochs=online_epochs,
        max_batches=max_batches,
        registry=registry,
    )
    result["dataset"] = name
    result["title"] = title
    result["n_features"] = int(X.shape[1])
    result["n_total"] = int(len(Y))
    result["feature_names"] = list(map(str, cols))
    result["n_new"] = int(batch)
    return result


def _get_size(by_size: dict, n_new: int) -> dict:
    if n_new in by_size:
        return by_size[n_new]
    return by_size[str(n_new)]


def attach_ma(cmp_: dict) -> dict:
    """Replay-safe: annotate MA on saved gate rows. No extra inference."""
    by_size = {int(k): v for k, v in cmp_["by_size"].items()}
    for n_new, rec in by_size.items():
        info = annotate_moving_average(rec["rows"], rec["po_base"], n_new=n_new)
        rec.update(info)
    cmp_["by_size"] = by_size
    return cmp_


def plot_size_compare(by_size: dict, title: str, out_dir: Path) -> str:
    """Raw PO-risk jitters; the moving average says whether all layers can backprop."""
    sizes = sorted(int(s) for s in by_size)
    fig, axes = plt.subplots(len(sizes), 1, figsize=(8.8, 2.2 * len(sizes)), sharex=False)
    if len(sizes) == 1:
        axes = [axes]
    for ax, n_new in zip(axes, sizes):
        rec = _get_size(by_size, n_new)
        rows = rec["rows"]
        ts = [r["t"] for r in rows]
        raw = np.asarray([r["po_stream"] for r in rows], dtype=float)
        w = int(rec.get("ma_window") or ma_window(n_new))
        ma = moving_average(raw, w)
        ax.plot(ts, raw, color="#9bb4cc", lw=0.7 if len(ts) > 80 else 1.0, label="raw")
        ax.plot(ts, ma, color="#1f4e79", lw=2.0, label=f"MA({w})")
        ax.axhline(rec["po_base"], color="#888", ls="--", lw=1.2, label="baseline")
        ax.axhline(2.0 * rec["po_base"], color="#b33", ls=":", lw=1.0, label="2× baseline")
        ax.set_ylabel("PO-risk")
        stable = rec.get("all_layer_backprop")
        if stable is None:
            stable = bool(np.nanmax(ma) < 2.0 * max(rec["po_base"], 1e-12))
        tag = "MA stable → all-layer backprop" if stable else "MA hop"
        ax.set_title(f"n_new={n_new}  {tag}", fontsize=10)
        if n_new == sizes[0]:
            ax.legend(fontsize=7, ncol=4)
    axes[-1].set_xlabel("incoming batch index (T=1)")
    fig.suptitle(title + " — moving average, no bootstrap", y=1.01)
    fig.tight_layout()
    path = out_dir / "batch_size_gate.png"
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path.name


def write_batch_size_md(cmp_: dict, out_path: Path) -> str:
    lines = [
        f"# MA of PO-risk — {cmp_['title']}",
        "",
        "Raw PO-risk jitters. Causal moving average is the stability readout.",
        "No online-bootstrap (repeated MLP inference cannot be afforded).",
        "MA below 2× ref-split baseline → all-layer backprop.",
        f"n_ref={cmp_.get('n_ref', REF_N)}.",
        "",
        "| n_new | batches | MA window | raw frac large | MA frac large | MA max / baseline | all-layer backprop |",
        "|---:|---:|---:|---:|---:|---:|---|",
    ]
    for n_new in sorted(int(s) for s in cmp_["by_size"]):
        rec = _get_size(cmp_["by_size"], n_new)
        base = max(float(rec["po_base"]), 1e-12)
        ma_max = rec.get("ma_max")
        if ma_max is None:
            raw = np.asarray([r["po_stream"] for r in rec["rows"]], dtype=float)
            ma_max = float(np.nanmax(moving_average(raw, ma_window(n_new))))
        ratio = float(ma_max) / base
        backprop = "yes" if rec.get("all_layer_backprop", ratio < 2.0) else "no"
        lines.append(
            f"| {n_new} | {rec['n_batches']} | {rec.get('ma_window', ma_window(n_new))} | "
            f"{rec['frac_large']:.2f} | {rec.get('frac_ma_large', 0.0):.2f} | {ratio:.2f} | {backprop} |"
        )
    text = "\n".join(lines) + "\n"
    out_path.write_text(text, encoding="utf-8")
    return text


def write_index(out: Path) -> None:
    bits = [
        "# PO-risk board",
        "",
        "When to **update** is business logic. PO-risk does not justify that.",
        "Raw PO-risk jitters; a causal moving average is the stability readout.",
        "**MA stable (below 2× baseline) → all-layer backprop.** No online-bootstrap.",
        "Large MA hop (and the business already wants a hop) → from which layer to freeze.",
        "",
    ]
    html_items = []
    for name in ("electricity", "covertype", "airlines"):
        sub = out / name
        if (sub / "board.html").exists():
            bits.append(f"- [{name} freeze board]({name}/board.html)")
        if (sub / "batch_size_gate.png").exists():
            bits.append(f"- [{name} MA gate]({name}/batch_size_gate.png)")
        if (sub / "board.html").exists() or (sub / "batch_size_gate.png").exists():
            html_items.append(
                f'<li><a href="{name}/board.html">{name}</a> · '
                f'<a href="{name}/batch_size_gate.png">MA</a></li>'
            )
    (out / "REPORT.md").write_text("\n".join(bits) + "\n", encoding="utf-8")
    html = """<!DOCTYPE html><meta charset="utf-8"><title>PO-risk board</title>
<h1>MA 稳定 → 全层 backprop</h1>
<p>raw PO-risk 会抖。causal moving average 不过 2× baseline，就全量每一层 back-propagate。
不做 online-bootstrap。何时 update 是业务逻辑。</p>
<ul>""" + "".join(html_items) + "</ul>"
    (out / "index.html").write_text(html, encoding="utf-8")


def run_size_compare(name: str, n_ref: int, sizes: list[int], stream_cap: int) -> dict:
    X, Y, title, cols = LOADERS[name]()
    X = _standardize_from_ref(X, n_ref)
    keep = min(len(Y), n_ref + stream_cap)
    X, Y = X[:keep], Y[:keep]
    by_size = {}
    for n_new in sizes:
        max_batches = max(1, stream_cap // n_new)
        rec = run_deviation_gate(X, Y, n_ref=n_ref, batch_size_stream=n_new, max_batches=max_batches)
        by_size[int(n_new)] = rec
        print(f"  n_new={n_new} batches={rec['n_batches']} frac_large={rec['frac_large']:.2f}", flush=True)
    return {"title": title, "dataset": name, "n_ref": n_ref, "by_size": by_size}


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def replay_saved_gates(names: list[str]) -> None:
    """Redraw MA boards from on-disk JSON. No MLP, no extra inference."""
    for name in names:
        path = OUT / name / "batch_size.json"
        if not path.exists():
            print(f"skip {name}: no {path}", flush=True)
            continue
        cmp_ = attach_ma(json.loads(path.read_text(encoding="utf-8")))
        sub = OUT / name
        plot_size_compare(cmp_["by_size"], cmp_["title"], sub)
        path.write_text(json.dumps(jsonable(cmp_), indent=2) + "\n", encoding="utf-8")
        text = write_batch_size_md(cmp_, sub / "BATCH_SIZE.md")
        print(text, flush=True)
        board = sub / "summary.json"
        if board.exists():
            result = json.loads(board.read_text(encoding="utf-8"))
            spec = {"name": name, "title": result.get("title", name)}
            images = plot_boards(result, spec["title"], sub)
            render_html(spec, result, images, sub / "board.html")


def main() -> int:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="electricity", choices=sorted(LOADERS) + ["both", "all"])
    p.add_argument("--n-ref", type=int, default=REF_N)
    p.add_argument("--batch", type=int, default=0, help="incoming T=1 size; 0 = per-dataset default (≥5000)")
    p.add_argument("--max-batches", type=int, default=8)
    p.add_argument("--hidden", default="64,32")
    p.add_argument("--online-epochs", type=int, default=2)
    p.add_argument("--skip-freeze", action="store_true", help="only the batch-size gate, no freeze prototype")
    p.add_argument("--sizes", default="500,1000,2000,5000")
    p.add_argument("--stream-cap", type=int, default=30000)
    p.add_argument(
        "--replay-json",
        action="store_true",
        help="redraw MA boards from saved batch_size.json; no extra inference",
    )
    args = p.parse_args()
    hidden = tuple(int(x) for x in args.hidden.split(",") if x.strip())
    if args.dataset == "both":
        names = ["electricity", "covertype"]
    elif args.dataset == "all":
        names = ["electricity", "covertype", "airlines"]
    else:
        names = [args.dataset]
    sizes = [int(x) for x in args.sizes.split(",") if x.strip()]
    OUT.mkdir(parents=True, exist_ok=True)
    if args.replay_json:
        replay_saved_gates(names)
        write_index(OUT)
        print("wrote", OUT)
        return 0
    for name in names:
        if not args.skip_freeze:
            batch = args.batch if args.batch > 0 else DEFAULT_STREAM[name]
            if batch < MIN_STREAM_N:
                print(
                    f"note: n_new={batch} < {MIN_STREAM_N}; "
                    "raw PO-risk jitters — read the moving average, do not bootstrap",
                    flush=True,
                )
            spec = {"name": name, "title": name}
            print(f"=== {name} n_ref={args.n_ref} n_new={batch} ===", flush=True)
            result = run_one(name, args.n_ref, batch, args.max_batches, hidden, args.online_epochs)
            spec["title"] = result["title"]
            sub = OUT / name
            sub.mkdir(parents=True, exist_ok=True)
            (sub / "summary.json").write_text(json.dumps(jsonable(result), indent=2) + "\n", encoding="utf-8")
            images = plot_boards(result, result["title"], sub)
            render_html(spec, result, images, sub / "board.html")
            report = render_report(spec, result)
            (sub / "REPORT.md").write_text(report + "\n", encoding="utf-8")
            print(report, flush=True)
        print(f"=== {name} batch-size gate (MA → all-layer backprop) ===", flush=True)
        cmp_ = attach_ma(run_size_compare(name, args.n_ref, sizes, args.stream_cap))
        sub = OUT / name
        sub.mkdir(parents=True, exist_ok=True)
        plot_size_compare(cmp_["by_size"], cmp_["title"], sub)
        (sub / "batch_size.json").write_text(json.dumps(jsonable(cmp_), indent=2) + "\n", encoding="utf-8")
        print(write_batch_size_md(cmp_, sub / "BATCH_SIZE.md"), flush=True)
    write_index(OUT)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
