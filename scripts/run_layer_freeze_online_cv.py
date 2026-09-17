#!/usr/bin/env python3
"""Layer-freeze board: read PO-risk to see where to start updating.

k+1 models; model_i starts updating from the top i layers.
Incoming batch is T=1, n_ref=10000. Keep n_new large so the board
can read PO-risk directly — a small batch would need online-bootstrap.

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
from layer_freeze_cv import run_layer_freeze_cv  # noqa: E402
from streaming_po_risk import MIN_STREAM_N, REF_N  # noqa: E402

OUT = ROOT / "results" / "layer_freeze_online_cv"
DEFAULT_STREAM = {
    "electricity": MIN_STREAM_N,
    "covertype": REF_N,
}


def _standardize_from_ref(X, n_ref: int):
    mu = X[:n_ref].mean(axis=0)
    sd = X[:n_ref].std(axis=0)
    sd = np.where(sd < 1e-8, 1.0, sd)
    return (X - mu) / sd


def load_electricity():
    from sklearn.datasets import fetch_openml

    bunch = fetch_openml("electricity", version=1, as_frame=True, parser="auto")
    df = bunch.data.copy()
    y_raw = bunch.target.astype(str).to_numpy()
    Y = (y_raw == "UP").astype(int)
    drop = [c for c in df.columns if str(c).lower() == "date"]
    X = df.drop(columns=drop).to_numpy(dtype=float)
    return X, Y, "electricity (NSW, ordered in time)", list(df.drop(columns=drop).columns)


def load_covertype():
    from sklearn.datasets import fetch_covtype

    bunch = fetch_covtype()
    X = np.asarray(bunch.data, dtype=float)
    Y = (np.asarray(bunch.target) == 2).astype(int)
    return X, Y, "covertype (geographic order, class 2 vs rest)", [f"x{j}" for j in range(X.shape[1])]


LOADERS = {
    "electricity": load_electricity,
    "covertype": load_covertype,
}


def plot_boards(result: dict, title: str, out_dir: Path) -> list[str]:
    rows = result["rows"]
    k = result["k"]
    names = result["layer_names"]
    ts = [r["t"] for r in rows]
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    fig, ax = plt.subplots(figsize=(8.6, 4.0))
    ax.plot(ts, [r["po_stream"] for r in rows], color="#1f4e79", lw=2.2, marker="o", label="stream PO-risk")
    ax.axhline(result["po_base"], color="#888", ls="--", lw=1.6, label="ref-split baseline")
    for r in rows:
        if r["large_deviation"]:
            ax.scatter([r["t"]], [r["po_stream"]], s=80, color="#b33", zorder=4)
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_ylabel("PO-risk")
    ax.set_title(title + " — deviation gate")
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
新 batch 是 <b>T=1</b>。直接读 PO-risk。没有大偏差就全开；有大偏差才看从第几层开始冻。
n_new 不宜过小，否则要 online-bootstrap。
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
        "No large deviation → all trainable. Large deviation → from which layer to freeze.",
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


def main() -> int:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="electricity", choices=sorted(LOADERS) + ["both"])
    p.add_argument("--n-ref", type=int, default=REF_N)
    p.add_argument("--batch", type=int, default=0, help="incoming T=1 size; 0 = per-dataset default (≥5000)")
    p.add_argument("--max-batches", type=int, default=8)
    p.add_argument("--hidden", default="64,32")
    p.add_argument("--online-epochs", type=int, default=2)
    args = p.parse_args()
    hidden = tuple(int(x) for x in args.hidden.split(",") if x.strip())
    names = ["electricity", "covertype"] if args.dataset == "both" else [args.dataset]
    OUT.mkdir(parents=True, exist_ok=True)
    index_bits = [
        "# PO-risk board",
        "",
        "No large deviation → all trainable. Large deviation → from which layer to freeze.",
        "Tabular PO-risk: separate outcome model + propensity. T=1, n_ref=10000, n_new ≥ 5000.",
        "",
    ]
    for name in names:
        batch = args.batch if args.batch > 0 else DEFAULT_STREAM[name]
        if batch < MIN_STREAM_N:
            print(
                f"warning: n_new={batch} < {MIN_STREAM_N}; "
                "PO-risk will jitter and you would need online-bootstrap",
                flush=True,
            )
        max_batches = args.max_batches
        spec = {"name": name, "title": name}
        print(f"=== {name} n_ref={args.n_ref} n_new={batch} ===", flush=True)
        result = run_one(name, args.n_ref, batch, max_batches, hidden, args.online_epochs)
        spec["title"] = result["title"]
        sub = OUT / name
        sub.mkdir(parents=True, exist_ok=True)
        (sub / "summary.json").write_text(json.dumps(jsonable(result), indent=2) + "\n", encoding="utf-8")
        images = plot_boards(result, result["title"], sub)
        render_html(spec, result, images, sub / "board.html")
        report = render_report(spec, result)
        (sub / "REPORT.md").write_text(report + "\n", encoding="utf-8")
        print(report, flush=True)
        index_bits.append(
            f"- [{result['title']}]({name}/board.html) — large batches: {result.get('n_large', 0)}; recommend **model_{result['recommend_i']}**"
        )
    (OUT / "REPORT.md").write_text("\n".join(index_bits) + "\n", encoding="utf-8")
    html_index = """<!DOCTYPE html><meta charset="utf-8"><title>PO-risk board</title>
<h1>除非有大 deviation，否则一直 trainable</h1>
<p>直接读 PO-risk。表格：单独的 outcome + propensity。</p>
<ul>""" + "".join(
        f'<li><a href="{n}/board.html">{n}</a></li>' for n in names
    ) + "</ul>"
    (OUT / "index.html").write_text(html_index, encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
