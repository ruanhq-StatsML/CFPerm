#!/usr/bin/env python3
"""Layer-freeze online CV on a time-varying tabular stream.

Pretrain AnyMLP on D_ref (n_ref=10000, T=0). Keep k+1 copies; model_i
trains the top i layers with CosineAnnealingLR on each incoming batch
(T=1). Streaming PO-risk says which freeze-depth is OOD-stable.

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
from streaming_po_risk import REF_N  # noqa: E402

OUT = ROOT / "results" / "layer_freeze_online_cv"


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
    ts = sorted({r["t"] for r in rows})
    by_i = {i: [r for r in rows if r["i"] == i] for i in range(k + 1)}
    colors = plt.cm.tab10(np.linspace(0, 0.8, k + 1))
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    stream = [by_i[0][t]["po_stream"] for t in ts] if by_i.get(0) else []
    if stream:
        ax.plot(ts, stream, ls="--", color="#888", lw=1.6, label="stream hop (lstsq μ)")
    for i in range(k + 1):
        ys = [r["po_fit"] for r in by_i[i]]
        ax.plot(ts, ys, color=colors[i], lw=2.0, marker="o", ms=3.5, label=names[i])
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_ylabel("streaming PO-risk")
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    p = out_dir / "po_risk_by_layer.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p.name)

    mat = np.array([[by_i[i][t]["po_fit"] for t in ts] for i in range(k + 1)], dtype=float)
    fig, ax = plt.subplots(figsize=(8.6, 3.6))
    im = ax.imshow(mat, aspect="auto", cmap="magma_r", origin="lower")
    ax.set_yticks(range(k + 1), labels=names)
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_title("PO-risk heatmap — darker = more leftover hop")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    p = out_dir / "po_risk_heatmap.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p.name)

    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    stars = [by_i[0][t]["i_star"] for t in ts]
    ax.step(ts, stars, where="mid", color="#1f4e79", lw=2.2)
    ax.scatter(ts, stars, color="#1f4e79", zorder=3)
    ax.set_yticks(range(k + 1), labels=names)
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_title("i* = argmin_i PO-risk  →  freeze below this layer")
    ax.set_ylim(-0.4, k + 0.4)
    fig.tight_layout()
    p = out_dir / "i_star.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p.name)

    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    for i in range(k + 1):
        ax.plot(
            ts,
            [r["brier_new"] for r in by_i[i]],
            color=colors[i],
            lw=1.8,
            label=f"{names[i]} new",
        )
        ax.plot(
            ts,
            [r["brier_ref"] for r in by_i[i]],
            color=colors[i],
            lw=1.2,
            ls=":",
            label=f"{names[i]} ref",
        )
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_ylabel("Brier")
    ax.set_title("Side gauge (not the CV): Brier on T=1 vs reference")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    p = out_dir / "brier_tradeoff.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p.name)
    return written


def render_html(spec: dict, result: dict, images: list[str], out_path: Path) -> None:
    rec_i = result["recommend_i"]
    k = result["k"]
    freeze_below = "all frozen (model_0)" if rec_i == 0 else f"freeze below top {rec_i} layer(s); train model_{rec_i}"
    cards = []
    for img in images:
        cards.append(f'<figure><img src="{img}" alt="{img}"><figcaption>{img}</figcaption></figure>')
    rows = result["rows"]
    last_t = max(r["t"] for r in rows) if rows else 0
    last = [r for r in rows if r["t"] == last_t]
    table = [
        "<table><thead><tr><th>model</th><th>PO-risk</th><th>Brier new</th><th>Brier ref</th></tr></thead><tbody>"
    ]
    for r in last:
        mark = " ★" if r["i"] == rec_i else ""
        table.append(
            f"<tr><td>{r['name']}{mark}</td><td>{r['po_fit']:.4g}</td>"
            f"<td>{r['brier_new']:.4g}</td><td>{r['brier_ref']:.4g}</td></tr>"
        )
    table.append("</tbody></table>")
    html = f"""<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<title>Layer-freeze online CV</title>
<style>
body {{ font-family: "IBM Plex Sans", "Noto Sans SC", sans-serif; margin: 24px; color: #122; background: #f7f5f0; }}
h1 {{ font-size: 1.5rem; }}
.note {{ max-width: 820px; line-height: 1.45; }}
.rec {{ background: #1f4e79; color: #fff; padding: 12px 16px; border-radius: 8px; display: inline-block; }}
figure {{ margin: 18px 0; }}
img {{ max-width: 100%; background: #fff; border: 1px solid #ddd; }}
table {{ border-collapse: collapse; background: #fff; }}
td, th {{ border: 1px solid #ccc; padding: 6px 10px; font-variant-numeric: tabular-nums; }}
code {{ background: #eee; padding: 1px 4px; }}
</style>
</head>
<body>
<h1>哪一层该 freeze：streaming PO-risk CV</h1>
<p class="note">
{spec["title"]}. n_ref={result["n_ref"]}, hidden={result["hidden_dims"]},
k={k} layer groups → <code>model_0</code>…<code>model_{k}</code>.
新 batch 是 <b>T=1</b>。μ 来自该 freeze-depth 的 AnyMLP；PO-risk 是 leftover hop。
i* = argmin<sub>i</sub> PO-risk：只更新最上面 i* 层，下面冻住。
</p>
<p class="rec">建议：{freeze_below}（median i* = {rec_i}）</p>
{"".join(cards)}
<h2>Last batch</h2>
{"".join(table)}
<p class="note">Brier 只是 side gauge，不是 CV。HH chosen 不是这张表的 Y。</p>
</body>
</html>
"""
    out_path.write_text(html, encoding="utf-8")


def render_report(spec: dict, result: dict) -> str:
    rec_i = result["recommend_i"]
    lines = [
        f"# Layer-freeze online CV — {spec['title']}",
        "",
        "k layers → k+1 models (`model_0` … `model_k`). `model_i` trains the top i layers.",
        "Incoming batch is T=1. Reference is T=0, n_ref=10000.",
        "CV statistic = streaming PO-risk of that model's μ. i* = argmin_i PO-risk.",
        "",
        f"- hidden_dims = `{result['hidden_dims']}`",
        f"- k = {result['k']} → models 0…{result['k']}",
        f"- n_ref = {result['n_ref']}, n_batches = {result['n_batches']}",
        f"- **recommend train top {rec_i} layer(s)** (median i*). Freeze below that.",
        "",
        "| t | stream hop | " + " | ".join(result["layer_names"]) + " | i* |",
        "|---:|---:|" + "|".join(["---:"] * (result["k"] + 1)) + "|---:|",
    ]
    ts = sorted({r["t"] for r in result["rows"]})
    by = {(r["t"], r["i"]): r for r in result["rows"]}
    for t in ts:
        po_s = by[(t, 0)]["po_stream"]
        cells = " | ".join(f"{by[(t, i)]['po_fit']:.3g}" for i in range(result["k"] + 1))
        lines.append(f"| {t} | {po_s:.3g} | {cells} | {by[(t, 0)]['i_star']} |")
    lines += [
        "",
        "Read: if `model_0` (frozen) PO-risk stays high while a shallow `model_i` drops, unfreeze that far.",
        "If a deep i spikes, you over-updated and washed the reference P(Y|X) — freeze those bottom layers.",
        "",
    ]
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
    p.add_argument("--batch", type=int, default=2500)
    p.add_argument("--max-batches", type=int, default=12)
    p.add_argument("--hidden", default="64,32")
    p.add_argument("--online-epochs", type=int, default=2)
    args = p.parse_args()
    hidden = tuple(int(x) for x in args.hidden.split(",") if x.strip())
    names = ["electricity", "covertype"] if args.dataset == "both" else [args.dataset]
    OUT.mkdir(parents=True, exist_ok=True)
    index_bits = ["# Layer-freeze online CV", "", "Incoming batch is **T=1**. n_ref defaults to 10000.", ""]
    for name in names:
        batch = args.batch
        max_batches = args.max_batches
        if name == "covertype" and args.dataset == "both":
            batch = 4000
            max_batches = min(args.max_batches, 10)
        spec = {"name": name, "title": name}
        print(f"=== {name} n_ref={args.n_ref} batch={batch} ===", flush=True)
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
        index_bits.append(f"- [{result['title']}]({name}/board.html) — recommend train top **{result['recommend_i']}**")
    (OUT / "REPORT.md").write_text("\n".join(index_bits) + "\n", encoding="utf-8")
    html_index = """<!DOCTYPE html><meta charset="utf-8"><title>Layer-freeze CV</title>
<h1>Layer-freeze online CV</h1><ul>""" + "".join(
        f'<li><a href="{n}/board.html">{n}</a></li>' for n in names
    ) + "</ul>"
    (OUT / "index.html").write_text(html_index, encoding="utf-8")
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
