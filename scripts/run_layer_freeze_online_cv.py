#!/usr/bin/env python3
"""Layer-freeze board: read PO-risk to see where to start updating.

k+1 models; model_i starts updating from the top i layers.
Incoming batch is T=1, n_ref=10000. Raw PO-risk jitters on a small
n_new; a causal moving average is the stability readout. No
online-bootstrap — repeated MLP inference cannot be afforded.
MA below 2× baseline → all-layer backprop.

PO-risk nuisances are Random Forests. RF PO-risk should not collapse
suddenly; serving MSE of the MLP is likelier to break first.

MMD口径 is MMD²(X_new, X_ref) vs the T=0 reference batch — not
pairwise vs history, not last-batch layer representations.

Closed loop: OnlineRFPerm (frozen RandomForestRegressor, predict(X_new),
T=MSE−E_ref, last-two hop) marks WHEN. PO × MSE × MMD says WHAT.

    PYTHONPATH=Python/src:. python3 scripts/run_layer_freeze_online_cv.py --justify
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
from layer_freeze_cv import (  # noqa: E402
    attach_layer_dicts,
    freeze_training_of,
    layer_key,
    run_deviation_gate,
    run_layer_freeze_cv,
)
from streaming_po_risk import (  # noqa: E402
    ACTION_FREEZE,
    ACTION_KEEP,
    ACTION_TRICKY,
    ACTION_WATCH,
    ACTION_XSHIFT,
    MIN_STREAM_N,
    REF_N,
    annotate_moving_average,
    ma_window,
    moving_average,
)

from stream_dgps import ONSET_BATCH, make_gradual_concept, make_gradual_covariate  # noqa: E402

OUT = ROOT / "results" / "layer_freeze_online_cv"
DEFAULT_STREAM = {
    "electricity": MIN_STREAM_N,
    "covertype": REF_N,
    "airlines": MIN_STREAM_N,
    "bankmarketing": MIN_STREAM_N,
    "eeg": 1500,
    "dgp_concept": 2000,
    "dgp_covariate": 2000,
}
ONSET_TRUE = {
    "dgp_concept": ONSET_BATCH,
    "dgp_covariate": ONSET_BATCH,
}
INDEX_NAMES = (
    "electricity",
    "covertype",
    "airlines",
    "dgp_concept",
    "dgp_covariate",
    "bankmarketing",
    "eeg",
)
JUSTIFY = ("dgp_concept", "dgp_covariate", "bankmarketing", "eeg")
JUSTIFY_SPECS = {
    "dgp_concept": {"n_ref": REF_N, "batch": 2000, "max_batches": 10},
    "dgp_covariate": {"n_ref": REF_N, "batch": 2000, "max_batches": 10},
    "bankmarketing": {"n_ref": REF_N, "batch": MIN_STREAM_N, "max_batches": 6},
    "eeg": {"n_ref": 6000, "batch": 1500, "max_batches": 5},
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


def load_bankmarketing():
    from sklearn.datasets import fetch_openml

    bunch = fetch_openml("bank-marketing", version=1, as_frame=True, parser="auto")
    return _frame_to_xy(bunch.data.copy(), bunch.target, "bank-marketing (campaign order)")


def load_eeg():
    from sklearn.datasets import fetch_openml

    bunch = fetch_openml("eeg-eye-state", version=1, as_frame=True, parser="auto")
    return _frame_to_xy(bunch.data.copy(), bunch.target, "eeg-eye-state (time-ordered EEG)")


def load_dgp_concept():
    X, Y, meta = make_gradual_concept()
    return X, Y, meta["title"], [f"x{j}" for j in range(X.shape[1])]


def load_dgp_covariate():
    X, Y, meta = make_gradual_covariate()
    return X, Y, meta["title"], [f"x{j}" for j in range(X.shape[1])]


LOADERS = {
    "electricity": load_electricity,
    "covertype": load_covertype,
    "airlines": load_airlines,
    "bankmarketing": load_bankmarketing,
    "eeg": load_eeg,
    "dgp_concept": load_dgp_concept,
    "dgp_covariate": load_dgp_covariate,
}


def _float_arr(v) -> np.ndarray:
    a = np.asarray(v, dtype=object).ravel()
    out = np.empty(len(a), dtype=float)
    for i, x in enumerate(a):
        out[i] = np.nan if x is None else float(x)
    return out


def _layer_series(d: dict | None, k: int, n: int = 0) -> dict[str, np.ndarray]:
    d = d or {}
    out = {}
    empty = np.full(n, np.nan)
    for i in range(k + 1):
        key = layer_key(i)
        if key in d:
            out[key] = _float_arr(d[key])
        elif str(i) in d:
            out[key] = _float_arr(d[str(i)])
        else:
            out[key] = empty.copy()
    return out


def _finite_xy(ts, ys):
    ts = np.asarray(ts)
    ys = np.asarray(ys, dtype=float)
    m = np.isfinite(ys)
    return ts[m], ys[m]


def plot_layer_dicts(result: dict, title: str, out_dir: Path) -> list[str]:
    """PO-risk and MSE per freeze-depth. NaN hops were not cloned (MA quiet)."""
    k = int(result["k"])
    names = result.get("layer_names") or [f"model_{i}" for i in range(k + 1)]
    ts = [r["t"] for r in result["rows"]]
    po = _layer_series(result.get("PO_Dict"), k, n=len(ts))
    mse = _layer_series(result.get("MSE_Dict"), k, n=len(ts))
    has_po = any(np.isfinite(v).any() for v in po.values())
    has_mse = any(np.isfinite(v).any() for v in mse.values())
    if not has_po and not has_mse:
        return []
    colors = plt.cm.tab10(np.linspace(0, 0.8, k + 1))
    n_ax = int(has_po) + int(has_mse)
    fig, axes = plt.subplots(n_ax, 1, figsize=(8.6, 3.6 * n_ax), sharex=True)
    if n_ax == 1:
        axes = [axes]
    ax_i = 0
    if has_po:
        ax = axes[ax_i]
        ax_i += 1
        for i in range(k + 1):
            key = layer_key(i)
            x, y = _finite_xy(ts, po[key])
            if len(y):
                ax.plot(x, y, color=colors[i], lw=2.0, marker="o", ms=4, label=f"{key} / {names[i]}")
        ax.set_ylabel("PO-risk")
        ax.set_title(title + " — PO_Dict (from which layer the hop starts)")
        ax.legend(fontsize=8, ncol=2)
    if has_mse:
        ax = axes[ax_i]
        for i in range(k + 1):
            key = layer_key(i)
            x, y = _finite_xy(ts, mse[key])
            if len(y):
                ax.plot(x, y, color=colors[i], lw=2.0, marker="o", ms=4, label=f"{key} / {names[i]}")
        ax.set_ylabel("MSE")
        ax.set_title("MSE_Dict — same freeze-depths, prediction error on the new batch")
        ax.legend(fontsize=8, ncol=2)
    axes[-1].set_xlabel("incoming batch (T=1), large deviation only")
    fig.tight_layout()
    p = out_dir / "po_mse_by_layer.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written = [p.name]
    if has_po:
        fig, ax = plt.subplots(figsize=(8.6, 4.2))
        for i in range(k + 1):
            key = layer_key(i)
            x, y = _finite_xy(ts, po[key])
            if len(y):
                ax.plot(x, y, color=colors[i], lw=2.0, marker="o", ms=3.5, label=names[i])
        ax.set_xlabel("incoming batch (T=1), large deviation only")
        ax.set_ylabel("PO-risk (conditional on freeze-depth)")
        ax.set_title("from which layer to freeze")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        p = out_dir / "po_risk_by_layer.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        written.append(p.name)
    if has_mse:
        fig, ax = plt.subplots(figsize=(8.6, 4.2))
        for i in range(k + 1):
            key = layer_key(i)
            x, y = _finite_xy(ts, mse[key])
            if len(y):
                ax.plot(x, y, color=colors[i], lw=2.0, marker="o", ms=3.5, label=names[i])
        ax.set_xlabel("incoming batch (T=1), large deviation only")
        ax.set_ylabel("MSE on the new batch")
        ax.set_title("MSE by freeze-depth")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        p = out_dir / "mse_by_layer.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        written.append(p.name)
    return written


ACTION_COLOR = {
    ACTION_KEEP: "#2a7d4f",
    ACTION_WATCH: "#d48b16",
    ACTION_XSHIFT: "#6b4ea0",
    ACTION_TRICKY: "#666",
    ACTION_FREEZE: "#b33",
}
ACTION_LABEL = {
    ACTION_KEEP: "keep training",
    ACTION_WATCH: "watch (concept?)",
    ACTION_XSHIFT: "X shift (MMD)",
    ACTION_TRICKY: "tricky",
    ACTION_FREEZE: "freeze layers",
}


def _mmd_val(row: dict):
    """Board MMD is vs the reference batch. Alias mmd_stream is the same number."""
    v = row.get("mmd_vs_ref")
    if v is None:
        v = row.get("mmd_stream")
    return v


def _action_of(row: dict) -> str:
    return str(row.get("action") or (ACTION_FREEZE if not row.get("all_trainable", True) else ACTION_KEEP))


def _mark_onset(ax, result: dict) -> None:
    """Labeled onset (gray dash) and OnlineRFPerm hop onset (orange)."""
    true = result.get("onset_true")
    hat = result.get("onset_hat")
    if true is not None:
        ax.axvline(float(true), color="#888", ls="--", lw=1.3, zorder=2)
    if hat is not None:
        ax.axvline(float(hat), color="#c45c26", ls="-", lw=1.5, zorder=3)


def plot_po_mse_trend(result: dict, title: str, out_dir: Path) -> str | None:
    """PO, MSE, MMD vs ref, OnlineRFPerm onset. Freeze only when PO and MSE both break."""
    rows = result["rows"]
    if not rows or any(r.get("mse_stream") is None for r in rows):
        return None
    ts = np.asarray([r["t"] for r in rows])
    po = np.asarray([r["po_stream"] for r in rows], dtype=float)
    mse = np.asarray([r["mse_stream"] for r in rows], dtype=float)
    has_mmd = all(_mmd_val(r) is not None for r in rows)
    has_rf = all(r.get("rfperm_T") is not None for r in rows)
    n_new = int(rows[0]["n_new"])
    w = int(result.get("ma_window") or ma_window(n_new))
    po_ma = np.asarray([r.get("po_ma", np.nan) for r in rows], dtype=float)
    mse_ma = np.asarray([r.get("mse_ma", np.nan) for r in rows], dtype=float)
    if not np.isfinite(po_ma).all():
        po_ma = moving_average(po, w)
    if not np.isfinite(mse_ma).all():
        mse_ma = moving_average(mse, w)
    po_base = float(result["po_base"])
    mse_base = float(result.get("mse_base") or np.nanmean(mse[: max(1, len(mse) // 4)]))
    n_ax = 3 + int(has_mmd) + int(has_rf)
    heights = [2.0] * (n_ax - 1) + [0.8]
    fig, axes = plt.subplots(n_ax, 1, figsize=(8.8, 2.4 * n_ax + 0.6), sharex=True, gridspec_kw={"height_ratios": heights})
    ax_po = axes[0]
    ax_mse = axes[1]
    i_ax = 2
    ax_mmd = None
    ax_rf = None
    if has_mmd:
        ax_mmd = axes[i_ax]
        i_ax += 1
    if has_rf:
        ax_rf = axes[i_ax]
        i_ax += 1
    ax_act = axes[-1]
    ax_po.plot(ts, po, color="#9bb4cc", lw=1.0, marker="o", ms=3.5, label="raw")
    ax_po.plot(ts, po_ma, color="#1f4e79", lw=2.2, label=f"MA({w})")
    ax_po.axhline(po_base, color="#888", ls="--", lw=1.4, label="baseline")
    ax_po.axhline(2.0 * po_base, color="#b33", ls=":", lw=1.0, label="2×")
    ax_mse.plot(ts, mse, color="#c9b08b", lw=1.0, marker="o", ms=3.5, label="raw")
    ax_mse.plot(ts, mse_ma, color="#7a4e1f", lw=2.2, label=f"MA({w})")
    ax_mse.axhline(mse_base, color="#888", ls="--", lw=1.4, label="baseline")
    ax_mse.axhline(2.0 * mse_base, color="#b33", ls=":", lw=1.0, label="2×")
    if ax_mmd is not None:
        mmd = np.asarray([_mmd_val(r) for r in rows], dtype=float)
        mmd_ma = np.asarray([r.get("mmd_ma", np.nan) for r in rows], dtype=float)
        if not np.isfinite(mmd_ma).all():
            mmd_ma = moving_average(mmd, w)
        mmd_base = float(result.get("mmd_base") or 0.0)
        ax_mmd.plot(ts, mmd, color="#c5b3d9", lw=1.0, marker="o", ms=3.5, label="raw")
        ax_mmd.plot(ts, mmd_ma, color="#6b4ea0", lw=2.2, label=f"MA({w})")
        ax_mmd.axhline(mmd_base, color="#888", ls="--", lw=1.4, label="baseline")
        ax_mmd.axhline(2.0 * max(mmd_base, 1e-12), color="#b33", ls=":", lw=1.0, label="2×")
        ax_mmd.set_ylabel("MMD²(X_new, X_ref)")
        ax_mmd.legend(fontsize=8, ncol=3)
    if ax_rf is not None:
        T = np.asarray([r["rfperm_T"] for r in rows], dtype=float)
        ax_rf.plot(ts, T, color="#c45c26", lw=2.0, marker="o", ms=3.5, label="T = MSE − E_ref")
        ax_rf.axhline(0.0, color="#888", ls="--", lw=1.0)
        hops = [r for r in rows if r.get("rfperm_hop")]
        if hops:
            ax_rf.scatter([r["t"] for r in hops], [r["rfperm_T"] for r in hops], s=70, color="#c45c26", zorder=4, label="hop")
        ax_rf.set_ylabel("OnlineRFPerm T")
        ax_rf.legend(fontsize=8, ncol=2)
    y_map = {ACTION_KEEP: 0, ACTION_WATCH: 1, ACTION_XSHIFT: 2, ACTION_TRICKY: 3, ACTION_FREEZE: 4}
    act_y = [y_map.get(_action_of(r), 0) for r in rows]
    ax_act.step(ts, act_y, where="mid", color="#333", lw=2.0)
    seen = set()
    for r, x, y in zip(rows, ts, act_y):
        act = _action_of(r)
        color = ACTION_COLOR.get(act, "#333")
        ax_po.scatter([x], [r["po_stream"]], s=55, color=color, zorder=4, label=ACTION_LABEL.get(act, act) if act not in seen else None)
        ax_mse.scatter([x], [r["mse_stream"]], s=55, color=color, zorder=4)
        if ax_mmd is not None:
            ax_mmd.scatter([x], [_mmd_val(r)], s=55, color=color, zorder=4)
        ax_act.scatter([x], [y], s=55, color=color, zorder=4)
        seen.add(act)
    for ax in (ax_po, ax_mse, ax_mmd, ax_rf, ax_act):
        if ax is not None:
            _mark_onset(ax, result)
    ax_po.set_ylabel("PO-risk")
    ax_po.set_title(title + " — PO vs MSE vs MMD²(X_new, X_ref); OnlineRFPerm onset")
    ax_po.legend(fontsize=8, ncol=3)
    ax_mse.set_ylabel("serving MSE")
    ax_mse.legend(fontsize=8, ncol=3)
    ax_act.set_yticks([0, 1, 2, 3, 4], labels=["keep", "watch", "X shift", "tricky", "freeze"])
    ax_act.set_ylim(-0.4, 4.4)
    ax_act.set_xlabel("incoming batch (T=1)")
    hat = result.get("onset_hat")
    true = result.get("onset_true")
    ax_act.set_title(
        f"onset_hat={hat} (OnlineRFPerm hop); labeled={true}. MSE likelier to break first than RF PO-risk."
    )
    fig.tight_layout()
    p = out_dir / "po_mse_trend.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p.name


def plot_boards(result: dict, title: str, out_dir: Path) -> list[str]:
    result = attach_layer_dicts(result)
    rows = result["rows"]
    k = result["k"]
    names = result["layer_names"]
    ts = [r["t"] for r in rows]
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []

    trend = plot_po_mse_trend(result, title, out_dir)
    if trend:
        written.append(trend)

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
        ax.scatter([r["t"]], [r["po_stream"]], s=70, color=ACTION_COLOR[_action_of(r)], zorder=4)
    ax.set_xlabel("incoming batch (T=1)")
    ax.set_ylabel("PO-risk")
    tag = ACTION_LABEL.get(result.get("board_action"), "read PO × MSE")
    ax.set_title(title + f" — PO-risk ({tag}; no bootstrap)")
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
    ax.set_title("which layers stop training (only when PO and MSE both break)")
    ax.set_ylim(-0.4, k + 0.4)
    fig.tight_layout()
    p = out_dir / "freeze_from.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p.name)

    written.extend(plot_layer_dicts(result, title, out_dir))
    return written


def render_html(spec: dict, result: dict, images: list[str], out_path: Path) -> None:
    result = attach_layer_dicts(result)
    rec_i = result["recommend_i"]
    k = result["k"]
    board = result.get("board_action") or ACTION_KEEP
    rec_zh = result.get("freeze_training_zh") or freeze_training_of(rec_i, k)["freeze_training_zh"]
    if board == ACTION_KEEP:
        rec = "不冻任何一层的 training — 接着 train"
    elif board == ACTION_WATCH:
        rec = "再观察 — 先不冻任何一层的 training"
    elif board == ACTION_XSHIFT:
        rec = "MSE 崩了但 PO 正常 — 不是 concept drift，看 MMD：X 在动"
    elif board == ACTION_TRICKY:
        rec = "MSE 崩了但 PO 和 MMD 都正常 — tricky，不是 P(Y|X) 也不是 X"
    else:
        rec = f"冻住 — {rec_zh}"
    n_new = result["rows"][0]["n_new"] if result["rows"] else ""
    cards = []
    for img in images:
        cards.append(f'<figure><img src="{img}" alt="{img}"><figcaption>{img}</figcaption></figure>')
    table = [
        "<table><thead><tr><th>t</th><th>PO-risk</th><th>PO MA</th><th>MSE</th><th>MSE MA</th><th>MMD vs ref</th><th>MMD MA</th><th>RFPerm T</th><th>hop</th><th>action</th><th>冻结哪一层的 training</th></tr></thead><tbody>"
    ]
    for r in result["rows"]:
        act = _action_of(r)
        po_ma = r.get("po_ma")
        mse = r.get("mse_stream")
        mse_ma = r.get("mse_ma")
        mmd = _mmd_val(r)
        mmd_ma = r.get("mmd_ma")
        po_ma_s = "" if po_ma is None else f"{float(po_ma):.4g}"
        mse_s = "" if mse is None else f"{float(mse):.4g}"
        mse_ma_s = "" if mse_ma is None else f"{float(mse_ma):.4g}"
        mmd_s = "" if mmd is None else f"{float(mmd):.4g}"
        mmd_ma_s = "" if mmd_ma is None else f"{float(mmd_ma):.4g}"
        T = r.get("rfperm_T")
        T_s = "" if T is None else f"{float(T):.4g}"
        hop_s = "●" if r.get("rfperm_hop") else ""
        freeze_s = r.get("freeze_training_zh") or ACTION_LABEL.get(act, act)
        table.append(
            f"<tr><td>{r['t']}</td><td>{r['po_stream']:.4g}</td>"
            f"<td>{po_ma_s}</td><td>{mse_s}</td><td>{mse_ma_s}</td>"
            f"<td>{mmd_s}</td><td>{mmd_ma_s}</td><td>{T_s}</td><td>{hop_s}</td>"
            f"<td>{ACTION_LABEL.get(act, act)}</td><td>{freeze_s}</td></tr>"
        )
    table.append("</tbody></table>")
    layer_tab = ""
    last_freeze = next((r for r in reversed(result["rows"]) if r.get("action") == ACTION_FREEZE and r.get("layers")), None)
    if last_freeze:
        layer_tab = (
            "<h2>Last freeze batch — PO-risk 和 MSE 同时读</h2>"
            "<table><thead><tr><th>layer</th><th>model</th><th>PO-risk</th><th>MSE</th></tr></thead><tbody>"
        )
        for x in last_freeze["layers"]:
            mark = " ★" if x["i"] == last_freeze["i_star"] else ""
            mse_s = f"{x['mse']:.4g}" if x.get("mse") is not None else ""
            layer_tab += (
                f"<tr><td>{x.get('layer', layer_key(x['i']))}</td>"
                f"<td>{x['name']}{mark}</td><td>{x['po_fit']:.4g}</td><td>{mse_s}</td></tr>"
            )
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
<h1>该冻结哪一层的 training</h1>
<p class="note">
{spec["title"]}. n_ref={result["n_ref"]}, n_new={n_new}, hidden={result["hidden_dims"]}.
看板的输出就是：<b>该冻结哪一层的 training</b>（仅当 PO 和 MSE 都崩）。
PO-risk 用 RF，不该突然崩；<b>更可能先崩的是 serving MSE</b>。
MSE 崩了但 PO 正常 → 不是 concept drift。MMD 口径是 <b>MMD²(X_new, X_ref)</b>，
跟 reference batch 比，不是跟历史所有 batch 的 pairwise 均值，也不是上个 batch 的层 representation。
闭环：<b>OnlineRFPerm</b> 冻参考窗 RF，T=MSE−E_ref，last-two hop 标 <b>shift-onset</b>（onset_hat={result.get("onset_hat")}，labeled={result.get("onset_true")}）。
看板说这是哪种 shift。PO 崩了模型没崩 → 再观察，先不冻。
两都崩 → 开 freeze-depth，<code>model_i</code> 只训 top i。
不做 online-bootstrap。
</p>
<p class="rec">{rec}</p>
{"".join(cards)}
<h2>每一段 — PO × MSE</h2>
{"".join(table)}
{layer_tab}
</body>
</html>
"""
    out_path.write_text(html, encoding="utf-8")


def _fmt_layer_row(layers, field, k):
    by = {x["i"]: x.get(field) for x in layers}
    cells = []
    for i in range(k + 1):
        v = by.get(i)
        cells.append("" if v is None else f"{float(v):.3g}")
    return " | ".join(cells)


def render_report(spec: dict, result: dict) -> str:
    result = attach_layer_dicts(result)
    rec_i = result["recommend_i"]
    n_new = result["rows"][0]["n_new"] if result["rows"] else ""
    k = result["k"]
    board = result.get("board_action") or ACTION_KEEP
    rec = ACTION_LABEL.get(board, board)
    if board == ACTION_FREEZE and rec_i < k:
        rec = f"freeze from model_{rec_i}"
    mse_base = result.get("mse_base")
    mse_base_s = "" if mse_base is None else f"{float(mse_base):.4g}"
    lines = [
        f"# PO × MSE board — {spec['title']}",
        "",
        "RF PO-risk should not collapse first; serving MSE is the series that breaks.",
        "MMD口径 is MMD²(X_new, X_ref) — vs the reference batch, not history, not layer reps.",
        "PO broken, MSE holds → watch. MSE broken, PO quiet → read MMD vs ref.",
        "OnlineRFPerm marks the shift-onset (frozen RF, T=MSE−E_ref, last-two hop).",
        f"T=1 on the incoming batch. n_ref={result['n_ref']}, n_new={n_new}.",
        f"po_base={result['po_base']:.4g}. mse_base={mse_base_s}. mmd_base={result.get('mmd_base')}.",
        f"onset_hat={result.get('onset_hat')}, onset_rank={result.get('onset_rank')}, labeled onset={result.get('onset_true')}.",
        "",
        f"- hidden_dims = `{result['hidden_dims']}`, k = {k}",
        f"- n_batches = {result['n_batches']}, n_watch = {result.get('n_watch', 0)}, n_x_shift = {result.get('n_x_shift', 0)}, n_tricky = {result.get('n_tricky', 0)}, n_freeze = {result.get('n_freeze', 0)}, n_rfperm_hop = {result.get('n_rfperm_hop', 0)}",
        f"- **{rec}**",
        "",
        "| t | PO-risk | PO MA | MSE | MSE MA | MMD | MMD MA | RFPerm T | hop | action | freeze training |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for r in result["rows"]:
        act = _action_of(r)
        po_ma = r.get("po_ma")
        mse = r.get("mse_stream")
        mse_ma = r.get("mse_ma")
        po_ma_s = "" if po_ma is None else f"{float(po_ma):.3g}"
        mse_s = "" if mse is None else f"{float(mse):.3g}"
        mse_ma_s = "" if mse_ma is None else f"{float(mse_ma):.3g}"
        mmd = _mmd_val(r)
        mmd_ma = r.get("mmd_ma")
        mmd_s = "" if mmd is None else f"{float(mmd):.3g}"
        mmd_ma_s = "" if mmd_ma is None else f"{float(mmd_ma):.3g}"
        freeze_s = r.get("freeze_training") or ""
        T = r.get("rfperm_T")
        T_s = "" if T is None else f"{float(T):.3g}"
        hop_s = "yes" if r.get("rfperm_hop") else ""
        lines.append(
            f"| {r['t']} | {r['po_stream']:.3g} | {po_ma_s} | {mse_s} | {mse_ma_s} | {mmd_s} | {mmd_ma_s} | {T_s} | {hop_s} | {ACTION_LABEL.get(act, act)} | {freeze_s} |"
        )
    freeze_rows = [r for r in result["rows"] if r.get("action") == ACTION_FREEZE and r.get("layers")]
    keys = result.get("layer_keys") or [layer_key(i) for i in range(k + 1)]
    if freeze_rows:
        lines += [
            "",
            "PO_Dict (conditional on freeze-depth):",
            "",
            "| t | " + " | ".join(keys) + " | freeze from |",
            "|---:|" + "|".join(["---:"] * (k + 1)) + "|---|",
        ]
        for r in freeze_rows:
            action = "all trainable" if r["all_trainable"] else f"freeze from {r['freeze_from']}"
            lines.append(f"| {r['t']} | {_fmt_layer_row(r['layers'], 'po_fit', k)} | {action} |")
        has_mse = any(x.get("mse") is not None for r in freeze_rows for x in r["layers"])
        if has_mse:
            lines += [
                "",
                "MSE_Dict (new-batch prediction error, same freeze-depths):",
                "",
                "| t | " + " | ".join(keys) + " | MSE-best |",
                "|---:|" + "|".join(["---:"] * (k + 1)) + "|---|",
            ]
            for r in freeze_rows:
                star = r.get("i_star_mse")
                star_s = "" if star is None else f"layer{star}"
                lines.append(f"| {r['t']} | {_fmt_layer_row(r['layers'], 'mse', k)} | {star_s} |")
    lines += ["", "Read the two trends. Freeze only when both break.", ""]
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
    if name in ONSET_TRUE:
        result["onset_true"] = int(ONSET_TRUE[name])
        result["dgp_kind"] = "concept" if name == "dgp_concept" else "covariate"
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
        "OnlineRFPerm (frozen RF.predict(X_new), T=MSE−E_ref) marks shift-onset.",
        "PO × MSE × MMD²(X_new, X_ref) says what kind of shift it was.",
        "PO+MSE both break → freeze that layer's training.",
        "RF PO-risk should not collapse first; serving MSE is likelier to break.",
        "MSE broken, PO quiet → not concept drift; read MMD vs the reference batch.",
        "No online-bootstrap.",
        "",
    ]
    html_items = []
    for name in INDEX_NAMES:
        sub = out / name
        if (sub / "board.html").exists():
            bits.append(f"- [{name} freeze board]({name}/board.html)")
            html_items.append(f'<li><a href="{name}/board.html">{name}</a>')
        else:
            html_items.append(f"<li>{name}")
        if (sub / "po_mse_trend.png").exists():
            bits.append(f"- [{name} PO×MSE trend]({name}/po_mse_trend.png)")
            html_items.append(f' · <a href="{name}/po_mse_trend.png">trend</a>')
        if (sub / "batch_size_gate.png").exists():
            bits.append(f"- [{name} MA gate]({name}/batch_size_gate.png)")
            html_items.append(f' · <a href="{name}/batch_size_gate.png">MA</a></li>')
        else:
            html_items.append("</li>")
    if (out / "dgp_contrast.png").exists():
        bits.append("- [DGP concept vs covariate](dgp_contrast.png)")
        html_items.append('<li><a href="dgp_contrast.png">DGP contrast</a></li>')
    if (out / "JUSTIFY.md").exists():
        bits.append("- [justify table](JUSTIFY.md)")
        html_items.append('<li><a href="JUSTIFY.md">justify</a></li>')
    (out / "REPORT.md").write_text("\n".join(bits) + "\n", encoding="utf-8")
    html = """<!DOCTYPE html><meta charset="utf-8"><title>PO × MSE board</title>
<h1>PO × MSE 对照</h1>
<p>闭环：OnlineRFPerm 标 onset（冻住的 RandomForestRegressor.predict(X_new)）。MMD 口径是 MMD²(X_new, X_ref)。RF PO-risk 不该先崩；MSE 先崩再看 MMD。PO 崩模型没崩 → 再观察。</p>
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
            result = attach_layer_dicts(json.loads(board.read_text(encoding="utf-8")))
            spec = {"name": name, "title": result.get("title", name)}
            images = plot_boards(result, spec["title"], sub)
            render_html(spec, result, images, sub / "board.html")
            (sub / "REPORT.md").write_text(render_report(spec, result) + "\n", encoding="utf-8")
            board.write_text(json.dumps(jsonable(result), indent=2) + "\n", encoding="utf-8")


def write_freeze_board(name: str, result: dict) -> None:
    spec = {"name": name, "title": result.get("title", name)}
    sub = OUT / name
    sub.mkdir(parents=True, exist_ok=True)
    (sub / "summary.json").write_text(json.dumps(jsonable(result), indent=2) + "\n", encoding="utf-8")
    images = plot_boards(result, result["title"], sub)
    render_html(spec, result, images, sub / "board.html")
    report = render_report(spec, result)
    (sub / "REPORT.md").write_text(report + "\n", encoding="utf-8")
    print(report, flush=True)


def plot_dgp_contrast(out: Path) -> str | None:
    """Side-by-side concept vs covariate: WHEN (OnlineRFPerm) vs WHAT (PO/MSE/MMD)."""
    paths = [out / "dgp_concept" / "summary.json", out / "dgp_covariate" / "summary.json"]
    if not all(p.exists() for p in paths):
        return None
    fig, axes = plt.subplots(4, 2, figsize=(11.2, 9.2), sharex=True)
    titles = []
    for col, path in enumerate(paths):
        result = json.loads(path.read_text(encoding="utf-8"))
        rows = result["rows"]
        ts = [r["t"] for r in rows]
        titles.append(result.get("title", path.parent.name))
        series = [
            (np.asarray([r["po_stream"] for r in rows], dtype=float), "PO-risk", result.get("po_base")),
            (np.asarray([r["mse_stream"] for r in rows], dtype=float), "serving MSE", result.get("mse_base")),
            (np.asarray([_mmd_val(r) for r in rows], dtype=float), "MMD²(X_new, X_ref)", result.get("mmd_base")),
            (np.asarray([0.0 if r.get("rfperm_T") is None else float(r["rfperm_T"]) for r in rows], dtype=float), "OnlineRFPerm T", 0.0),
        ]
        for row, (ys, ylab, base) in enumerate(series):
            ax = axes[row, col]
            ax.plot(ts, ys, color="#1f4e79", lw=2.0, marker="o", ms=3.5)
            if base is not None:
                ax.axhline(float(base), color="#888", ls="--", lw=1.0)
                if row < 3:
                    ax.axhline(2.0 * max(float(base), 1e-12), color="#b33", ls=":", lw=1.0)
            _mark_onset(ax, result)
            if col == 0:
                ax.set_ylabel(ylab)
            hops = [r["t"] for r in rows if r.get("rfperm_hop")]
            if hops and row == 3:
                ax.scatter(hops, [ys[t] for t in hops], s=70, color="#c45c26", zorder=4)
        axes[0, col].set_title(titles[col], fontsize=10)
        axes[-1, col].set_xlabel("incoming batch (T=1)")
    fig.suptitle("labeled onset = gray dash; OnlineRFPerm hop = orange. Concept: MMD quiet. Covariate: MMD fires.", y=1.01)
    fig.tight_layout()
    dest = out / "dgp_contrast.png"
    fig.savefig(dest, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return dest.name


def write_justify_md(results: list[dict], out_path: Path) -> str:
    lines = [
        "# Justify: OnlineRFPerm onset × PO / MSE / MMD vs ref",
        "",
        "WHEN = frozen `RandomForestRegressor().predict(X_new)`, T = MSE − E_ref, last-two hop.",
        "WHAT = PO-risk (RF nuisances) × serving MSE × MMD²(X_new, X_ref).",
        "Concept DGP: P(X) fixed, β rotates after labeled onset → MMD stays quiet, MSE/PO move.",
        "Covariate DGP: P(Y|X) fixed, μ(X) walks → MMD fires, PO stays quieter; MSE-only is X shift.",
        "",
        "| dataset | n_new | onset_true | onset_hat | onset_rank | board | n_watch | n_x_shift | n_tricky | n_freeze | n_hop |",
        "|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        lines.append(
            f"| {r.get('dataset')} | {r.get('n_new')} | {r.get('onset_true')} | {r.get('onset_hat')} | "
            f"{r.get('onset_rank')} | {r.get('board_action')} | {r.get('n_watch', 0)} | "
            f"{r.get('n_x_shift', 0)} | {r.get('n_tricky', 0)} | {r.get('n_freeze', 0)} | {r.get('n_rfperm_hop', 0)} |"
        )
    lines += ["", "Freeze only when PO and MSE both break. No online-bootstrap.", ""]
    text = "\n".join(lines)
    out_path.write_text(text + "\n", encoding="utf-8")
    return text


def run_justify(hidden, online_epochs: int) -> list[dict]:
    """DGP contrast + extra datasets. Does not touch the size-grid boards."""
    results = []
    for name in JUSTIFY:
        spec = JUSTIFY_SPECS[name]
        print(
            f"=== justify {name} n_ref={spec['n_ref']} n_new={spec['batch']} max_batches={spec['max_batches']} ===",
            flush=True,
        )
        result = run_one(name, spec["n_ref"], spec["batch"], spec["max_batches"], hidden, online_epochs)
        write_freeze_board(name, result)
        results.append(result)
        print(
            f"done {name}: board={result.get('board_action')} onset_hat={result.get('onset_hat')} "
            f"labeled={result.get('onset_true')} n_x_shift={result.get('n_x_shift')} n_watch={result.get('n_watch')}",
            flush=True,
        )
    plot_dgp_contrast(OUT)
    print(write_justify_md(results, OUT / "JUSTIFY.md"), flush=True)
    write_index(OUT)
    return results


def main() -> int:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--dataset", default="electricity", choices=sorted(LOADERS) + ["both", "all", "justify"])
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
    p.add_argument(
        "--justify",
        action="store_true",
        help="gradual concept vs covariate DGPs + extra datasets; OnlineRFPerm onset. No size grid.",
    )
    args = p.parse_args()
    hidden = tuple(int(x) for x in args.hidden.split(",") if x.strip())
    if args.justify or args.dataset == "justify":
        OUT.mkdir(parents=True, exist_ok=True)
        run_justify(hidden, args.online_epochs)
        print("wrote", OUT)
        return 0
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
            if name in JUSTIFY_SPECS:
                js = JUSTIFY_SPECS[name]
                n_ref, batch, max_batches = js["n_ref"], js["batch"], js["max_batches"]
            else:
                n_ref = args.n_ref
                batch = args.batch if args.batch > 0 else DEFAULT_STREAM[name]
                max_batches = args.max_batches
            if batch < MIN_STREAM_N:
                print(
                    f"note: n_new={batch} < {MIN_STREAM_N}; "
                    "raw PO-risk jitters — read the moving average, do not bootstrap",
                    flush=True,
                )
            print(f"=== {name} n_ref={n_ref} n_new={batch} ===", flush=True)
            result = run_one(name, n_ref, batch, max_batches, hidden, args.online_epochs)
            write_freeze_board(name, result)
        if name in JUSTIFY:
            continue
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
