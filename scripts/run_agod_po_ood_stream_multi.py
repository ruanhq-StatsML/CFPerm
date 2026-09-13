#!/usr/bin/env python3
"""Streaming real-data: PO-risk as OOD score (√PO) vs logistic DRE.

batch_size=100, preserve gradual order, no synth.

  PYTHONPATH=. python3 scripts/run_agod_po_ood_stream_multi.py \\
    --batch-size 100 --n-batches 20
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import accuracy_score, mean_squared_error

from agod.po_iptw import dre_weights, instance_po_risk, po_iptw_weights

warnings.filterwarnings("ignore", category=UserWarning)

MODES = ("uniform", "prop", "sqrt", "inv", "dre")
# real only — continuous=mse, discrete=acc
DATASETS = {
    "affec": "mse",
    "tencent": "mse",
    "msrvtt": "acc",
    "coco_time": "acc",
    "fashion_iq": "acc",
    "indiana_cxr": "acc",
}


def batch_po_risk(X0, y0, X1, y1, *, seed: int, task: str) -> float:
    if len(X0) < 16 or len(X1) < 16:
        return 0.0
    if task == "acc":
        clf = RandomForestClassifier(
            n_estimators=20, max_depth=4, min_samples_leaf=2, random_state=seed, n_jobs=1
        )
        clf.fit(X0, y0.astype(int))
        e0 = 1.0 - float(accuracy_score(y0.astype(int), clf.predict(X0)))
        e1 = 1.0 - float(accuracy_score(y1.astype(int), clf.predict(X1)))
        return max(e1 - e0, 0.0)
    rf = RandomForestRegressor(
        n_estimators=20, max_depth=4, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    rf.fit(X0, y0)
    e0 = float(np.mean(np.abs(y0 - rf.predict(X0))))
    e1 = float(np.mean(np.abs(y1 - rf.predict(X1))))
    return max(e1 - e0, 0.0)


def _pca(X: np.ndarray, d: int, seed: int) -> np.ndarray:
    if X.shape[1] <= d:
        return X.astype(np.float32)
    return PCA(n_components=d, random_state=seed).fit_transform(X).astype(np.float32)


def _coarse_labels(y: np.ndarray, top_k: int = 4) -> np.ndarray:
    """Map rare classes → other so Acc is measurable in bs=100 streams."""
    y = np.asarray(y).ravel().astype(int)
    u, c = np.unique(y, return_counts=True)
    keep = set(u[np.argsort(-c)[:top_k]].tolist())
    other = top_k
    out = np.empty_like(y)
    remap = {lab: i for i, lab in enumerate(sorted(keep))}
    for i, lab in enumerate(y):
        out[i] = remap[lab] if lab in keep else other
    return out


def load_affec(root: Path, max_n: int, seed: int, pca_d: int):
    cache = root / "results/affec_fsds/affec_fsds_xyw_cache.npz"
    if not cache.is_file():
        return None
    z = np.load(cache, allow_pickle=True)
    X, y = z["X"].astype(np.float32), z["Y"].astype(np.float64)
    n = min(max_n, len(X))
    return _pca(X[:n], pca_d, seed), y[:n]


def load_msrvtt(root: Path, max_n: int, seed: int, pca_d: int):
    d = root / "data/msrvtt/packed"
    if not (d / "video_feat.npy").is_file():
        return None
    X = np.concatenate(
        [
            np.load(d / "video_feat.npy"),
            np.load(d / "audio_feat.npy"),
            np.load(d / "text_feat.npy"),
        ],
        axis=1,
    )
    y = np.load(d / "labelsmsr.npy").astype(np.int64)
    Xp = _pca(X, pca_d, seed)
    # raw pack is class-sorted; PC1 order → gradual feature drift + mixed labels
    order = np.argsort(Xp[:, 0])
    Xp, y = Xp[order], y[order]
    n = min(max_n, len(Xp))
    return Xp[:n], y[:n]


def load_img_txt(root: Path, name: str, max_n: int, seed: int, pca_d: int):
    d = root / "data/img_txt" / name
    if not (d / "img_feats.npy").is_file():
        return None
    X = np.concatenate(
        [np.load(d / "img_feats.npy"), np.load(d / "txt_feats.npy")], axis=1
    )
    y = _coarse_labels(np.load(d / "labels.npy"), top_k=4)
    n = min(max_n, len(X))
    # keep file order (coco_time_order = temporal gradual shift)
    return _pca(X[:n], pca_d, seed), y[:n]


def load_tencent(root: Path, max_n: int, seed: int, pca_d: int):
    p = root / "results/tencent_gr/user_feats_lean.parquet"
    if not p.is_file():
        return None
    import pandas as pd

    df = pd.read_parquet(p)
    y_col = "life_ctcvr" if "life_ctcvr" in df.columns else "arpu_mean"
    drop = {"user_id", y_col}
    num = [c for c in df.columns if c not in drop and np.issubdtype(df[c].dtype, np.number)]
    X = df[num].fillna(0).to_numpy(np.float32)
    y = df[y_col].fillna(0).to_numpy(np.float64)
    n = min(max_n, len(X))
    return _pca(X[:n], pca_d, seed), y[:n]


def load_dataset(name: str, root: Path, max_n: int, seed: int, pca_d: int):
    if name == "affec":
        return load_affec(root, max_n, seed, pca_d)
    if name == "msrvtt":
        return load_msrvtt(root, max_n, seed, pca_d)
    if name == "coco_time":
        return load_img_txt(root, "coco_time_order", max_n, seed, pca_d)
    if name == "fashion_iq":
        return load_img_txt(root, "fashion_iq", max_n, seed, pca_d)
    if name == "indiana_cxr":
        return load_img_txt(root, "indiana_cxr", max_n, seed, pca_d)
    if name == "tencent":
        return load_tencent(root, max_n, seed, pca_d)
    raise ValueError(name)


def make_stream(X, y, batch: int):
    return [(X[i : i + batch], y[i : i + batch]) for i in range(0, len(X) - batch + 1, batch)]


def fit_model(X, y, w, seed, task: str):
    if task == "acc":
        m = RandomForestClassifier(
            n_estimators=40, max_depth=8, min_samples_leaf=2, random_state=seed, n_jobs=1
        )
        m.fit(X, y.astype(int), sample_weight=w)
        return m
    m = RandomForestRegressor(
        n_estimators=30, max_depth=6, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    m.fit(X, y, sample_weight=w)
    return m


def eval_next(model, X, y, task: str) -> float:
    if task == "acc":
        return float(accuracy_score(y.astype(int), model.predict(X)))
    return float(mean_squared_error(y, model.predict(X)))


def instance_score(model, X, y, *, batch_po: float, task: str) -> np.ndarray:
    """PO-risk OOD score from prev probe on current batch."""
    if task == "acc":
        y = np.asarray(y).ravel().astype(int)
        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(X)
            classes = list(model.classes_)
            p_true = np.zeros(len(y), float)
            for i, yi in enumerate(y):
                if yi in classes:
                    p_true[i] = float(proba[i, classes.index(yi)])
            po = 1.0 - p_true  # high PO = OOD / hard under probe
        else:
            wrong = (y != model.predict(X).astype(int)).astype(float)
            po = 0.25 + wrong
        return instance_po_risk(po, np.zeros_like(po), batch_po=batch_po, mix=0.5)
    pred = model.predict(X)
    return instance_po_risk(y, pred, batch_po=batch_po, mix=0.5)


def _pack(xs: List[float]) -> dict:
    a = np.asarray(xs, float)
    return {
        "mean": float(a.mean()) if len(a) else float("nan"),
        "std": float(a.std()) if len(a) else float("nan"),
        "p50": float(np.median(a)) if len(a) else float("nan"),
        "p90": float(np.percentile(a, 90)) if len(a) else float("nan"),
        "traj": [float(v) for v in a],
    }


def run_mode(stream, mode: str, seed: int, task: str) -> dict:
    cur_m: List[float] = []
    next_m: List[float] = []
    batch_po: List[float] = []
    X0, y0 = stream[0]
    probe = fit_model(X0, y0, np.ones(len(y0)), seed, task)

    for t in range(1, len(stream)):
        Xp, yp = stream[t - 1]
        Xc, yc = stream[t]
        po_b = batch_po_risk(Xp, yp, Xc, yc, seed=seed + t, task=task)
        batch_po.append(po_b)

        if mode == "dre":
            w = dre_weights(Xp, Xc, seed=seed + t)
        else:
            po_row = instance_score(probe, Xc, yc, batch_po=po_b, task=task)
            w = po_iptw_weights(po_row, mode=mode)  # type: ignore[arg-type]

        model = fit_model(Xc, yc, w, seed + 17 * t, task)
        cur_m.append(eval_next(model, Xc, yc, task))
        if t + 1 < len(stream):
            Xn, yn = stream[t + 1]
            next_m.append(eval_next(model, Xn, yn, task))
        probe = fit_model(Xc, yc, w, seed + 31 * t, task)

    key = "acc" if task == "acc" else "mse"
    return {
        "mode": mode,
        "task": task,
        "metric": key,
        "batch_po": batch_po,
        "cur": _pack(cur_m),
        "next": _pack(next_m),
    }


def table_md(all_res: Dict[str, Dict[str, dict]], tasks: Dict[str, str]) -> str:
    lines = [
        "# PO-risk as OOD score — streaming real-data (bs=100)",
        "",
        "No synth. Continuous → next **MSE** (↓); discrete → next **Acc** (↑).",
        "PO-√ soft IPTW vs logistic density-ratio (DRE).",
        "",
        "| dataset | task | uniform | prop | sqrt | inv | dre | best | √PO vs DRE |",
        "|---|---|---:|---:|---:|---:|---:|---|---|",
    ]
    wins = {m: 0 for m in MODES}
    sqrt_beats_dre = 0
    n_cmp = 0
    for ds, results in all_res.items():
        task = tasks[ds]
        means = {m: results[m]["next"]["mean"] for m in MODES}
        if task == "acc":
            best = max(means, key=means.get)
            beat = means["sqrt"] > means["dre"]
            gap = means["sqrt"] - means["dre"]
            gap_s = f"+{gap:.4f} Acc" if beat else f"{gap:.4f} Acc"
        else:
            best = min(means, key=means.get)
            beat = means["sqrt"] < means["dre"]
            gap = means["dre"] - means["sqrt"]
            gap_s = f"−{gap:.4f} MSE" if beat else f"+{-gap:.4f} MSE"
        if beat:
            sqrt_beats_dre += 1
        n_cmp += 1
        wins[best] += 1
        cells = [f"{means[m]:.4f}" for m in MODES]
        for i, m in enumerate(MODES):
            if m == best:
                cells[i] = f"**{cells[i]}**"
            elif m == "sqrt":
                cells[i] = f"*{cells[i]}*"
        lines.append(
            f"| `{ds}` | {task} | " + " | ".join(cells) + f" | `{best}` | {gap_s} |"
        )
    lines += [
        "",
        f"**Wins:** " + ", ".join(f"`{m}`={wins[m]}" for m in MODES),
        f"**sqrt vs dre (head-to-head):** `{sqrt_beats_dre}/{n_cmp}` favor PO-√ OOD over DRE.",
        "",
        "```python",
        "w_i = np.sqrt(PO-risk(X_i, Y_i, T_i=1))  # OOD score → soft IPTW",
        "model.fit(X, y, sample_weight=w / w.mean())",
        "```",
        "",
        "DRE baseline: logistic `w ∝ p(cur|x)/p(ref|x)` on X only — ignores label risk.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--datasets", nargs="+", default=list(DATASETS.keys()))
    ap.add_argument("--n-batches", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--max-n", type=int, default=2500)
    ap.add_argument("--pca-d", type=int, default=64)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/agod_po_ood_stream"))
    args = ap.parse_args()

    need = args.n_batches * args.batch_size
    all_res: Dict[str, Dict[str, dict]] = {}
    tasks: Dict[str, str] = {}
    skipped: List[str] = []

    for name in args.datasets:
        if name == "synth":
            print("[skip] synth: real-data only")
            skipped.append(name)
            continue
        task = DATASETS.get(name)
        if task is None:
            print(f"[skip] {name}: unknown")
            skipped.append(name)
            continue
        packed = load_dataset(name, args.root, max(need, args.max_n), args.seed, args.pca_d)
        if packed is None:
            print(f"[skip] {name}: missing data")
            skipped.append(name)
            continue
        X, y = packed
        if len(X) < need:
            print(f"[skip] {name}: need>={need}, got {len(X)}")
            skipped.append(name)
            continue
        if task == "acc" and len(np.unique(y[:need])) < 2:
            print(f"[skip] {name}: <2 classes")
            skipped.append(name)
            continue
        stream = make_stream(X[:need], y[:need], args.batch_size)
        print(
            f"=== {name} task={task} batches={len(stream)} bs={args.batch_size} "
            f"d={X.shape[1]} nuniq={len(np.unique(y[:need]))} ===",
            flush=True,
        )
        results = {}
        for mode in MODES:
            print(f"  [{mode}] ...", flush=True)
            results[mode] = run_mode(stream, mode, args.seed, task)
            n = results[mode]["next"]
            print(f"    next {results[mode]['metric']} mean={n['mean']:.4f} std={n['std']:.4f}")
        all_res[name] = results
        tasks[name] = task

    args.out.mkdir(parents=True, exist_ok=True)
    payload = {
        "batch_size": args.batch_size,
        "n_batches": args.n_batches,
        "modes": list(MODES),
        "skipped": skipped,
        "tasks": tasks,
        "results": all_res,
    }
    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = table_md(all_res, tasks)
    (args.out / "PO_OOD_STREAM_REPORT.md").write_text(md, encoding="utf-8")
    Path("docs/agod").mkdir(parents=True, exist_ok=True)
    Path("docs/agod/AGOD_po_ood_stream_multi.md").write_text(md, encoding="utf-8")
    print(md)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
