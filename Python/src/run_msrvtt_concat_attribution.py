#!/usr/bin/env python3
"""Run Layer-1 attribution on the already-concatenated MSR-VTT windows.

On (N_videos * 100, 2048) = video[768] ⊕ audio[512] ⊕ text[768]:

  1. RF VIMP modality shares (the FSDS block shares already used in extract).
  2. Leave-one-group-out AUC for predicting the batch label W.
  3. PO-risk (R-learner residual) attributed to features, then summed to blocks.

Two batches only: pooled early vs late windows, and a 50/50 video mixture.
Y for PO-risk is the MSR-VTT category id. Not online learning.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

VIDEO_DIM = 768
AUDIO_DIM = 512
TEXT_DIM = 768
CONCAT_DIM = VIDEO_DIM + AUDIO_DIM + TEXT_DIM
GROUPS = {
    "video": slice(0, VIDEO_DIM),
    "audio": slice(VIDEO_DIM, VIDEO_DIM + AUDIO_DIM),
    "text": slice(VIDEO_DIM + AUDIO_DIM, CONCAT_DIM),
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--feat-dir", default="/workspace/data/msrvtt/features_window")
    p.add_argument("--caption-json", default="/workspace/data/hf_ann/msrvtt_train_7k.json")
    p.add_argument("--out-dir", default="experiments/msrvtt")
    p.add_argument("--n-estimators", type=int, default=40)
    p.add_argument("--n-splits", type=int, default=3)
    p.add_argument("--n-bootstrap", type=int, default=3)
    p.add_argument("--max-depth", type=int, default=12)
    p.add_argument("--min-samples-leaf", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--clip-e", type=float, default=0.01)
    return p.parse_args()


def _rf_kwargs(seed: int, n_estimators: int, max_depth: int, min_samples_leaf: int) -> dict:
    return dict(
        n_estimators=n_estimators,
        max_features="sqrt",
        max_depth=max_depth,
        min_samples_leaf=min_samples_leaf,
        random_state=seed,
        n_jobs=-1,
    )


def modality_shares(vimp: np.ndarray) -> dict[str, float]:
    v = float(np.clip(vimp[:VIDEO_DIM], 0, None).sum())
    a = float(np.clip(vimp[VIDEO_DIM : VIDEO_DIM + AUDIO_DIM], 0, None).sum())
    t = float(np.clip(vimp[VIDEO_DIM + AUDIO_DIM :], 0, None).sum())
    total = v + a + t + 1e-12
    return {"video_share": v / total, "audio_share": a / total, "text_share": t / total}


def fit_rf_vimp(X: np.ndarray, y: np.ndarray, *, seed: int, n_estimators: int, max_depth: int, min_samples_leaf: int) -> np.ndarray:
    rf = RandomForestClassifier(**_rf_kwargs(seed, n_estimators, max_depth, min_samples_leaf))
    rf.fit(X, y)
    return rf.feature_importances_.astype(np.float64)


def bootstrap_shares(X, y, *, n_bootstrap, seed, n_estimators, max_depth, min_samples_leaf) -> dict:
    rng = np.random.RandomState(seed)
    rows = []
    b = 0
    attempts = 0
    n = len(y)
    while len(rows) < n_bootstrap and attempts < n_bootstrap * 5:
        attempts += 1
        idx = rng.randint(0, n, size=n)
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            continue
        vimp = fit_rf_vimp(
            X[idx], yb, seed=seed + b, n_estimators=n_estimators,
            max_depth=max_depth, min_samples_leaf=min_samples_leaf,
        )
        rows.append(modality_shares(vimp))
        b += 1
    if not rows:
        vimp = fit_rf_vimp(
            X, y, seed=seed, n_estimators=n_estimators,
            max_depth=max_depth, min_samples_leaf=min_samples_leaf,
        )
        rows.append(modality_shares(vimp))
    out = {"n_bootstrap": float(len(rows))}
    for k in ("video_share", "audio_share", "text_share"):
        arr = np.asarray([r[k] for r in rows], dtype=np.float64)
        mean = float(arr.mean())
        var = float(arr.var(ddof=1)) if arr.size > 1 else 0.0
        out[k] = mean
        out[f"{k}_mean"] = mean
        out[f"{k}_var"] = var
        out[f"{k}_std"] = float(np.sqrt(var))
    return out


def load_arrays(feat_dir: Path, caption_json: Path):
    feat_dir = Path(feat_dir)
    X = np.load(feat_dir / "concat_feat.npy").astype(np.float32)
    labels = np.load(feat_dir / "video_labels.npy").astype(int)
    windows = np.load(feat_dir / "window_index.npy").astype(int)
    meta = json.loads((feat_dir / "meta.json").read_text())
    video_ids = list(meta["video_ids"])
    n_windows = int(meta["n_windows"])
    n_videos = int(meta["n_videos"])
    if X.shape != (n_videos * n_windows, CONCAT_DIM):
        raise ValueError(f"concat shape {X.shape} != {(n_videos * n_windows, CONCAT_DIM)}")

    cap = {rec["video_id"]: int(rec["category"]) for rec in json.loads(Path(caption_json).read_text())}
    missing = [vid for vid in video_ids if vid not in cap]
    if missing:
        raise KeyError(f"missing categories for {missing[:8]}")
    y_cat = np.repeat(np.asarray([cap[v] for v in video_ids], dtype=np.float64), n_windows)

    mid_w = n_windows // 2
    w_temporal = (windows >= mid_w).astype(int)
    mid_v = n_videos // 2
    w_mixture = (labels >= mid_v).astype(int)
    return {
        "X": X,
        "labels": labels,
        "windows": windows,
        "video_ids": video_ids,
        "n_videos": n_videos,
        "n_windows": n_windows,
        "y_cat": y_cat,
        "contrasts": {
            "pooled_temporal": w_temporal,
            "mixture_shift": w_mixture,
        },
        "meta": meta,
    }


def cv_auc(X: np.ndarray, w: np.ndarray, *, n_splits: int, seed: int, n_estimators: int, max_depth: int, min_samples_leaf: int) -> dict:
    X = np.ascontiguousarray(X)
    w = np.asarray(w).astype(int)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for fold, (tr, te) in enumerate(skf.split(X, w)):
        rf = RandomForestClassifier(**_rf_kwargs(seed + fold, n_estimators, max_depth, min_samples_leaf))
        rf.fit(X[tr], w[tr])
        p = rf.predict_proba(X[te])[:, 1]
        aucs.append(float(roc_auc_score(w[te], p)))
    arr = np.asarray(aucs, dtype=np.float64)
    return {"mean": float(arr.mean()), "std": float(arr.std(ddof=1) if arr.size > 1 else 0.0), "folds": aucs}


def leave_one_group_out_auc(X, w, **rfkw) -> dict:
    full = cv_auc(X, w, **rfkw)
    out = {"full": full}
    for name, sl in GROUPS.items():
        keep = np.ones(X.shape[1], dtype=bool)
        keep[sl] = False
        dropped = cv_auc(X[:, keep], w, **rfkw)
        only = cv_auc(X[:, sl], w, **rfkw)
        out[f"without_{name}"] = dropped
        out[f"only_{name}"] = only
        out[f"delta_without_{name}"] = float(full["mean"] - dropped["mean"])
    return out


def crossfit_nuisances(X, Y, W, *, n_splits, seed, n_estimators, max_depth, min_samples_leaf, clip_e):
    n = X.shape[0]
    mu = np.zeros(n, dtype=np.float64)
    e = np.zeros(n, dtype=np.float64)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(skf.split(X, W)):
        kw = _rf_kwargs(seed + fold, n_estimators, max_depth, min_samples_leaf)
        mu[te] = RandomForestRegressor(**kw).fit(X[tr], Y[tr]).predict(X[te])
        e[te] = RandomForestClassifier(**kw).fit(X[tr], W[tr]).predict_proba(X[te])[:, 1]
    e = np.clip(e, clip_e, 1.0 - clip_e)
    return mu, e


def po_risk_attribution(X, Y, W, *, n_splits, seed, n_estimators, max_depth, min_samples_leaf, clip_e) -> dict:
    """R-learner φ = (Y-μ)(W-e); τ from RF on φ; impurity = feature-level attribution."""
    X = np.ascontiguousarray(X)
    Y = np.asarray(Y, dtype=np.float64)
    W = np.asarray(W, dtype=int)
    mu, e = crossfit_nuisances(
        X, Y, W, n_splits=n_splits, seed=seed, n_estimators=n_estimators,
        max_depth=max_depth, min_samples_leaf=min_samples_leaf, clip_e=clip_e,
    )
    residual_y = Y - mu
    residual_t = W.astype(np.float64) - e
    phi = residual_y * residual_t
    tau_model = RandomForestRegressor(**_rf_kwargs(seed + 200, n_estimators, max_depth, min_samples_leaf))
    tau_model.fit(X, phi)
    tau = tau_model.predict(X)
    po_risk = float(np.mean(tau ** 2))
    r_risk = float(np.mean((residual_y - tau * residual_t) ** 2))
    vimp = tau_model.feature_importances_.astype(np.float64)
    shares = modality_shares(vimp)

    logo_r_risk = {}
    logo_po_risk = {}
    for name, sl in GROUPS.items():
        keep = np.ones(X.shape[1], dtype=bool)
        keep[sl] = False
        tau_m = RandomForestRegressor(**_rf_kwargs(seed + 300 + {"video": 1, "audio": 2, "text": 3}[name], n_estimators, max_depth, min_samples_leaf))
        tau_m.fit(X[:, keep], phi)
        tau_hat = tau_m.predict(X[:, keep])
        rr = float(np.mean((residual_y - tau_hat * residual_t) ** 2))
        pr = float(np.mean(tau_hat ** 2))
        logo_r_risk[name] = {"r_risk": rr, "delta_vs_full": rr - r_risk}
        logo_po_risk[name] = {"po_risk": pr, "delta_vs_full": po_risk - pr}

    return {
        "po_risk": po_risk,
        "r_risk": r_risk,
        "phi_mean": float(phi.mean()),
        "phi_std": float(phi.std()),
        "blockwise_shares": shares,
        "leave_one_group_r_risk": logo_r_risk,
        "leave_one_group_po_risk": logo_po_risk,
        "feature_vimp": vimp,
    }


def top_features(vimp: np.ndarray, k: int = 30) -> list[dict]:
    order = np.argsort(-vimp)[:k]
    rows = []
    for rank, j in enumerate(order, start=1):
        if j < VIDEO_DIM:
            group, local = "video", int(j)
        elif j < VIDEO_DIM + AUDIO_DIM:
            group, local = "audio", int(j - VIDEO_DIM)
        else:
            group, local = "text", int(j - VIDEO_DIM - AUDIO_DIM)
        rows.append({"rank": rank, "index": int(j), "group": group, "local_index": local, "vimp": float(vimp[j])})
    return rows


def plot_results(payload: dict, out_png: Path) -> None:
    import matplotlib.pyplot as plt

    contrasts = list(payload["contrasts"].keys())
    x = np.arange(len(contrasts))
    w = 0.25

    def share_of(contrast, key):
        rec = payload["contrasts"][contrast]["fsds_shares"]
        return float(rec.get(f"{key}_mean", rec.get(key, 0.0)))

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.4), dpi=150)

    ax = axes[0]
    ax.bar(x - w, [share_of(c, "video_share") for c in contrasts], w, color="#E07A3D", label="video")
    ax.bar(x, [share_of(c, "audio_share") for c in contrasts], w, color="#2C4A6E", label="audio")
    ax.bar(x + w, [share_of(c, "text_share") for c in contrasts], w, color="#2F6B4F", label="text")
    ax.set_xticks(x)
    ax.set_xticklabels(["early vs late", "video mix A vs B"])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("block share")
    ax.set_title("RF VIMP shares")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    full = [payload["contrasts"][c]["logo_auc"]["full"]["mean"] for c in contrasts]
    ov = [payload["contrasts"][c]["logo_auc"]["only_video"]["mean"] for c in contrasts]
    oa = [payload["contrasts"][c]["logo_auc"]["only_audio"]["mean"] for c in contrasts]
    ot = [payload["contrasts"][c]["logo_auc"]["only_text"]["mean"] for c in contrasts]
    ax.bar(x - 1.5 * 0.2, full, 0.2, color="#1A2332", label="full")
    ax.bar(x - 0.5 * 0.2, ov, 0.2, color="#E07A3D", label="video only")
    ax.bar(x + 0.5 * 0.2, oa, 0.2, color="#2C4A6E", label="audio only")
    ax.bar(x + 1.5 * 0.2, ot, 0.2, color="#2F6B4F", label="text only")
    ax.axhline(0.5, color="#888", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(["early vs late", "video mix A vs B"])
    ax.set_ylim(0.0, 1.08)
    ax.set_ylabel("CV AUC")
    ax.set_title("AUC using one group only")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[2]
    def po_share(contrast, key):
        return float(payload["contrasts"][contrast]["po_risk"]["blockwise_shares"][key])
    ax.bar(x - w, [po_share(c, "video_share") for c in contrasts], w, color="#E07A3D", label="video")
    ax.bar(x, [po_share(c, "audio_share") for c in contrasts], w, color="#2C4A6E", label="audio")
    ax.bar(x + w, [po_share(c, "text_share") for c in contrasts], w, color="#2F6B4F", label="text")
    ax.set_xticks(x)
    ax.set_xticklabels(["early vs late", "video mix A vs B"])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("block share")
    ax.set_title("PO-risk feature VIMP (block sums)")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("MSR-VTT concat 2048-d  ·  39 clips × 100 windows", fontsize=11)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def json_safe(obj):
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items() if k != "feature_vimp"}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = load_arrays(Path(args.feat_dir), Path(args.caption_json))
    X = data["X"]
    print(
        f"concat {X.shape}  videos={data['n_videos']}  windows={data['n_windows']}",
        flush=True,
    )
    rfkw = dict(
        n_splits=args.n_splits,
        seed=args.seed,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        min_samples_leaf=args.min_samples_leaf,
    )
    payload = {
        "n_videos": data["n_videos"],
        "n_windows": data["n_windows"],
        "concat_shape": list(X.shape),
        "layout": "video[768] | audio[512] | text[768] = 2048",
        "y": "msrvtt category id",
        "n_estimators": args.n_estimators,
        "n_splits": args.n_splits,
        "n_bootstrap": args.n_bootstrap,
        "contrasts": {},
    }
    vimp_store = {}
    for name, W in data["contrasts"].items():
        print(f"\n=== {name}  W mean={float(W.mean()):.3f} ===", flush=True)
        print("  FSDS RF VIMP shares...", flush=True)
        shares = bootstrap_shares(
            X, W, n_bootstrap=args.n_bootstrap, seed=args.seed,
            n_estimators=args.n_estimators, max_depth=args.max_depth,
            min_samples_leaf=args.min_samples_leaf,
        )
        print(
            f"    video={shares['video_share']:.3f}  audio={shares['audio_share']:.3f}  "
            f"text={shares['text_share']:.3f}",
            flush=True,
        )
        print("  leave-one-group-out AUC...", flush=True)
        logo = leave_one_group_out_auc(X, W, **rfkw)
        print(
            f"    full={logo['full']['mean']:.3f}  "
            f"-video={logo['without_video']['mean']:.3f} (Δ={logo['delta_without_video']:.3f})  "
            f"-audio={logo['without_audio']['mean']:.3f} (Δ={logo['delta_without_audio']:.3f})  "
            f"-text={logo['without_text']['mean']:.3f} (Δ={logo['delta_without_text']:.3f})",
            flush=True,
        )
        print("  PO-risk feature attribution...", flush=True)
        po = po_risk_attribution(
            X, data["y_cat"], W, clip_e=args.clip_e, **rfkw,
        )
        print(
            f"    PO-risk={po['po_risk']:.4g}  R-risk={po['r_risk']:.4g}  "
            f"shares video={po['blockwise_shares']['video_share']:.3f}  "
            f"audio={po['blockwise_shares']['audio_share']:.3f}  "
            f"text={po['blockwise_shares']['text_share']:.3f}",
            flush=True,
        )
        vimp_store[name] = po["feature_vimp"]
        payload["contrasts"][name] = {
            "fsds_shares": shares,
            "logo_auc": logo,
            "po_risk": {
                **{k: v for k, v in po.items() if k != "feature_vimp"},
                "top_features": top_features(po["feature_vimp"], 30),
            },
        }

    np.savez_compressed(
        out_dir / "po_risk_feature_vimp.npz",
        pooled_temporal=vimp_store["pooled_temporal"],
        mixture_shift=vimp_store["mixture_shift"],
    )
    (out_dir / "concat_attribution.json").write_text(json.dumps(json_safe(payload), indent=2))
    plot_results(payload, out_dir / "concat_attribution.png")
    print(f"\nwrote {out_dir / 'concat_attribution.json'}", flush=True)
    print(f"wrote {out_dir / 'concat_attribution.png'}", flush=True)
    print(f"wrote {out_dir / 'po_risk_feature_vimp.npz'}", flush=True)


if __name__ == "__main__":
    main()
