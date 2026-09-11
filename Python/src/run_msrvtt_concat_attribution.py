#!/usr/bin/env python3
"""Attribution on concatenated MSR-VTT windows (3900 x 2048).

  1. RF VIMP modality shares.
  2. Leave-one-group-out AUC for batch label W.
  3. Leave-one-group-out PO-risk: drop a modality block, refit (μ, e, τ),
     recompute E[τ²] / R-risk.
  4. After the top modality is named, subset localization: high-|τ| videos
     and top coordinates in that block, with a simple W=0 vs W=1 mean split.

Two batches: pooled early vs late windows, and a 50/50 video mixture.
Y for PO-risk is the MSR-VTT category id.
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
GROUP_SEED = {"video": 1, "audio": 2, "text": 3}


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
    p.add_argument("--top-k-features", type=int, default=20)
    p.add_argument("--top-n-videos", type=int, default=12)
    p.add_argument("--top-row-frac", type=float, default=0.1)
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
        "categories": [int(cap[v]) for v in video_ids],
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


def fit_po_risk(X, Y, W, *, n_splits, seed, n_estimators, max_depth, min_samples_leaf, clip_e) -> dict:
    """Refit μ, e, τ on this X. PO-risk = mean(τ²); R-risk = mean((Y-μ-τ(W-e))²)."""
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
    return {
        "po_risk": float(np.mean(tau ** 2)),
        "r_risk": float(np.mean((residual_y - tau * residual_t) ** 2)),
        "phi": phi,
        "tau": tau,
        "vimp": tau_model.feature_importances_.astype(np.float64),
        "phi_mean": float(phi.mean()),
        "phi_std": float(phi.std()),
    }


def leave_one_group_out_po_risk(X, Y, W, **kw) -> dict:
    """Drop (or keep-only) a modality, refit the whole PO-risk pipeline."""
    full = fit_po_risk(X, Y, W, **kw)
    out = {
        "full": {"po_risk": full["po_risk"], "r_risk": full["r_risk"]},
        "fit": full,
    }
    for name, sl in GROUPS.items():
        keep = np.ones(X.shape[1], dtype=bool)
        keep[sl] = False
        dropped = fit_po_risk(X[:, keep], Y, W, **{**kw, "seed": kw["seed"] + 10 * GROUP_SEED[name]})
        only = fit_po_risk(X[:, sl], Y, W, **{**kw, "seed": kw["seed"] + 20 * GROUP_SEED[name]})
        out[f"without_{name}"] = {
            "po_risk": dropped["po_risk"],
            "r_risk": dropped["r_risk"],
            "delta_po_risk": full["po_risk"] - dropped["po_risk"],
            "delta_r_risk": dropped["r_risk"] - full["r_risk"],
        }
        out[f"only_{name}"] = {
            "po_risk": only["po_risk"],
            "r_risk": only["r_risk"],
        }
    # Drill into the group whose removal raises R-risk the most (LOCO).
    r_deltas = {g: out[f"without_{g}"]["delta_r_risk"] for g in GROUPS}
    po_inflate = {g: out[f"without_{g}"]["po_risk"] - out["full"]["po_risk"] for g in GROUPS}
    out["top_modality"] = max(r_deltas, key=r_deltas.get)
    out["rank_by_r_risk_loco"] = r_deltas
    out["rank_by_po_risk_inflate"] = po_inflate
    return out


def _feature_meta(j: int) -> tuple[str, int]:
    if j < VIDEO_DIM:
        return "video", int(j)
    if j < VIDEO_DIM + AUDIO_DIM:
        return "audio", int(j - VIDEO_DIM)
    return "text", int(j - VIDEO_DIM - AUDIO_DIM)


def subset_localize(
    X: np.ndarray,
    W: np.ndarray,
    Y: np.ndarray,
    tau: np.ndarray,
    phi: np.ndarray,
    vimp_full: np.ndarray,
    *,
    labels: np.ndarray,
    windows: np.ndarray,
    video_ids: list[str],
    categories: list[int],
    top_modality: str,
    top_k_features: int,
    top_n_videos: int,
    top_row_frac: float,
) -> dict:
    """Layer-2 subset after the modality is named: videos/windows and coordinates."""
    sl = GROUPS[top_modality]
    tau_sq = tau ** 2
    n = len(tau)
    k_rows = max(1, int(round(top_row_frac * n)))
    row_order = np.argsort(-tau_sq)

    top_rows = []
    for rank, i in enumerate(row_order[:k_rows], start=1):
        vid = video_ids[int(labels[i])]
        top_rows.append({
            "rank": rank,
            "row": int(i),
            "video_id": vid,
            "window": int(windows[i]),
            "W": int(W[i]),
            "Y": float(Y[i]),
            "tau": float(tau[i]),
            "tau_sq": float(tau_sq[i]),
            "phi": float(phi[i]),
        })

    per_video = []
    for vi, vid in enumerate(video_ids):
        mask = labels == vi
        per_video.append({
            "video_id": vid,
            "category": int(categories[vi]),
            "n": int(mask.sum()),
            "mean_tau_sq": float(tau_sq[mask].mean()),
            "mean_abs_phi": float(np.abs(phi[mask]).mean()),
            "mean_W": float(W[mask].mean()),
        })
    per_video.sort(key=lambda r: -r["mean_tau_sq"])

    # Coordinates inside the named modality, ranked by full-model VIMP on those columns.
    idx = np.arange(CONCAT_DIM)[sl]
    local_vimp = vimp_full[sl]
    order = np.argsort(-local_vimp)[:top_k_features]
    feat_rows = []
    for rank, loc in enumerate(order, start=1):
        j = int(idx[loc])
        x0 = X[W == 0, j]
        x1 = X[W == 1, j]
        feat_rows.append({
            "rank": rank,
            "index": j,
            "group": top_modality,
            "local_index": int(loc),
            "vimp": float(local_vimp[loc]),
            "mean_W0": float(x0.mean()) if x0.size else 0.0,
            "mean_W1": float(x1.mean()) if x1.size else 0.0,
            "diff_W1_minus_W0": float(x1.mean() - x0.mean()) if x0.size and x1.size else 0.0,
        })

    # Discretize the top coordinate in that modality (median split + quartiles).
    bins = None
    if feat_rows:
        j0 = feat_rows[0]["index"]
        col = X[:, j0]
        qs = np.quantile(col, [0.25, 0.5, 0.75])
        bin_id = np.digitize(col, qs, right=True)  # 0..3
        bins = []
        for b in range(4):
            m = bin_id == b
            bins.append({
                "bin": b,
                "n": int(m.sum()),
                "n_W1": int((W[m] == 1).sum()) if m.any() else 0,
                "mean_Y": float(Y[m].mean()) if m.any() else 0.0,
                "mean_tau_sq": float(tau_sq[m].mean()) if m.any() else 0.0,
                "frac_W1": float(W[m].mean()) if m.any() else 0.0,
            })
        median = float(qs[1])
        low, high = col <= median, col > median
        median_split = {
            "feature_index": int(j0),
            "median": median,
            "low": {
                "n": int(low.sum()),
                "mean_Y": float(Y[low].mean()) if low.any() else 0.0,
                "mean_tau_sq": float(tau_sq[low].mean()) if low.any() else 0.0,
                "frac_W1": float(W[low].mean()) if low.any() else 0.0,
            },
            "high": {
                "n": int(high.sum()),
                "mean_Y": float(Y[high].mean()) if high.any() else 0.0,
                "mean_tau_sq": float(tau_sq[high].mean()) if high.any() else 0.0,
                "frac_W1": float(W[high].mean()) if high.any() else 0.0,
            },
        }
    else:
        median_split = None

    mass = float(tau_sq.sum()) + 1e-12
    top_video_mass = float(sum(r["mean_tau_sq"] * r["n"] for r in per_video[:top_n_videos])) / mass
    top_row_mass = float(tau_sq[row_order[:k_rows]].sum()) / mass

    return {
        "top_modality": top_modality,
        "top_row_frac": top_row_frac,
        "top_row_mass": top_row_mass,
        "top_video_mass": top_video_mass,
        "top_rows": top_rows[: min(50, len(top_rows))],
        "per_video": per_video,
        "top_videos": per_video[:top_n_videos],
        "top_features_in_modality": feat_rows,
        "quartile_bins_top_feature": bins,
        "median_split_top_feature": median_split,
    }


def top_features(vimp: np.ndarray, k: int = 30) -> list[dict]:
    order = np.argsort(-vimp)[:k]
    rows = []
    for rank, j in enumerate(order, start=1):
        group, local = _feature_meta(int(j))
        rows.append({"rank": rank, "index": int(j), "group": group, "local_index": local, "vimp": float(vimp[j])})
    return rows


def plot_results(payload: dict, out_png: Path) -> None:
    import matplotlib.pyplot as plt

    contrasts = list(payload["contrasts"].keys())
    x = np.arange(len(contrasts))
    w = 0.18
    labels = ["early vs late", "video mix A vs B"]

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.6), dpi=150)

    def share_of(contrast, key):
        rec = payload["contrasts"][contrast]["fsds_shares"]
        return float(rec.get(f"{key}_mean", rec.get(key, 0.0)))

    ax = axes[0, 0]
    ww = 0.25
    ax.bar(x - ww, [share_of(c, "video_share") for c in contrasts], ww, color="#E07A3D", label="video")
    ax.bar(x, [share_of(c, "audio_share") for c in contrasts], ww, color="#2C4A6E", label="audio")
    ax.bar(x + ww, [share_of(c, "text_share") for c in contrasts], ww, color="#2F6B4F", label="text")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("block share")
    ax.set_title("RF VIMP shares")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[0, 1]
    full = [payload["contrasts"][c]["logo_po_risk"]["full"]["po_risk"] for c in contrasts]
    dv = [payload["contrasts"][c]["logo_po_risk"]["without_video"]["po_risk"] for c in contrasts]
    da = [payload["contrasts"][c]["logo_po_risk"]["without_audio"]["po_risk"] for c in contrasts]
    dt = [payload["contrasts"][c]["logo_po_risk"]["without_text"]["po_risk"] for c in contrasts]
    ax.bar(x - 1.5 * w, full, w, color="#1A2332", label="full")
    ax.bar(x - 0.5 * w, dv, w, color="#E07A3D", label="−video")
    ax.bar(x + 0.5 * w, da, w, color="#2C4A6E", label="−audio")
    ax.bar(x + 1.5 * w, dt, w, color="#2F6B4F", label="−text")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r"PO-risk  $E[\hat\tau^2]$")
    ax.set_title("Leave-one-group-out PO-risk")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1, 0]
    ov = [payload["contrasts"][c]["logo_po_risk"]["only_video"]["po_risk"] for c in contrasts]
    oa = [payload["contrasts"][c]["logo_po_risk"]["only_audio"]["po_risk"] for c in contrasts]
    ot = [payload["contrasts"][c]["logo_po_risk"]["only_text"]["po_risk"] for c in contrasts]
    ax.bar(x - ww, ov, ww, color="#E07A3D", label="video only")
    ax.bar(x, oa, ww, color="#2C4A6E", label="audio only")
    ax.bar(x + ww, ot, ww, color="#2F6B4F", label="text only")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r"PO-risk  $E[\hat\tau^2]$")
    ax.set_title("PO-risk using one group only")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1, 1]
    colors = ["#E07A3D", "#2C4A6E"]
    y_pos = None
    names = None
    # show subset videos for pooled_temporal (more within-clip structure)
    rec = payload["contrasts"]["pooled_temporal"]["subset"]
    names = [r["video_id"] for r in rec["top_videos"]][::-1]
    vals = [r["mean_tau_sq"] for r in rec["top_videos"]][::-1]
    y_pos = np.arange(len(names))
    ax.barh(y_pos, vals, color=colors[0])
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    top_m = rec["top_modality"]
    ax.set_xlabel(r"mean $\tau^2$")
    ax.set_title(f"Subset: top videos by |τ|  (early/late, after {top_m})")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("MSR-VTT concat 2048-d  ·  LOGO PO-risk then subset", fontsize=11)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def json_safe(obj):
    skip = {"feature_vimp", "phi", "tau", "vimp", "fit"}
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items() if k not in skip}
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
    pokw = {**rfkw, "clip_e": args.clip_e}
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
    tau_store = {}
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
        logo_auc = leave_one_group_out_auc(X, W, **rfkw)
        print(
            f"    full={logo_auc['full']['mean']:.3f}  "
            f"-video={logo_auc['without_video']['mean']:.3f} (Δ={logo_auc['delta_without_video']:.3f})  "
            f"-audio={logo_auc['without_audio']['mean']:.3f} (Δ={logo_auc['delta_without_audio']:.3f})  "
            f"-text={logo_auc['without_text']['mean']:.3f} (Δ={logo_auc['delta_without_text']:.3f})",
            flush=True,
        )
        print("  leave-one-group-out PO-risk (refit μ,e,τ)...", flush=True)
        logo_po = leave_one_group_out_po_risk(X, data["y_cat"], W, **pokw)
        full_fit = logo_po["fit"]
        print(
            f"    full PO-risk={logo_po['full']['po_risk']:.4g} R-risk={logo_po['full']['r_risk']:.4g}\n"
            f"    -video PO={logo_po['without_video']['po_risk']:.4g} "
            f"RΔ={logo_po['without_video']['delta_r_risk']:.4g}  "
            f"-audio PO={logo_po['without_audio']['po_risk']:.4g} "
            f"RΔ={logo_po['without_audio']['delta_r_risk']:.4g}  "
            f"-text PO={logo_po['without_text']['po_risk']:.4g} "
            f"RΔ={logo_po['without_text']['delta_r_risk']:.4g}",
            flush=True,
        )
        print(
            f"    only: video={logo_po['only_video']['po_risk']:.4g}  "
            f"audio={logo_po['only_audio']['po_risk']:.4g}  "
            f"text={logo_po['only_text']['po_risk']:.4g}  "
            f"top_modality={logo_po['top_modality']}",
            flush=True,
        )
        print(f"  subset localization inside {logo_po['top_modality']}...", flush=True)
        subset = subset_localize(
            X, W, data["y_cat"], full_fit["tau"], full_fit["phi"], full_fit["vimp"],
            labels=data["labels"], windows=data["windows"], video_ids=data["video_ids"],
            categories=data["categories"], top_modality=logo_po["top_modality"],
            top_k_features=args.top_k_features, top_n_videos=args.top_n_videos,
            top_row_frac=args.top_row_frac,
        )
        print(
            f"    top videos: {', '.join(r['video_id'] for r in subset['top_videos'][:5])}  "
            f"row-mass={subset['top_row_mass']:.3f}  video-mass={subset['top_video_mass']:.3f}",
            flush=True,
        )
        vimp_store[name] = full_fit["vimp"]
        tau_store[name] = full_fit["tau"]
        payload["contrasts"][name] = {
            "fsds_shares": shares,
            "logo_auc": logo_auc,
            "logo_po_risk": {k: v for k, v in logo_po.items() if k != "fit"},
            "po_risk": {
                "po_risk": full_fit["po_risk"],
                "r_risk": full_fit["r_risk"],
                "blockwise_shares": modality_shares(full_fit["vimp"]),
                "top_features": top_features(full_fit["vimp"], 30),
            },
            "subset": subset,
        }

    np.savez_compressed(
        out_dir / "po_risk_feature_vimp.npz",
        pooled_temporal=vimp_store["pooled_temporal"],
        mixture_shift=vimp_store["mixture_shift"],
        tau_pooled_temporal=tau_store["pooled_temporal"],
        tau_mixture_shift=tau_store["mixture_shift"],
    )
    (out_dir / "concat_attribution.json").write_text(json.dumps(json_safe(payload), indent=2))
    plot_results(payload, out_dir / "concat_attribution.png")
    print(f"\nwrote {out_dir / 'concat_attribution.json'}", flush=True)
    print(f"wrote {out_dir / 'concat_attribution.png'}", flush=True)
    print(f"wrote {out_dir / 'po_risk_feature_vimp.npz'}", flush=True)


if __name__ == "__main__":
    main()
