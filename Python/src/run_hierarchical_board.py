#!/usr/bin/env python3
"""Group rankings and a second hierarchical board (Diabetes readmission).

Rankings for a feature partition (nothing beyond drop-and-fit / permute-and-fit):

  rf_vimp     sum of RF impurity for W ~ X
  logo_auc    AUC drop after dropping the group and fitting again
  logo_r      R-risk rise after dropping the group and fitting again
  po_inflate  PO-risk rise after dropping the group and fitting again
  only_po     PO-risk using only that group
  perm_auc    AUC drop after permuting the group's columns (rows shuffled together)

Second dataset: source vs target DiabetesReadmission in datasets/datasets.zip.
Same three-layer board: group ranking → drop-and-fit PO-risk → high-|τ| subset.
"""
from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

VIDEO_DIM, AUDIO_DIM, TEXT_DIM = 768, 512, 768
MSRVTT_GROUPS = {
    "video": np.arange(0, VIDEO_DIM),
    "audio": np.arange(VIDEO_DIM, VIDEO_DIM + AUDIO_DIM),
    "text": np.arange(VIDEO_DIM + AUDIO_DIM, VIDEO_DIM + AUDIO_DIM + TEXT_DIM),
}

DIABETES_GROUPS = {
    "demographics": [
        "gender", "race_Caucasian", "race_AfricanAmerican", "race_Hispanic", "race_Asian", "age>=70",
    ],
    "utilization": [
        "time_in_hospital", "num_lab_procedures", "num_procedures", "num_medications",
        "number_outpatient", "number_emergency", "number_inpatient", "number_diagnoses",
    ],
    "labs": [
        "max_glu_serum", "A1Cresult", "max_glu_serum>200", "max_glu_serum>300", "A1Cresult>7", "A1Cresult>8",
    ],
    "meds": [
        "change", "diabetesMed", "metformin_Up", "metformin_Down", "metformin_Steady",
        "insulin_Up", "insulin_Down", "insulin_Steady",
    ],
    "admission": [
        "Transferred", "Home", "Emergency_admission", "Elective_admission", "Urgent_admission",
    ],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--feat-dir", default="/workspace/data/msrvtt/features_window")
    p.add_argument("--caption-json", default="/workspace/data/hf_ann/msrvtt_train_7k.json")
    p.add_argument("--msrvtt-json", default="experiments/msrvtt/concat_attribution.json")
    p.add_argument("--diabetes-zip", default="datasets/datasets.zip")
    p.add_argument("--out-dir", default="experiments/hierarchical")
    p.add_argument("--n-estimators", type=int, default=40)
    p.add_argument("--n-splits", type=int, default=3)
    p.add_argument("--max-depth", type=int, default=12)
    p.add_argument("--min-samples-leaf", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--clip-e", type=float, default=0.01)
    p.add_argument("--diabetes-per-batch", type=int, default=6000)
    p.add_argument("--skip-msrvtt-refit", action="store_true",
                   help="Reuse concat_attribution.json for MSR-VTT drop-and-fit numbers.")
    return p.parse_args()


def _rf_kwargs(seed, n_estimators, max_depth, min_samples_leaf):
    return dict(
        n_estimators=n_estimators, max_features="sqrt", max_depth=max_depth,
        min_samples_leaf=min_samples_leaf, random_state=seed, n_jobs=-1,
    )


def cv_auc(X, w, *, n_splits, seed, n_estimators, max_depth, min_samples_leaf):
    X, w = np.ascontiguousarray(X), np.asarray(w).astype(int)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for fold, (tr, te) in enumerate(skf.split(X, w)):
        rf = RandomForestClassifier(**_rf_kwargs(seed + fold, n_estimators, max_depth, min_samples_leaf))
        rf.fit(X[tr], w[tr])
        aucs.append(float(roc_auc_score(w[te], rf.predict_proba(X[te])[:, 1])))
    a = np.asarray(aucs)
    return float(a.mean())


def fit_po(X, Y, W, *, n_splits, seed, n_estimators, max_depth, min_samples_leaf, clip_e):
    X = np.ascontiguousarray(X)
    Y = np.asarray(Y, dtype=np.float64)
    W = np.asarray(W, dtype=int)
    n = X.shape[0]
    mu = np.zeros(n)
    e = np.zeros(n)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(skf.split(X, W)):
        kw = _rf_kwargs(seed + fold, n_estimators, max_depth, min_samples_leaf)
        mu[te] = RandomForestRegressor(**kw).fit(X[tr], Y[tr]).predict(X[te])
        e[te] = RandomForestClassifier(**kw).fit(X[tr], W[tr]).predict_proba(X[te])[:, 1]
    e = np.clip(e, clip_e, 1.0 - clip_e)
    ry, rt = Y - mu, W.astype(np.float64) - e
    phi = ry * rt
    tau_m = RandomForestRegressor(**_rf_kwargs(seed + 200, n_estimators, max_depth, min_samples_leaf))
    tau_m.fit(X, phi)
    tau = tau_m.predict(X)
    return {
        "po_risk": float(np.mean(tau ** 2)),
        "r_risk": float(np.mean((ry - tau * rt) ** 2)),
        "tau": tau,
        "phi": phi,
        "vimp": tau_m.feature_importances_.astype(np.float64),
    }


def ranks_from_scores(scores: dict[str, float], higher_is_more_important: bool) -> dict[str, int]:
    items = sorted(scores.items(), key=lambda kv: kv[1], reverse=higher_is_more_important)
    return {name: i + 1 for i, (name, _) in enumerate(items)}


def group_rankings(X, Y, W, groups: dict[str, np.ndarray], *, rfkw, pokw, seed: int, skip_po: bool = False) -> dict:
    """Several block rankings. Each is drop-and-fit or permute-and-fit."""
    names = list(groups)
    p = X.shape[1]
    rf = RandomForestClassifier(**_rf_kwargs(seed, rfkw["n_estimators"], rfkw["max_depth"], rfkw["min_samples_leaf"]))
    rf.fit(X, W)
    vimp = rf.feature_importances_
    rf_share = {}
    for g, idx in groups.items():
        rf_share[g] = float(np.clip(vimp[idx], 0, None).sum())
    s = sum(rf_share.values()) + 1e-12
    rf_share = {g: v / s for g, v in rf_share.items()}

    full_auc = cv_auc(X, W, **rfkw)
    logo_auc, perm_auc = {}, {}
    rng = np.random.RandomState(seed + 7)
    for g, idx in groups.items():
        keep = np.ones(p, dtype=bool)
        keep[idx] = False
        logo_auc[g] = full_auc - cv_auc(X[:, keep], W, **rfkw)
        Xp = X.copy()
        perm = rng.permutation(X.shape[0])
        Xp[:, idx] = X[perm][:, idx]
        perm_auc[g] = full_auc - cv_auc(Xp, W, **{**rfkw, "seed": rfkw["seed"] + 50})

    methods = {
        "rf_vimp": (rf_share, True),
        "logo_auc": (logo_auc, True),
        "perm_auc": (perm_auc, True),
    }
    full_po = None
    if not skip_po:
        full_po = fit_po(X, Y, W, **pokw)
        logo_r, po_inflate, only_po = {}, {}, {}
        for i, (g, idx) in enumerate(groups.items()):
            keep = np.ones(p, dtype=bool)
            keep[idx] = False
            dropped = fit_po(X[:, keep], Y, W, **{**pokw, "seed": pokw["seed"] + 10 * (i + 1)})
            only = fit_po(X[:, idx], Y, W, **{**pokw, "seed": pokw["seed"] + 20 * (i + 1)})
            logo_r[g] = dropped["r_risk"] - full_po["r_risk"]
            po_inflate[g] = dropped["po_risk"] - full_po["po_risk"]
            only_po[g] = only["po_risk"]
        methods["logo_r"] = (logo_r, True)
        methods["po_inflate"] = (po_inflate, True)
        methods["only_po"] = (only_po, True)

    scores = {m: sc for m, (sc, _) in methods.items()}
    ranks = {m: ranks_from_scores(sc, hi) for m, (sc, hi) in methods.items()}
    kendall = {}
    mnames = list(methods)
    for a in mnames:
        for b in mnames:
            ra = [ranks[a][g] for g in names]
            rb = [ranks[b][g] for g in names]
            tau, _ = kendalltau(ra, rb)
            kendall[f"{a}__{b}"] = float(tau) if tau == tau else 0.0
    out = {
        "groups": names,
        "scores": scores,
        "ranks": ranks,
        "kendall": kendall,
        "full_auc": full_auc,
    }
    if full_po is not None:
        out["top_group_logo_r"] = min(ranks["logo_r"], key=ranks["logo_r"].get)
        out["full_po_risk"] = full_po["po_risk"]
        out["full_r_risk"] = full_po["r_risk"]
        out["fit"] = full_po
    return out


def kendall_matrix(names, kendall):
    k = len(names)
    M = np.zeros((k, k))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            M[i, j] = kendall[f"{a}__{b}"]
    return M


def subset_from_tau(X, W, Y, tau, vimp, groups, top_g, *, top_k=8, top_frac=0.1, feature_names=None):
    idx = groups[top_g]
    tau_sq = tau ** 2
    n = len(tau)
    k = max(1, int(round(top_frac * n)))
    order = np.argsort(-tau_sq)
    local = vimp[idx]
    loc_order = np.argsort(-local)[:top_k]
    feats = []
    for rank, loc in enumerate(loc_order, start=1):
        j = int(idx[loc])
        name = feature_names[j] if feature_names is not None else f"x{j}"
        feats.append({
            "rank": rank, "index": j, "name": name, "vimp": float(local[loc]),
            "mean_W0": float(X[W == 0, j].mean()),
            "mean_W1": float(X[W == 1, j].mean()),
            "diff": float(X[W == 1, j].mean() - X[W == 0, j].mean()),
        })
    median_split = None
    if feats:
        j0 = feats[0]["index"]
        med = float(np.median(X[:, j0]))
        low, high = X[:, j0] <= med, X[:, j0] > med
        median_split = {
            "feature": feats[0]["name"],
            "median": med,
            "low": {"n": int(low.sum()), "mean_tau_sq": float(tau_sq[low].mean()), "frac_W1": float(W[low].mean()), "mean_Y": float(Y[low].mean())},
            "high": {"n": int(high.sum()), "mean_tau_sq": float(tau_sq[high].mean()), "frac_W1": float(W[high].mean()), "mean_Y": float(Y[high].mean())},
        }
    return {
        "top_group": top_g,
        "top_row_mass": float(tau_sq[order[:k]].sum() / (tau_sq.sum() + 1e-12)),
        "top_rows": [{"row": int(i), "W": int(W[i]), "Y": float(Y[i]), "tau_sq": float(tau_sq[i])} for i in order[: min(30, k)]],
        "top_features": feats,
        "median_split": median_split,
    }


def load_msrvtt(feat_dir, caption_json):
    feat_dir = Path(feat_dir)
    X = np.load(feat_dir / "concat_feat.npy").astype(np.float32)
    windows = np.load(feat_dir / "window_index.npy")
    labels = np.load(feat_dir / "video_labels.npy")
    meta = json.loads((feat_dir / "meta.json").read_text())
    n_windows = int(meta["n_windows"])
    video_ids = list(meta["video_ids"])
    cap = {rec["video_id"]: int(rec["category"]) for rec in json.loads(Path(caption_json).read_text())}
    Y = np.repeat(np.asarray([cap[v] for v in video_ids], dtype=np.float64), n_windows)
    W = (windows >= n_windows // 2).astype(int)
    return X, Y, W, MSRVTT_GROUPS, None


def load_diabetes(zip_path: Path, per_batch: int, seed: int):
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open("source_DiabetesReadmission.csv") as f:
            src = pd.read_csv(f)
        with zf.open("target_DiabetesReadmission.csv") as f:
            tgt = pd.read_csv(f)
    rng = np.random.RandomState(seed)
    if len(src) > per_batch:
        src = src.iloc[rng.choice(len(src), per_batch, replace=False)]
    if len(tgt) > per_batch:
        tgt = tgt.iloc[rng.choice(len(tgt), per_batch, replace=False)]
    df = pd.concat([src, tgt], axis=0, ignore_index=True)
    y = df["readmitted"].to_numpy(dtype=np.float64)
    w = np.concatenate([np.zeros(len(src), dtype=int), np.ones(len(tgt), dtype=int)])
    feat_cols = [c for c in df.columns if c != "readmitted"]
    X = df[feat_cols].to_numpy(dtype=np.float32)
    groups = {}
    for g, cols in DIABETES_GROUPS.items():
        groups[g] = np.asarray([feat_cols.index(c) for c in cols], dtype=int)
    return X, y, w, groups, feat_cols, {"n_source": int(len(src)), "n_target": int(len(tgt))}


def plot_rank_table(payload, out_png: Path, title: str):
    import matplotlib.pyplot as plt

    groups = payload["groups"]
    methods = list(payload["ranks"])
    cell = [[payload["ranks"][m][g] for g in groups] for m in methods]
    fig, ax = plt.subplots(figsize=(1.4 * (len(groups) + 2), 0.55 * (len(methods) + 3)), dpi=150)
    ax.axis("off")
    tbl = ax.table(
        cellText=cell,
        rowLabels=methods,
        colLabels=groups,
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.15, 1.35)
    for i, m in enumerate(methods):
        for j, g in enumerate(groups):
            if payload["ranks"][m][g] == 1:
                tbl[i + 1, j].set_facecolor("#F4C790")
    ax.set_title(title, fontsize=11, pad=18)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_diabetes_board(rank_payload, subset, out_png: Path):
    import matplotlib.pyplot as plt

    groups = rank_payload["groups"]
    x = np.arange(len(groups))
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.8), dpi=150)
    ax = axes[0]
    scores = rank_payload["scores"]["logo_r"]
    ax.bar(x, [scores[g] for g in groups], color="#2C4A6E")
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=25, ha="right")
    ax.set_ylabel("R-risk rise")
    ax.set_title("Drop group, fit again")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    methods = ["rf_vimp", "logo_auc", "logo_r", "po_inflate"]
    for i, g in enumerate(groups):
        ax.scatter(
            [i] * len(methods),
            [rank_payload["ranks"][m][g] for m in methods],
            c=["#E07A3D", "#2C4A6E", "#2F6B4F", "#7A4E8A"][: len(methods)],
            s=40,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=25, ha="right")
    ax.set_yticks([1, 2, 3, 4, 5])
    ax.invert_yaxis()
    ax.set_ylabel("rank (1 = top)")
    ax.set_title("Rankings agree or not")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[2]
    feats = subset["top_features"][:8]
    ax.barh(range(len(feats))[::-1], [abs(f["diff"]) for f in feats][::-1], color="#E07A3D")
    ax.set_yticks(range(len(feats)))
    ax.set_yticklabels([f["name"] for f in feats][::-1], fontsize=8)
    ax.set_xlabel("|mean W1 − W0|")
    ax.set_title(f"Subset coords in {subset['top_group']}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle("Diabetes readmission  ·  source vs target", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def json_safe(obj):
    skip = {"fit", "tau", "phi", "vimp"}
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items() if k not in skip}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def main():
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rfkw = dict(n_splits=args.n_splits, seed=args.seed, n_estimators=args.n_estimators,
                max_depth=args.max_depth, min_samples_leaf=args.min_samples_leaf)
    pokw = {**rfkw, "clip_e": args.clip_e}

    print("=== MSR-VTT early vs late: several group rankings ===", flush=True)
    X, Y, W, groups, _ = load_msrvtt(args.feat_dir, args.caption_json)
    json_path = Path(args.msrvtt_json)
    if json_path.exists():
        print("  rf_vimp / logo_auc / perm_auc on concat; PO drop-and-fit from saved json", flush=True)
        msrvtt = group_rankings(X, Y, W, groups, rfkw=rfkw, pokw=pokw, seed=args.seed, skip_po=True)
        prev = json.loads(json_path.read_text())["contrasts"]["pooled_temporal"]
        logo = prev["logo_po_risk"]
        msrvtt["scores"]["logo_r"] = {g: logo[f"without_{g}"]["delta_r_risk"] for g in groups}
        msrvtt["scores"]["po_inflate"] = {g: logo[f"without_{g}"]["po_risk"] - logo["full"]["po_risk"] for g in groups}
        msrvtt["scores"]["only_po"] = {g: logo[f"only_{g}"]["po_risk"] for g in groups}
        msrvtt["ranks"]["logo_r"] = ranks_from_scores(msrvtt["scores"]["logo_r"], True)
        msrvtt["ranks"]["po_inflate"] = ranks_from_scores(msrvtt["scores"]["po_inflate"], True)
        msrvtt["ranks"]["only_po"] = ranks_from_scores(msrvtt["scores"]["only_po"], True)
        msrvtt["top_group_logo_r"] = min(msrvtt["ranks"]["logo_r"], key=msrvtt["ranks"]["logo_r"].get)
        names = list(msrvtt["ranks"])
        kendall = {}
        for a in names:
            for b in names:
                ra = [msrvtt["ranks"][a][g] for g in groups]
                rb = [msrvtt["ranks"][b][g] for g in groups]
                tau, _ = kendalltau(ra, rb)
                kendall[f"{a}__{b}"] = float(tau) if tau == tau else 0.0
        msrvtt["kendall"] = kendall
    else:
        msrvtt = group_rankings(X, Y, W, groups, rfkw=rfkw, pokw=pokw, seed=args.seed)
        msrvtt.pop("fit", None)
    print("  ranks", json.dumps(msrvtt["ranks"], indent=2), flush=True)
    plot_rank_table(msrvtt, out / "msrvtt_group_ranks.png", "MSR-VTT group ranks (1 = top)")
    (out / "msrvtt_group_ranks.json").write_text(json.dumps(json_safe(msrvtt), indent=2))

    print("\n=== Diabetes source vs target: hierarchical board ===", flush=True)
    Xd, Yd, Wd, gdiab, feat_cols, meta = load_diabetes(Path(args.diabetes_zip), args.diabetes_per_batch, args.seed)
    print(f"  n={len(Yd)} p={Xd.shape[1]} source/target={meta}", flush=True)
    diab = group_rankings(Xd, Yd, Wd, gdiab, rfkw=rfkw, pokw=pokw, seed=args.seed)
    fit = diab.pop("fit")
    subset = subset_from_tau(
        Xd, Wd, Yd, fit["tau"], fit["vimp"], gdiab, diab["top_group_logo_r"],
        feature_names=feat_cols,
    )
    print("  ranks", json.dumps(diab["ranks"], indent=2), flush=True)
    print(f"  top group (logo_r) = {diab['top_group_logo_r']}", flush=True)
    print(f"  top features: {[f['name'] for f in subset['top_features'][:5]]}", flush=True)
    plot_rank_table(diab, out / "diabetes_group_ranks.png", "Diabetes group ranks (1 = top)")
    plot_diabetes_board(diab, subset, out / "diabetes_board.png")
    (out / "diabetes_hierarchical.json").write_text(json.dumps(json_safe({
        "meta": meta,
        "rankings": diab,
        "subset": subset,
    }), indent=2))
    print(f"\nwrote {out}", flush=True)


if __name__ == "__main__":
    main()
