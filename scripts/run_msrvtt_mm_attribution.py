#!/usr/bin/env python3
"""MSR-VTT multimodal attribution board.

Packed features: [video 768 | audio 512 | text 768], W=labelsmsr ∈ {0,1}.

Delivers:
  1) Leave-one-group-out (LOGO) PO-risk contribution by modality
  2) Feature-specific PO / RF risk contribution (τ-VIMP & RF-VIMP)
  3) Video/sample-granularity variation via block bootstrap
  4) RF-Domain (AUC/VIMP) + MMD-LOGO (group VIMP) ~ PO-risk LOGO
  5) Resampling SI: bootstrap CIs that modality shares differ
  6) Per-subset / per-modality AUC visualization + LaTeX tables

  python3 scripts/run_msrvtt_mm_attribution.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "msrvtt" / "packed"
OUT = ROOT / "results" / "msrvtt_mm_attr"
DOCS = ROOT / "docs" / "method"

SEED = 2026
ALPHA_CLIP = 0.05
N_BOOT = 10  # mean/std over 10 resamples (keep light)
N_SUB = 2000
MMD_MAX_N = 300
RF_TREES = 160
PO_TREES = 60

BLOCKS = {
    "video": (0, 768),
    "audio": (768, 768 + 512),
    "text": (768 + 512, 768 + 512 + 768),
}
MODS = ["video", "audio", "text"]


def load_pack():
    V = np.load(DATA / "video_feat.npy").astype(np.float64)
    A = np.load(DATA / "audio_feat.npy").astype(np.float64)
    T = np.load(DATA / "text_feat.npy").astype(np.float64)
    W = np.load(DATA / "labelsmsr.npy").astype(int)
    assert V.shape == (len(W), 768) and A.shape[1] == 512 and T.shape[1] == 768
    X = np.hstack([V, A, T])
    # continuous outcome for PO path: PC1 of X, rank-quantile → [0,1]
    pc1 = PCA(n_components=1, random_state=SEED).fit_transform(X).ravel()
    order = np.argsort(pc1, kind="mergesort")
    ranks = np.empty(len(pc1), float)
    ranks[order] = np.linspace(0.0, 1.0, len(pc1))
    Y = ranks
    return X, Y, W


def subsample(X, Y, W, *, n_max=N_SUB, seed=SEED):
    rng = np.random.default_rng(seed)
    i0, i1 = np.where(W == 0)[0], np.where(W == 1)[0]
    n0 = min(len(i0), n_max // 2)
    n1 = min(len(i1), n_max // 2)
    take = np.concatenate([rng.choice(i0, n0, False), rng.choice(i1, n1, False)])
    rng.shuffle(take)
    return X[take], Y[take], W[take], take


def rf_domain(X, W, *, seed=SEED):
    depth = max(3, int(round(np.sqrt(X.shape[1]))))
    leaf = max(1, int(round(np.sqrt(len(X)) // 2)))
    clf = RandomForestClassifier(
        n_estimators=RF_TREES,
        max_depth=depth,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    Xtr, Xte, Wtr, Wte = train_test_split(
        X, W, test_size=0.25, random_state=seed, stratify=W
    )
    clf2 = RandomForestClassifier(
        n_estimators=120,
        max_depth=depth,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed + 1,
    )
    clf2.fit(Xtr, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return vimp, auc


def modality_mass(vimp: np.ndarray) -> dict:
    mass = {}
    for m, (a, b) in BLOCKS.items():
        mass[m] = float(vimp[a:b].sum())
    s = sum(mass.values()) + 1e-12
    return {m: mass[m] / s for m in MODS}


def po_risk_fit(X, Y, W, *, seed=SEED):
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    Yf = Y.astype(float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=PO_TREES,
            max_depth=8,
            min_samples_leaf=4,
            random_state=seed + fold,
            n_jobs=-1,
        )
        e = RandomForestClassifier(
            n_estimators=PO_TREES,
            max_depth=8,
            min_samples_leaf=4,
            random_state=seed + 40 + fold,
            n_jobs=-1,
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1 - ALPHA_CLIP)
    po = (Yf - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=120,
        max_depth=10,
        min_samples_leaf=4,
        random_state=seed + 7,
        n_jobs=-1,
    )
    tau.fit(X, po)
    tau_hat = tau.predict(X)
    risk = float(np.mean(tau_hat**2))
    return dict(risk=risk, tau=tau, tau_hat=tau_hat, po=po, vimp=tau.feature_importances_.astype(float))


def drop_block(X, mod: str) -> np.ndarray:
    a, b = BLOCKS[mod]
    return np.hstack([X[:, :a], X[:, b:]]) if a > 0 else X[:, b:]


def logo_po_contributions(X, Y, W, *, seed=SEED):
    full = po_risk_fit(X, Y, W, seed=seed)
    out = {"full_po_risk": full["risk"], "delta": {}, "rel_delta": {}, "tau_vimp_mass": modality_mass(full["vimp"])}
    for m in MODS:
        Xm = drop_block(X, m)
        rm = po_risk_fit(Xm, Y, W, seed=seed + 11 + MODS.index(m))["risk"]
        d = full["risk"] - rm
        out["delta"][m] = float(d)
        out["rel_delta"][m] = float(d / (full["risk"] + 1e-12))
    # normalize positive parts as share (for SI compare)
    pos = {m: max(out["delta"][m], 0.0) for m in MODS}
    s = sum(pos.values()) + 1e-12
    out["share"] = {m: pos[m] / s for m in MODS}
    out["feature_vimp_top"] = []
    order = np.argsort(-full["vimp"])
    for r, j in enumerate(order[:15]):
        j = int(j)
        mod = "video" if j < 768 else ("audio" if j < 768 + 512 else "text")
        local = j if j < 768 else (j - 768 if j < 768 + 512 else j - 768 - 512)
        out["feature_vimp_top"].append(
            {"rank": r + 1, "global_index": j, "modality": mod, "local_index": local, "vimp": float(full["vimp"][j])}
        )
    return out, full


def rbf_mmd2(X0, X1, *, max_n=MMD_MAX_N, rng=None):
    rng = rng or np.random.default_rng(0)
    if len(X0) > max_n:
        X0 = X0[rng.choice(len(X0), max_n, False)]
    if len(X1) > max_n:
        X1 = X1[rng.choice(len(X1), max_n, False)]
    Z = np.vstack([X0, X1])
    # median heuristic on subsample
    s = Z[rng.choice(len(Z), min(400, len(Z)), False)]
    d2 = ((s[:, None, :] - s[None, :, :]) ** 2).sum(-1)
    med = np.median(d2[np.triu_indices_from(d2, 1)])
    gamma = 1.0 / max(med, 1e-8)
    def k(A, B):
        return np.exp(-gamma * ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))
    K00, K11, K01 = k(X0, X0), k(X1, X1), k(X0, X1)
    n0, n1 = len(X0), len(X1)
    return float(
        (K00.sum() - np.trace(K00)) / max(n0 * (n0 - 1), 1)
        + (K11.sum() - np.trace(K11)) / max(n1 * (n1 - 1), 1)
        - 2.0 * K01.mean()
    )


def mmd_logo(X, W, *, seed=SEED):
    rng = np.random.default_rng(seed)
    X0, X1 = X[W == 0], X[W == 1]
    full = rbf_mmd2(X0, X1, rng=rng)
    delta = {}
    for m in MODS:
        Xm = drop_block(X, m)
        d = full - rbf_mmd2(Xm[W == 0], Xm[W == 1], rng=rng)
        delta[m] = float(d)
    pos = {m: max(delta[m], 0.0) for m in MODS}
    s = sum(pos.values()) + 1e-12
    share = {m: pos[m] / s for m in MODS}
    return {"mmd_full": float(full), "delta": delta, "share": share}


def modality_only_auc(X, W, *, seed=SEED):
    out = {}
    for m, (a, b) in BLOCKS.items():
        _, auc = rf_domain(X[:, a:b], W, seed=seed + MODS.index(m))
        out[m] = auc
    return out


def bootstrap_si(X, Y, W, *, n_boot=N_BOOT, seed=SEED):
    """Resample rows; recompute RF/PO/MMD modality shares; CI on pairwise diffs."""
    rng = np.random.default_rng(seed)
    keys = ["rf", "po", "mmd"]
    stores = {k: {m: [] for m in MODS} for k in keys}
    aucs = []
    mod_aucs = {m: [] for m in MODS}
    for b in range(n_boot):
        # stratified bootstrap
        i0, i1 = np.where(W == 0)[0], np.where(W == 1)[0]
        take = np.concatenate(
            [rng.choice(i0, len(i0), True), rng.choice(i1, len(i1), True)]
        )
        Xb, Yb, Wb = X[take], Y[take], W[take]
        # light subsample for speed inside bootstrap
        Xs, Ys, Ws, _ = subsample(Xb, Yb, Wb, n_max=1600, seed=seed + b)
        if len(np.unique(Ws)) < 2:
            continue
        vimp, auc = rf_domain(Xs, Ws, seed=seed + b)
        rf_mass = modality_mass(vimp)
        for m in MODS:
            stores["rf"][m].append(rf_mass[m])
        aucs.append(auc)
        ma = modality_only_auc(Xs, Ws, seed=seed + 3 + b)
        for m in MODS:
            mod_aucs[m].append(ma[m])
        # PO share via τ-VIMP mass (fast); LOGO reserved for point estimate
        po_fit = po_risk_fit(Xs, Ys, Ws, seed=seed + 5 + b)
        po_mass = modality_mass(po_fit["vimp"])
        for m in MODS:
            stores["po"][m].append(po_mass[m])
        mmd = mmd_logo(Xs, Ws, seed=seed + 7 + b)
        for m in MODS:
            stores["mmd"][m].append(mmd["share"][m])
        if (b + 1) % 5 == 0:
            print(f"  bootstrap {b+1}/{n_boot}", flush=True)

    def summarize(arr):
        a = np.asarray(arr, float)
        return {
            "mean": float(a.mean()),
            "std": float(a.std()),
            "ci95": [float(np.quantile(a, 0.025)), float(np.quantile(a, 0.975))],
        }

    share_sum = {}
    for k in keys:
        share_sum[k] = {m: summarize(stores[k][m]) for m in MODS}
    # pairwise differences
    diffs = {}
    for k in keys:
        for i, m1 in enumerate(MODS):
            for m2 in MODS[i + 1 :]:
                d = np.asarray(stores[k][m1]) - np.asarray(stores[k][m2])
                lo, hi = float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))
                diffs[f"{k}:{m1}-{m2}"] = {
                    "mean": float(d.mean()),
                    "std": float(d.std()),
                    "ci95": [lo, hi],
                    "rejects_equal_0": bool(abs(float(d.mean())) > 2 * float(d.std() + 1e-12)),
                }
    return {
        "n_boot": n_boot,
        "share": share_sum,
        "pairwise_diff": diffs,
        "auc": summarize(aucs),
        "modality_auc": {m: summarize(mod_aucs[m]) for m in MODS},
        "raw_aucs": [float(x) for x in aucs],
        "raw_modality_aucs": {m: [float(x) for x in mod_aucs[m]] for m in MODS},
    }


def video_granularity_variation(X, Y, W, *, n_groups=25, seed=SEED):
    """Partition samples into contiguous groups (proxy video blocks); per-group RF mass."""
    rng = np.random.default_rng(seed)
    n = len(W)
    # shuffle then chunk — variation across subsets
    perm = rng.permutation(n)
    sizes = np.full(n_groups, n // n_groups)
    sizes[: n % n_groups] += 1
    rows = []
    start = 0
    for g, sz in enumerate(sizes):
        idx = perm[start : start + sz]
        start += sz
        if len(np.unique(W[idx])) < 2 or sz < 80:
            continue
        Xs, Ys, Ws, _ = subsample(X[idx], Y[idx], W[idx], n_max=min(800, sz), seed=seed + g)
        if len(np.unique(Ws)) < 2:
            continue
        vimp, auc = rf_domain(Xs, Ws, seed=seed + g)
        mass = modality_mass(vimp)
        rows.append({"group": g, "n": int(len(idx)), "auc": auc, **mass})
    df = pd.DataFrame(rows)
    summary = {
        "n_groups": int(len(df)),
        "mean": {m: float(df[m].mean()) for m in MODS} if len(df) else {},
        "std": {m: float(df[m].std()) for m in MODS} if len(df) else {},
        "auc_mean": float(df["auc"].mean()) if len(df) else None,
        "auc_std": float(df["auc"].std()) if len(df) else None,
        "groups": rows,
    }
    return summary


def make_plots(point, boot, gran, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    # 1) modality share comparison
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), dpi=140)
    methods = [("rf", "RF-Domain VIMP"), ("po", "PO-risk LOGO"), ("mmd", "MMD-LOGO")]
    colors = {"video": "#2B6CB0", "audio": "#C05621", "text": "#276749"}
    for ax, (key, title) in zip(axes, methods):
        means = [point[key]["share"][m] for m in MODS]
        if key in boot["share"]:
            stds = [boot["share"][key][m]["std"] for m in MODS]
            yerr = np.asarray(stds)
        else:
            yerr = None
        ax.bar(MODS, means, color=[colors[m] for m in MODS], yerr=yerr, capsize=4, ecolor="#333")
        ax.set_ylim(0, 1)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("modality share")
    fig.suptitle("MSR-VTT multimodal contribution (bootstrap mean±std, B=10)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "modality_share_compare.png", bbox_inches="tight")
    plt.close(fig)

    # 2) AUC distributions
    fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=140)
    ax.hist(boot["raw_aucs"], bins=min(10, max(5, len(boot["raw_aucs"]))), color="#4A5568", alpha=0.85, edgecolor="white")
    ax.axvline(point["rf_auc"], color="#E53E3E", lw=2, label=f"point AUC={point['rf_auc']:.3f}")
    mu, sd = boot["auc"]["mean"], boot["auc"]["std"]
    ax.axvline(mu, color="#DD6B20", ls="--", label=f"boot {mu:.3f}±{sd:.3f}")
    ax.set_xlabel("RF-Domain AUC")
    ax.set_ylabel("bootstrap count")
    ax.set_title("Resampled RF-Domain AUC (B=10)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "rf_auc_bootstrap.png", bbox_inches="tight")
    plt.close(fig)

    # 3) per-modality AUC
    fig, ax = plt.subplots(figsize=(6.5, 3.8), dpi=140)
    means = [boot["modality_auc"][m]["mean"] for m in MODS]
    stds = [boot["modality_auc"][m]["std"] for m in MODS]
    ax.bar(MODS, means, color=[colors[m] for m in MODS], yerr=stds, capsize=4)
    # point estimates
    for i, m in enumerate(MODS):
        ax.scatter([i], [point["modality_auc"][m]], color="black", zorder=3, s=28)
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel("RF AUC (modality alone)")
    ax.set_title("Per-modality subset AUC (point + boot mean±std)")
    fig.tight_layout()
    fig.savefig(out_dir / "modality_subset_auc.png", bbox_inches="tight")
    plt.close(fig)

    # 4) video/group granularity variation
    if gran["groups"]:
        df = pd.DataFrame(gran["groups"])
        fig, ax = plt.subplots(figsize=(7.5, 3.8), dpi=140)
        x = np.arange(len(df))
        ax.plot(x, df["video"], label="video", color=colors["video"], lw=1.5)
        ax.plot(x, df["audio"], label="audio", color=colors["audio"], lw=1.5)
        ax.plot(x, df["text"], label="text", color=colors["text"], lw=1.5)
        ax.set_xlabel("subset / video-proxy group")
        ax.set_ylabel("RF modality share")
        ax.set_title("Video-granularity variation of modality mass")
        ax.legend(fontsize=8)
        ax.set_ylim(0, 1)
        fig.tight_layout()
        fig.savefig(out_dir / "video_granularity_mass.png", bbox_inches="tight")
        plt.close(fig)


def write_latex(point, boot, gran, out_paths):
    def esc(s):
        return str(s).replace("_", "\\_")

    lines = [
        "% MSR-VTT multimodal attribution — RF / MMD-LOGO / PO-risk LOGO + bootstrap SI\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{MSR-VTT multimodal attribution on packed features "
        "$X=[\\mathrm{video}_{768}\\,\\|\\,\\mathrm{audio}_{512}\\,\\|\\,\\mathrm{text}_{768}]$, "
        "$W=\\texttt{labelsmsr}\\in\\{0,1\\}$. "
        "Point estimates with bootstrap mean$\\pm$std ($B=" + str(boot["n_boot"]) + "$).}\n",
        "\\label{tab:msrvtt-mm-share}\n\\small\n",
        "\\begin{tabular}{@{}l ccc ccc@{}}\n\\toprule\n",
        "Method & \\multicolumn{3}{c}{Modality share (point)} & \\multicolumn{3}{c}{Boot mean$\\pm$std} \\\\\n",
        "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
        " & video & audio & text & video & audio & text \\\\\n\\midrule\n",
    ]
    for key, name in [("rf", "RF-Domain VIMP"), ("po", "PO-risk LOGO"), ("mmd", "MMD-LOGO")]:
        sh = point[key]["share"]
        ci = boot["share"][key]
        lines.append(
            f"{name} & ${sh['video']:.3f}$ & ${sh['audio']:.3f}$ & ${sh['text']:.3f}$ & "
            f"${ci['video']['mean']:.3f}\\pm{ci['video']['std']:.3f}$ & "
            f"${ci['audio']['mean']:.3f}\\pm{ci['audio']['std']:.3f}$ & "
            f"${ci['text']['mean']:.3f}\\pm{ci['text']['std']:.3f}$ \\\\\n"
        )
    lines += [
        "\\midrule\n",
        f"RF-Domain AUC & \\multicolumn{{6}}{{c}}{{${point['rf_auc']:.3f}$ "
        f"(boot ${boot['auc']['mean']:.3f}\\pm{boot['auc']['std']:.3f}$)}} \\\\\n",
        f"PO-risk (full) & \\multicolumn{{6}}{{c}}{{${point['po']['full_po_risk']:.6f}$}} \\\\\n",
        f"MMD$^2$ (full) & \\multicolumn{{6}}{{c}}{{${point['mmd']['mmd_full']:.6f}$}} \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
    ]

    lines += [
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Leave-one-group-out PO-risk contributions "
        "$\\Delta_g=\\widehat R-\\widehat R^{(-g)}$ and relative $\\Delta_g/\\widehat R$.}\n",
        "\\label{tab:msrvtt-logo-po}\n\\small\n",
        "\\begin{tabular}{@{}l c c c@{}}\n\\toprule\n",
        "Modality & $\\Delta$ PO-risk & Rel. $\\Delta$ & Share \\\\\n\\midrule\n",
    ]
    for m in MODS:
        lines.append(
            f"{m} & ${point['po']['delta'][m]:.6f}$ & ${point['po']['rel_delta'][m]:.3f}$ & "
            f"${point['po']['share'][m]:.3f}$ \\\\\n"
        )
    lines += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n\n"]

    lines += [
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Statistical inference on modality-share differences "
        "(bootstrap $B=" + str(boot["n_boot"]) + "$, mean$\\pm$std of diffs; "
        "reject equal if $|\\mathrm{mean}| > 2\\,\\mathrm{std}$).}\n",
        "\\label{tab:msrvtt-si-diff}\n\\small\n",
        "\\begin{tabular}{@{}l l c c c@{}}\n\\toprule\n",
        "Method & Contrast & Mean diff & Std & Reject equal \\\\\n\\midrule\n",
    ]
    for k, v in boot["pairwise_diff"].items():
        method, contrast = k.split(":", 1)
        method = {"rf": "RF", "po": "PO-LOGO", "mmd": "MMD-LOGO"}[method]
        # also attach std if present
        std = v.get("std", abs(v["ci95"][1] - v["ci95"][0]) / 4.0)
        rej = "yes" if abs(v["mean"]) > 2 * max(std, 1e-12) else "no"
        lines.append(
            f"{method} & {esc(contrast)} & ${v['mean']:.3f}$ & "
            f"${std:.3f}$ & {rej} \\\\\n"
        )
    lines += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n\n"]

    lines += [
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Per-modality RF-Domain AUC when using only that block "
        "(point estimate and bootstrap mean$\\pm$std, $B=" + str(boot["n_boot"]) + "$).}\n",
        "\\label{tab:msrvtt-mod-auc}\n\\small\n",
        "\\begin{tabular}{@{}l c c@{}}\n\\toprule\n",
        "Modality & Point AUC & Boot mean$\\pm$std \\\\\n\\midrule\n",
    ]
    for m in MODS:
        b = boot["modality_auc"][m]
        lines.append(
            f"{m} & ${point['modality_auc'][m]:.3f}$ & "
            f"${b['mean']:.3f}\\pm{b['std']:.3f}$ \\\\\n"
        )
    lines += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n\n"]

    lines += [
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Video-/subset-granularity variation of RF modality mass "
        "(n\\_groups=" + str(gran["n_groups"]) + ").}\n",
        "\\label{tab:msrvtt-granularity}\n\\small\n",
        "\\begin{tabular}{@{}l c c@{}}\n\\toprule\n",
        "Modality & Mean share & Std across groups \\\\\n\\midrule\n",
    ]
    for m in MODS:
        lines.append(
            f"{m} & ${gran['mean'].get(m, float('nan')):.3f}$ & "
            f"${gran['std'].get(m, float('nan')):.3f}$ \\\\\n"
        )
    lines += [
        "\\midrule\n",
        f"Group AUC & ${gran['auc_mean']:.3f}$ & ${gran['auc_std']:.3f}$ \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
    ]

    # feature-specific top table
    lines += [
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Feature-specific PO $\\widehat\\tau$-VIMP top ranks "
        "(global index in concat $X$).}\n",
        "\\label{tab:msrvtt-feat-vimp}\n\\small\n",
        "\\begin{tabular}{@{}r l r c@{}}\n\\toprule\n",
        "Rank & Modality & Local idx & VIMP \\\\\n\\midrule\n",
    ]
    for r in point["po"]["feature_vimp_top"][:10]:
        lines.append(
            f"{r['rank']} & {r['modality']} & {r['local_index']} & ${r['vimp']:.5f}$ \\\\\n"
        )
    lines += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n"]

    text = "".join(lines)
    for p in out_paths:
        p.write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    print("loading pack…", flush=True)
    X, Y, W = load_pack()
    print(f"X={X.shape} Y∈[{Y.min():.3f},{Y.max():.3f}] W={np.bincount(W)}", flush=True)

    Xs, Ys, Ws, _ = subsample(X, Y, W, n_max=N_SUB, seed=SEED)
    print(f"fit subsample n={len(Ws)}", flush=True)

    print("RF-Domain…", flush=True)
    rf_vimp, rf_auc = rf_domain(Xs, Ws, seed=SEED)
    rf_share = modality_mass(rf_vimp)
    print(f"  AUC={rf_auc:.3f} share={rf_share}", flush=True)

    print("PO-risk LOGO…", flush=True)
    po, po_full = logo_po_contributions(Xs, Ys, Ws, seed=SEED)
    print(f"  risk={po['full_po_risk']:.6f} share={po['share']}", flush=True)

    print("MMD-LOGO…", flush=True)
    mmd = mmd_logo(Xs, Ws, seed=SEED)
    print(f"  mmd2={mmd['mmd_full']:.6f} share={mmd['share']}", flush=True)

    print("per-modality AUC…", flush=True)
    mod_auc = modality_only_auc(Xs, Ws, seed=SEED)
    print(mod_auc, flush=True)

    print("video/subset granularity…", flush=True)
    gran = video_granularity_variation(X, Y, W, n_groups=25, seed=SEED)
    print(f"  groups={gran['n_groups']} mean_share={gran['mean']}", flush=True)

    print(f"bootstrap SI B={N_BOOT}…", flush=True)
    boot = bootstrap_si(X, Y, W, n_boot=N_BOOT, seed=SEED)

    point = {
        "rf_auc": rf_auc,
        "rf": {"share": rf_share, "vimp_mass_raw": {m: float(rf_vimp[a:b].sum()) for m, (a, b) in BLOCKS.items()}},
        "po": po,
        "mmd": mmd,
        "modality_auc": mod_auc,
        "n": int(len(W)),
        "n_fit": int(len(Ws)),
        "layout": "video[0:768]|audio[768:1280]|text[1280:2048]|W=labelsmsr",
        "y_def": "rank-quantile of PCA-PC1(X) → [0,1]",
    }

    make_plots(point, boot, gran, OUT)
    write_latex(
        point,
        boot,
        gran,
        [OUT / "MSR_VTT_Multimodal_Attribution_tables_only.tex", DOCS / "MSR_VTT_Multimodal_Attribution_tables_only.tex"],
    )

    payload = {"point": point, "bootstrap": boot, "granularity": gran}
    # shrink raw auc lists already small
    (OUT / "msrvtt_mm_attribution_board.json").write_text(json.dumps(payload, indent=2))

    md = [
        "# MSR-VTT multimodal attribution\n\n",
        f"- n={point['n']}, fit_n={point['n_fit']}, layout=`{point['layout']}`\n",
        f"- RF AUC={rf_auc:.3f}; PO-risk={po['full_po_risk']:.6f}; MMD²={mmd['mmd_full']:.6f}\n\n",
        "## Modality shares\n\n",
        "| method | video | audio | text |\n|---|---:|---:|---:|\n",
        f"| RF | {rf_share['video']:.3f} | {rf_share['audio']:.3f} | {rf_share['text']:.3f} |\n",
        f"| PO-LOGO | {po['share']['video']:.3f} | {po['share']['audio']:.3f} | {po['share']['text']:.3f} |\n",
        f"| MMD-LOGO | {mmd['share']['video']:.3f} | {mmd['share']['audio']:.3f} | {mmd['share']['text']:.3f} |\n",
    ]
    (OUT / "README.md").write_text("".join(md))
    print("wrote", OUT, flush=True)
    print("SI rejects:", {k: v["rejects_equal_0"] for k, v in boot["pairwise_diff"].items()}, flush=True)


if __name__ == "__main__":
    main()
