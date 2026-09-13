#!/usr/bin/env python3
"""AGOD smoke test on MSR-VTT packed features.

Attribution-Guided Online Distillation:
  Monitor MSG (RF-AUC × VIMP + γ·PO) → Softmax route α → weight modality distill losses.

Baselines:
  B1 static α=1/3
  B2 covariate-only α ∝ AUC
  B3 AGOD α ∝ AUC·VIMP + γ·PO

Proxy distill loss (no full student train): per-modality ridge reconstruction of
teacher features on the drifted batch; report weighted loss under each α policy.

  python3 scripts/run_agod_msrvtt_smoke.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "msrvtt" / "packed"
OUT = ROOT / "results" / "agod"
ART = Path("/opt/cursor/artifacts/msrvtt_dashboard")
DOCS = ROOT / "docs" / "agod"

SEED = 2026
T_STEPS = 8
GAMMA = 1.0
TAU = 0.35
N_PER = 350  # samples per window side
RF_TREES = 80
PO_TREES = 40
ALPHA_CLIP = 0.05

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
    X = np.hstack([V, A, T])
    pc1 = PCA(n_components=1, random_state=SEED).fit_transform(X).ravel()
    order = np.argsort(pc1, kind="mergesort")
    ranks = np.empty(len(pc1), float)
    ranks[order] = np.linspace(0.0, 1.0, len(pc1))
    return X, ranks, W


def make_windows(X, Y, W, *, t_steps=T_STEPS, seed=SEED):
    """Reference = W==0 pool; online steps = stratified chunks of W==1 (and mixed)."""
    rng = np.random.default_rng(seed)
    i0, i1 = np.where(W == 0)[0], np.where(W == 1)[0]
    ref = rng.choice(i0, min(N_PER, len(i0)), False)
    # stream late + some early contamination as evolving batches
    perm1 = rng.permutation(i1)
    chunks = np.array_split(perm1, t_steps)
    windows = []
    for t, ch in enumerate(chunks):
        # mix a shrinking fraction of reference to simulate gradual shift
        n_ref = max(20, int(0.25 * (1 - t / max(t_steps - 1, 1)) * len(ch)))
        extra = rng.choice(i0, min(n_ref, len(i0)), False)
        idx = np.concatenate([ch, extra])
        rng.shuffle(idx)
        windows.append(idx)
    return ref, windows


def rf_auc_vimp(X0, X1, *, seed):
    X = np.vstack([X0, X1])
    W = np.array([0] * len(X0) + [1] * len(X1))
    if len(np.unique(W)) < 2 or len(X) < 40:
        return 0.5, np.zeros(X.shape[1])
    clf = RandomForestClassifier(
        n_estimators=RF_TREES,
        max_depth=max(3, int(np.sqrt(X.shape[1]))),
        min_samples_leaf=3,
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    Xtr, Xte, Wtr, Wte = train_test_split(
        X, W, test_size=0.3, random_state=seed, stratify=W
    )
    clf2 = RandomForestClassifier(
        n_estimators=60,
        max_depth=max(3, int(np.sqrt(X.shape[1]))),
        min_samples_leaf=3,
        n_jobs=-1,
        random_state=seed + 1,
    )
    clf2.fit(Xtr, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return auc, vimp


def po_risk_block(X0, Y0, X1, Y1, *, seed):
    X = np.vstack([X0, X1])
    Y = np.concatenate([Y0, Y1]).astype(float)
    W = np.array([0] * len(X0) + [1] * len(X1))
    if len(np.unique(W)) < 2:
        return 0.0
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=PO_TREES, max_depth=6, min_samples_leaf=4,
            random_state=seed + fold, n_jobs=-1,
        )
        e = RandomForestClassifier(
            n_estimators=PO_TREES, max_depth=6, min_samples_leaf=4,
            random_state=seed + 20 + fold, n_jobs=-1,
        )
        m.fit(X[tr], Y[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1 - ALPHA_CLIP)
    po = (Y - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=80, max_depth=8, min_samples_leaf=4,
        random_state=seed + 7, n_jobs=-1,
    )
    tau.fit(X, po)
    return float(np.mean(tau.predict(X) ** 2))


def softmax(z, tau=TAU):
    z = np.asarray(z, float) / max(tau, 1e-6)
    z = z - z.max()
    e = np.exp(z)
    return e / e.sum()


def normalize_gap(g):
    g = np.asarray(g, float)
    g = np.maximum(g, 0.0)
    s = g.sum()
    return g / s if s > 0 else np.ones_like(g) / len(g)


def proxy_modality_loss(X_ref, X_t, *, seed):
    """Teacher=identity features; student=Ridge fit on ref, MSE on drifted batch."""
    rng = np.random.default_rng(seed)
    # subsample for speed
    n0 = min(400, len(X_ref))
    n1 = min(400, len(X_t))
    A = X_ref[rng.choice(len(X_ref), n0, False)]
    B = X_t[rng.choice(len(X_t), n1, False)]
    # map ref→ref (trivial teacher match), evaluate on B
    model = Ridge(alpha=1.0, random_state=seed)
    model.fit(A, A)
    pred = model.predict(B)
    return float(np.mean((pred - B) ** 2))


def compute_msg(X, Y, ref, idx_t, *, seed):
    aucs, vimps, pos, gaps = {}, {}, {}, {}
    for m, (a, b) in BLOCKS.items():
        auc, vimp = rf_auc_vimp(X[ref, a:b], X[idx_t, a:b], seed=seed + MODS.index(m))
        vimp_mass = float(vimp.sum())  # already block-only
        # use mean vimp as scale-free mass proxy within block
        vimp_score = float(vimp.mean()) if len(vimp) else 0.0
        por = po_risk_block(
            X[ref, a:b], Y[ref], X[idx_t, a:b], Y[idx_t], seed=seed + 10 + MODS.index(m)
        )
        aucs[m] = auc
        vimps[m] = vimp_score
        pos[m] = por
        gaps[m] = auc * vimp_score + GAMMA * por
    g_vec = normalize_gap([gaps[m] for m in MODS])
    return {
        "auc": aucs,
        "vimp": vimps,
        "po": pos,
        "gap_raw": gaps,
        "g": {m: float(g_vec[i]) for i, m in enumerate(MODS)},
    }


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    print("load pack…", flush=True)
    X, Y, W = load_pack()
    ref, windows = make_windows(X, Y, W)
    print(f"ref={len(ref)} steps={len(windows)}", flush=True)

    traj = []
    losses = {"B1": [], "B2": [], "B3": []}
    alpha_hist = {"B1": [], "B2": [], "B3": []}

    for t, idx in enumerate(windows):
        print(f"[t={t+1}/{len(windows)}] n={len(idx)}", flush=True)
        msg = compute_msg(X, Y, ref, idx, seed=SEED + t)
        auc_vec = np.array([msg["auc"][m] for m in MODS])
        g_vec = np.array([msg["g"][m] for m in MODS])

        a_b1 = np.ones(3) / 3.0
        a_b2 = softmax(np.maximum(auc_vec - 0.5, 0.0) + 1e-6, tau=TAU)
        a_b3 = softmax(g_vec + 1e-6, tau=TAU)

        # modality proxy losses on this step
        Lm = {}
        for m, (a, b) in BLOCKS.items():
            Lm[m] = proxy_modality_loss(X[ref, a:b], X[idx, a:b], seed=SEED + 100 + t)

        def weighted(alpha):
            return float(sum(alpha[i] * Lm[MODS[i]] for i in range(3)))

        row = {
            "t": t,
            "msg": msg,
            "L_m": Lm,
            "alpha": {"B1": a_b1.tolist(), "B2": a_b2.tolist(), "B3": a_b3.tolist()},
            "L": {"B1": weighted(a_b1), "B2": weighted(a_b2), "B3": weighted(a_b3)},
        }
        traj.append(row)
        for k in losses:
            losses[k].append(row["L"][k])
            alpha_hist[k].append(row["alpha"][k])
        print(
            f"  AUC={ {m: round(msg['auc'][m],3) for m in MODS} } "
            f"g={ {m: round(msg['g'][m],3) for m in MODS} } "
            f"α_AGOD={np.round(a_b3,3).tolist()} "
            f"L B1/B2/B3={row['L']['B1']:.4f}/{row['L']['B2']:.4f}/{row['L']['B3']:.4f}",
            flush=True,
        )

    summary = {
        "T": len(windows),
        "gamma": GAMMA,
        "tau": TAU,
        "mean_loss": {k: float(np.mean(v)) for k, v in losses.items()},
        "final_alpha_B3": traj[-1]["alpha"]["B3"],
        "mods": MODS,
        "trajectory": traj,
    }
    (OUT / "agod_msrvtt_smoke.json").write_text(json.dumps(summary, indent=2))

    # ---- plots ----
    colors = {"video": "#2B6CB0", "audio": "#C05621", "text": "#276749"}
    # 1) alpha trajectories AGOD
    fig, ax = plt.subplots(figsize=(8.2, 4.0), dpi=140)
    A3 = np.asarray(alpha_hist["B3"])
    for i, m in enumerate(MODS):
        ax.plot(np.arange(1, len(A3) + 1), A3[:, i], marker="o", label=m, color=colors[m], lw=2)
    ax.set_xlabel("online step t")
    ax.set_ylabel(r"$\alpha_m^{(t)}$")
    ax.set_title("AGOD routing weights (Softmax MSG)")
    ax.set_ylim(0, 1)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "agod_alpha_trajectory.png", bbox_inches="tight")
    fig.savefig(ART / "agod_alpha_trajectory.png", bbox_inches="tight")
    plt.close(fig)

    # 2) MSG components over time
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), dpi=140)
    for ax, key, title in zip(
        axes,
        ["auc", "g", "po"],
        ["RF-Domain AUC", "MSG g (normalized)", "PO-risk"],
    ):
        for m in MODS:
            ys = [traj[t]["msg"][key][m] for t in range(len(traj))]
            ax.plot(np.arange(1, len(ys) + 1), ys, marker="o", label=m, color=colors[m], lw=1.8)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("t")
    axes[0].legend(fontsize=8)
    fig.suptitle("Modality-Specific Gap monitor on MSR-VTT stream", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "agod_msg_monitor.png", bbox_inches="tight")
    fig.savefig(ART / "agod_msg_monitor.png", bbox_inches="tight")
    plt.close(fig)

    # 3) baseline loss comparison
    fig, ax = plt.subplots(figsize=(7.2, 3.8), dpi=140)
    x = np.arange(1, len(losses["B1"]) + 1)
    ax.plot(x, losses["B1"], "o-", label="B1 static", color="#718096")
    ax.plot(x, losses["B2"], "s-", label="B2 AUC-only", color="#DD6B20")
    ax.plot(x, losses["B3"], "D-", label="B3 AGOD", color="#C53030", lw=2)
    ax.set_xlabel("online step t")
    ax.set_ylabel("proxy weighted distill loss")
    ax.set_title("Attribution-guided vs static / covariate-only routing")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "agod_loss_compare.png", bbox_inches="tight")
    fig.savefig(ART / "agod_loss_compare.png", bbox_inches="tight")
    plt.close(fig)

    # 4) product dashboard collage
    fig = plt.figure(figsize=(12.5, 8.5), dpi=150)
    fig.patch.set_facecolor("#F7F5F1")
    fig.suptitle(
        "MSR-VTT Multimodal Attribution → AGOD Online Routing Board",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    # top: existing attribution share if available
    share_p = ROOT / "results" / "msrvtt_mm_attr" / "modality_share_compare.png"
    auc_p = ROOT / "results" / "msrvtt_mm_attr" / "modality_subset_auc.png"
    ax1 = fig.add_axes([0.04, 0.52, 0.58, 0.40])
    ax2 = fig.add_axes([0.64, 0.52, 0.32, 0.40])
    ax3 = fig.add_axes([0.04, 0.08, 0.45, 0.38])
    ax4 = fig.add_axes([0.54, 0.08, 0.42, 0.38])
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xticks([])
        ax.set_yticks([])
    if share_p.exists():
        ax1.imshow(plt.imread(share_p))
        ax1.set_title("Offline attribution (RF / PO-LOGO / MMD)", fontsize=9)
    if auc_p.exists():
        ax2.imshow(plt.imread(auc_p))
        ax2.set_title("Per-modality subset AUC", fontsize=9)
    ax3.imshow(plt.imread(OUT / "agod_alpha_trajectory.png"))
    ax3.set_title("Online AGOD α routing", fontsize=9)
    ax4.imshow(plt.imread(OUT / "agod_loss_compare.png"))
    ax4.set_title("B1 / B2 / B3 proxy loss", fontsize=9)
    fig.text(
        0.5,
        0.01,
        "Philosophy: train on the shift — MSG = clever covariate → Softmax route distillation capacity",
        ha="center",
        fontsize=9,
        color="#4A5568",
    )
    dash = OUT / "AGOD_MSR_VTT_Product_Dashboard.png"
    fig.savefig(dash, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(ART / "AGOD_MSR_VTT_Product_Dashboard.png", bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)

    # latex mini table
    tex = [
        "% AGOD smoke on MSR-VTT\n% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{AGOD smoke test on MSR-VTT packed stream ($T="
        + str(len(windows))
        + "$). Proxy weighted distill loss under static / AUC-only / AGOD routing.}\n",
        "\\label{tab:agod-smoke}\n\\small\n",
        "\\begin{tabular}{@{}l c@{}}\n\\toprule\n",
        "Policy & Mean proxy loss \\\\\n\\midrule\n",
        f"B1 Static $\\alpha=1/3$ & ${summary['mean_loss']['B1']:.4f}$ \\\\\n",
        f"B2 Covariate-only (AUC) & ${summary['mean_loss']['B2']:.4f}$ \\\\\n",
        f"B3 AGOD (AUC$\\cdot$VIMP+$\\gamma$ PO) & ${summary['mean_loss']['B3']:.4f}$ \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n",
    ]
    (OUT / "AGOD_smoke_tables_only.tex").write_text("".join(tex))
    (DOCS / "AGOD_smoke_tables_only.tex").write_text("".join(tex))

    md = [
        "# AGOD smoke · MSR-VTT\n\n",
        f"- mean loss B1/B2/B3 = "
        f"{summary['mean_loss']['B1']:.4f} / {summary['mean_loss']['B2']:.4f} / "
        f"{summary['mean_loss']['B3']:.4f}\n",
        f"- final α_AGOD = {summary['final_alpha_B3']}\n",
        "- dashboard: `AGOD_MSR_VTT_Product_Dashboard.png`\n",
    ]
    (OUT / "README.md").write_text("".join(md))
    print("DONE", summary["mean_loss"], flush=True)
    print("dashboard", dash, flush=True)


if __name__ == "__main__":
    run()
