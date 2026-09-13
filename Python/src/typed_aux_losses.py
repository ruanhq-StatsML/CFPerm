"""Auxiliary losses next to TSS: balance (conventional), typed correlation, effective rank.

Balance loss is the usual multimodal regularizer: equalize unimodal CE / gradient
energy. That target is *symmetric*. TSS is not: large covariate intensity must
not chase X, and a quiet head (clip-level text) must not be up-weighted just
because its unimodal CE is high.

Correlation with W is Cohen's d in correlation units — not a new intensity.
The redesign is the *cross-modal* Gram: penalize Corr(z^m, z^{m'}) when m is
moving and m' is quiet, so the covariate head does not drag the quiet one.
Do not run an unweighted Barlow / CCA that maximizes every pair.

Effective rank (Roy–Vetterli) is a ranking, not a third intensity next to
(ĉ, δ̂). Cov-Δerank is location-invariant, so it does not see the synthetic
mean hop. The ranking that is not already ĉ is leave-one-out Δerank of the
Amazon cosine Gram: Gift Cards and Subscription Boxes collapse erank(R).
erank of a 2×2 hop Gram is a monotone rewrite of 1−cos, not a new hop score.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from msrvtt_continuous_trainer import SeparateHeadProbe, _split_modalities
from msrvtt_multimodal_attribution import GROUP_NAMES, GROUPS, SEED, write_json
from typed_shift_stepsize import _ce, covariate_intensity, make_typed_stream, modality_mean

TAU = 1.0


def effective_rank(vals, eps=1e-12):
    """Roy–Vetterli erank: exp(entropy of the normalized spectrum)."""
    v = np.asarray(vals, dtype=float).reshape(-1)
    v = np.clip(v, 0.0, None)
    s = float(v.sum())
    if s <= eps or v.size == 0:
        return 0.0
    p = v / s
    p = p[p > eps]
    if p.size == 0:
        return 0.0
    return float(np.exp(-np.sum(p * np.log(p))))


def spectrum_psd(A):
    """Eigenvalues of a symmetric matrix, clipped to ≥0."""
    A = np.asarray(A, dtype=float)
    if A.size == 0:
        return np.zeros(0)
    if A.ndim == 1:
        return np.clip(A, 0.0, None)
    A = 0.5 * (A + A.T)
    w = np.linalg.eigvalsh(A)
    return np.clip(w, 0.0, None)


def covariance_erank(X, max_cols=64):
    """erank of the column-covariance of X (random column subset if wide)."""
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[0] < 2 or X.shape[1] < 1:
        return 0.0
    rng = np.random.default_rng(X.shape[0] * 17 + X.shape[1])
    d = X.shape[1]
    if d > int(max_cols):
        cols = rng.choice(d, size=int(max_cols), replace=False)
        X = X[:, cols]
    Xc = X - X.mean(axis=0)
    gram = (Xc.T @ Xc) / max(X.shape[0] - 1, 1)
    return effective_rank(spectrum_psd(gram))


def heatmap_erank(R):
    """erank of a cosine(Bi, Bj) Gram (already a correlation-like kernel)."""
    R = np.asarray(R, dtype=float)
    return effective_rank(spectrum_psd(R))


def _corr(a, b):
    a = np.asarray(a, dtype=float).reshape(-1)
    b = np.asarray(b, dtype=float).reshape(-1)
    if a.size < 3 or b.size != a.size:
        return 0.0
    a = a - a.mean()
    b = b - b.mean()
    da = float(np.sqrt(np.dot(a, a)))
    db = float(np.sqrt(np.dot(b, b)))
    if da < 1e-12 or db < 1e-12:
        return 0.0
    return float(np.clip(np.dot(a, b) / (da * db), -1.0, 1.0))


def unimodal_ce(probe, Xs, y):
    return {g: _ce(Xs[g] @ probe.W[g], y) for g in GROUP_NAMES}


def balance_loss(uni_ce, kind="var"):
    """Conventional modality-balance scalar.

    ``var``: sum_m (L_m - mean)^2.
    ``inv``: softmax(L_m) weights (up-weight the lagging head).
    """
    vals = np.array([float(uni_ce[g]) for g in GROUP_NAMES], dtype=float)
    mean = float(vals.mean())
    gap = {g: float(uni_ce[g] - mean) for g in GROUP_NAMES}
    if kind == "inv":
        z = vals / float(TAU)
        z = z - z.max()
        w = np.exp(z)
        w = w / (w.sum() + 1e-12)
        weights = {g: float(w[i]) for i, g in enumerate(GROUP_NAMES)}
        value = float(np.dot(w, vals))
    else:
        value = float(np.mean((vals - mean) ** 2))
        weights = {g: 1.0 / 3.0 for g in GROUP_NAMES}
    return {"value": value, "mean": mean, "gap": gap, "weights": weights, "kind": kind}


def gradient_energy(probe, Xs, y):
    """OGM-style per-head gradient energy ‖X^{m⊤} (P−e_y)‖_F."""
    from msrvtt_continuous_trainer import _softmax

    y = np.asarray(y, dtype=int)
    P = _softmax(probe.logits(Xs))
    n = max(len(y), 1)
    G = P.copy()
    G[np.arange(len(y)), y] -= 1.0
    G /= n
    energy = {}
    for g in GROUP_NAMES:
        grad = Xs[g].T @ G
        energy[g] = float(np.linalg.norm(grad))
    tot = sum(energy.values()) + 1e-12
    share = {g: energy[g] / tot for g in GROUP_NAMES}
    return {"energy": energy, "share": share}


def corr_with_batch(X0, X1):
    """Point-biserial Corr(z^m, W). Same object as ĉ, in correlation units."""
    out = {}
    for name, sl in GROUPS.items():
        z0 = modality_mean(X0[:, sl])
        z1 = modality_mean(X1[:, sl])
        z = np.concatenate([z0, z1])
        w = np.concatenate([np.zeros(len(z0)), np.ones(len(z1))])
        out[name] = _corr(z, w)
    return out


def cross_modal_corr(X):
    """3×3 Corr(z^m, z^{m'}) on one batch (coordinate-mean scalars)."""
    z = {name: modality_mean(X[:, sl]) for name, sl in GROUPS.items()}
    C = np.eye(len(GROUP_NAMES))
    for i, a in enumerate(GROUP_NAMES):
        for j, b in enumerate(GROUP_NAMES):
            C[i, j] = _corr(z[a], z[b])
    return C


def typed_correlation_loss(X0, X1, c, tau_c=0.25):
    """Cross-modal leakage weighted by covariate intensity.

    L = sum_{m ≠ m'} c_m * (1 - c_{m'}/(sum c)) * C_{mm'}^2
    on the pooled windows. A moving head is not allowed to drag a quiet head.
    Pairwise |Corr(z^m, W)| is reported, not added — it duplicates ĉ.
    """
    X = np.vstack([X0, X1])
    C = cross_modal_corr(X)
    c = {g: float(c.get(g, 0.0)) for g in GROUP_NAMES}
    tot = sum(c.values()) + 1e-12
    rho = corr_with_batch(X0, X1)
    leak = 0.0
    terms = {}
    for i, m in enumerate(GROUP_NAMES):
        for j, mp in enumerate(GROUP_NAMES):
            if i == j:
                continue
            quiet = 1.0 - c[mp] / tot
            w = c[m] * quiet if c[m] >= tau_c else 0.0
            term = float(w * (C[i, j] ** 2))
            terms["%s>%s" % (m, mp)] = term
            leak += term
    return {
        "value": float(leak),
        "C": C,
        "rho_W": rho,
        "terms": terms,
        "c": c,
    }


def erank_by_modality(X, max_cols=48):
    return {name: covariance_erank(X[:, sl], max_cols=max_cols) for name, sl in GROUPS.items()}


def two_by_two_gram_erank(cosine):
    """Roy–Vetterli erank of [[1, ρ], [ρ, 1]]. Decreasing in |ρ|; ≡ 1−cos ranking."""
    c = float(np.clip(cosine, -1.0, 1.0))
    return heatmap_erank(np.array([[1.0, c], [c, 1.0]]))


def pair_mean_cosine(X0, X1):
    """Cosine of batch-mean vectors and the heatmap hop 1−cos, per modality."""
    out = {}
    for name, sl in GROUPS.items():
        m0 = np.asarray(X0[:, sl], dtype=float).mean(axis=0)
        m1 = np.asarray(X1[:, sl], dtype=float).mean(axis=0)
        n0 = float(np.linalg.norm(m0) + 1e-12)
        n1 = float(np.linalg.norm(m1) + 1e-12)
        c = float(np.clip(np.dot(m0, m1) / (n0 * n1), -1.0, 1.0))
        out[name] = {"cosine": c, "c_heat": float(max(0.0, 1.0 - c)), "erank": two_by_two_gram_erank(c)}
    return out


def pair_mean_gram_erank(X0, X1):
    """erank of the 2×2 cosine Gram of batch means.

    A pure location shift does not change Cov(X). Consecutive means that stay
    collinear (the synthetic all-ones drift) keep cosine ≈ 1 and erank ≈ 1, so
    this is *not* a ranking of ĉ on the oracle DGP. It is the same geometry as
    an Amazon heatmap hop, where Gift vs Music is an angular split.
    """
    rec = pair_mean_cosine(X0, X1)
    return {name: rec[name]["erank"] for name in GROUP_NAMES}


def erank_delta(X0, X1, max_cols=48):
    """Within-batch Cov erank and its hop difference.

    Location-invariant: X ↦ X + a 1^⊤ does not change centered covariance, so
    a covariate mean hop can leave Δerank at sampling noise (and can rank a
    quiet head first). Do not use this as ĉ.
    """
    r0 = erank_by_modality(X0, max_cols=max_cols)
    r1 = erank_by_modality(X1, max_cols=max_cols)
    return r0, r1, {g: float(r1[g] - r0[g]) for g in GROUP_NAMES}


def rank_modalities(scores, reverse=True):
    """Stable descending rank by default (larger score first)."""
    items = [(g, float(scores[g])) for g in GROUP_NAMES]
    items.sort(key=lambda kv: ((-kv[1] if reverse else kv[1]), kv[0]))
    return [g for g, _ in items]


def pair_report(X0, X1, y0, y1, probe=None, seed=SEED):
    """One consecutive pair: balance, typed corr, erank ranking, ĉ."""
    c, share, _ = covariate_intensity(X0, X1)
    Xs1 = _split_modalities(X1)
    if probe is None:
        n_classes = int(max(int(y0.max() if len(y0) else 0), int(y1.max() if len(y1) else 0))) + 1
        probe = SeparateHeadProbe(n_classes=n_classes, seed=seed)
        warm = {g: 0.03 for g in GROUP_NAMES}
        Xs0 = _split_modalities(X0)
        for g in GROUP_NAMES:
            probe.step(Xs0, y0, warm, active=g)
    uni = unimodal_ce(probe, Xs1, y1)
    bal = balance_loss(uni, kind="var")
    inv = balance_loss(uni, kind="inv")
    ge = gradient_energy(probe, Xs1, y1)
    corr = typed_correlation_loss(X0, X1, c)
    r0, r1, dr = erank_delta(X0, X1)
    mean_hop = pair_mean_cosine(X0, X1)
    gram = {g: mean_hop[g]["erank"] for g in GROUP_NAMES}
    hop_c = {g: mean_hop[g]["c_heat"] for g in GROUP_NAMES}
    return {
        "c": c,
        "pi": share,
        "uni_ce": uni,
        "balance_var": bal,
        "balance_inv": inv,
        "grad_energy": ge,
        "typed_corr": corr,
        "erank_prev": r0,
        "erank_curr": r1,
        "erank_delta": dr,
        "erank_mean_gram": gram,
        "mean_hop": mean_hop,
        "rank_c": rank_modalities(c),
        "rank_balance_gap": rank_modalities(bal["gap"]),
        "rank_inv_weight": rank_modalities(inv["weights"]),
        "rank_erank_delta": rank_modalities(dr),
        "rank_mean_gram": rank_modalities(gram),
        "rank_mean_hop": rank_modalities(hop_c),
        "rank_rho_W": rank_modalities({g: abs(corr["rho_W"][g]) for g in GROUP_NAMES}),
    }


def report_typed_streams(n_per=48, n_batches=8, seed=SEED):
    """cov-only vs concept-only: does conventional balance invert the quiet head?"""
    out = {}
    specs = {
        "cov_only": dict(cov={"video": 0.12, "audio": 0.04, "text": 0.0}, concept={}),
        "concept_only": dict(
            cov={},
            concept={"video": 1.0},
            concept_at=4,
        ),
    }
    for name, spec in specs.items():
        stream = make_typed_stream(n_batches=n_batches, n_per=n_per, seed=seed, **spec)
        t = int(spec.get("concept_at", n_batches // 2))
        t = min(max(t, 1), int(stream.batch.max()))
        i0 = np.flatnonzero(stream.batch == t - 1)
        i1 = np.flatnonzero(stream.batch == t)
        out[name] = pair_report(stream.X[i0], stream.X[i1], stream.y[i0], stream.y[i1], seed=seed)
        out[name]["pair"] = (t - 1, t)
    return out


def amazon_heatmap_ranks(R, labels=None, hops=None):
    """Rank hops on a category Gram: 1−cos vs leave-one-out Δerank."""
    R = np.asarray(R, dtype=float)
    k = R.shape[0]
    labels = list(labels or ["B%d" % i for i in range(k)])
    base = heatmap_erank(R)
    hop_c = []
    for t in range(1, k):
        hop_c.append(
            {
                "hop": t,
                "from": labels[t - 1],
                "to": labels[t],
                "cosine": float(R[t - 1, t]),
                "c_heat": float(hops[t - 1]) if hops is not None else float(max(0.0, 1.0 - R[t - 1, t])),
                "gram_erank": two_by_two_gram_erank(R[t - 1, t]),
            }
        )
    loo = []
    for i in range(k):
        mask = np.ones(k, dtype=bool)
        mask[i] = False
        Rsub = R[np.ix_(mask, mask)]
        loo.append(
            {
                "i": i,
                "label": labels[i],
                "erank_loo": heatmap_erank(Rsub),
                "delta_erank": float(base - heatmap_erank(Rsub)),
            }
        )
    loo.sort(key=lambda r: -r["delta_erank"])
    hop_c.sort(key=lambda r: -r["c_heat"])
    return {
        "erank": base,
        "k": k,
        "hops_by_c": hop_c,
        "categories_by_loo_erank": loo,
    }


def write_tex(report, amazon, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cov = report["cov_only"]
    con = report["concept_only"]

    def row(rec, key, fmt="%.3f"):
        if key == "uni":
            return " / ".join(fmt % rec["uni_ce"][g] for g in GROUP_NAMES)
        if key == "gap":
            return " / ".join(fmt % rec["balance_var"]["gap"][g] for g in GROUP_NAMES)
        if key == "inv":
            return " / ".join(fmt % rec["balance_inv"]["weights"][g] for g in GROUP_NAMES)
        if key == "c":
            return " / ".join(fmt % rec["c"][g] for g in GROUP_NAMES)
        if key == "er":
            return " / ".join(fmt % rec["erank_delta"][g] for g in GROUP_NAMES)
        if key == "gram":
            return " / ".join(fmt % rec["erank_mean_gram"][g] for g in GROUP_NAMES)
        if key == "rho":
            return " / ".join(fmt % abs(rec["typed_corr"]["rho_W"][g]) for g in GROUP_NAMES)
        return ""

    lines = [
        r"% Conventional balance vs typed correlation vs effective-rank. Auto-generated.",
        r"\begin{table}[ht]\centering",
        r"\caption{Unimodal CE balance is conventional and symmetric: the lagging head gets the largest inverse-CE weight.",
        r"On a covariate-only hop that is the wrong sign for clip-level text.",
        r"$\mathrm{Corr}(z^{(m)},W)$ duplicates $\hat c_m$.",
        r"Typed correlation penalizes $C_{mm'}^2$ only when $m$ is moving and $m'$ is quiet.",
        r"$\Delta\mathrm{erank}(\mathrm{Cov}\,X^{(m)})$ is location-invariant and does not rank the mean hop.",
        r"Amazon leave-one-out $\Delta\mathrm{erank}(R)$ ranks domains; a $2\times 2$ hop Gram is a monotone rewrite of $1-\cos$.}",
        r"\label{tab:typed-aux-losses}",
        r"\small",
        r"\begin{tabular}{@{}ll cc@{}}\toprule",
        r"& & covariate-only hop & concept hop \\",
        r"\midrule",
        r"$\hat c_{v/a/t}$ & & %s & %s \\" % (row(cov, "c"), row(con, "c")),
        r"unimodal CE & & %s & %s \\" % (row(cov, "uni"), row(con, "uni")),
        r"balance gap $L_m-\bar L$ & & %s & %s \\" % (row(cov, "gap"), row(con, "gap")),
        r"inv-CE weights & & %s & %s \\" % (row(cov, "inv"), row(con, "inv")),
        r"$|\mathrm{Corr}(z^{(m)},W)|$ & & %s & %s \\" % (row(cov, "rho"), row(con, "rho")),
        r"typed corr $L$ & & $%.3f$ & $%.3f$ \\" % (cov["typed_corr"]["value"], con["typed_corr"]["value"]),
        r"$\Delta\mathrm{erank}(\mathrm{Cov})_{v/a/t}$ & & %s & %s \\" % (row(cov, "er"), row(con, "er")),
        r"erank $2\times2$ mean-Gram & & %s & %s \\" % (row(cov, "gram"), row(con, "gram")),
        r"rank by $\hat c$ & & %s & %s \\" % (" $>$ ".join(cov["rank_c"]), " $>$ ".join(con["rank_c"])),
        r"rank by inv-CE & & %s & %s \\" % (" $>$ ".join(cov["rank_inv_weight"]), " $>$ ".join(con["rank_inv_weight"])),
        r"rank by $\Delta\mathrm{erank}(\mathrm{Cov})$ & & %s & %s \\"
        % (" $>$ ".join(cov["rank_erank_delta"]), " $>$ ".join(con["rank_erank_delta"])),
        r"rank by $2\times 2$ mean-Gram erank & & %s & %s \\"
        % (" $>$ ".join(cov["rank_mean_gram"]), " $>$ ".join(con["rank_mean_gram"])),
        r"\bottomrule",
        r"\end{tabular}\\[0.4em]",
    ]
    if amazon:
        top = amazon["categories_by_loo_erank"][:3]
        hops = amazon["hops_by_c"][:3]
        lines.append(
            r"{\footnotesize Amazon heatmap $\mathrm{erank}(R)=%.2f$ of $K=%d$."
            r" Leave-one-out $\Delta\mathrm{erank}$ leaders: %s."
            r" Largest hops $1-\cos$: %s.}"
            % (
                amazon["erank"],
                amazon["k"],
                ", ".join("%s ($%.3f$)" % (r["label"].replace("_", r"\_"), r["delta_erank"]) for r in top),
                ", ".join("%s$\\to$%s ($%.3f$)" % (h["from"].replace("_", r"\_"), h["to"].replace("_", r"\_"), h["c_heat"]) for h in hops),
            )
        )
    lines.append(r"\end{table}")
    path.write_text("\n".join(lines) + "\n")
    return path


def plot_aux_comparison(report, path, amazon=None):
    import matplotlib.pyplot as plt
    from msrvtt_attribution_plots import AUDIO_C, GRID, INK, MUTED, TEXT_C, VIDEO_C, _save, _style

    _style()
    colors = {"video": VIDEO_C, "audio": AUDIO_C, "text": TEXT_C}
    cov, con = report["cov_only"], report["concept_only"]
    n_rows = 2 if amazon else 1
    fig, axes = plt.subplots(n_rows, 3, figsize=(12.2, 3.7 * n_rows))
    if n_rows == 1:
        axes = np.asarray(axes).reshape(1, 3)
    x = np.arange(len(GROUP_NAMES))
    width = 0.36

    def grouped(ax, a, b, title, ylabel=""):
        ax.bar(x - width / 2, a, width=width, color=[colors[g] for g in GROUP_NAMES], alpha=0.95, label="covariate hop")
        ax.bar(x + width / 2, b, width=width, color=[colors[g] for g in GROUP_NAMES], alpha=0.40, label="concept hop")
        ax.set_xticks(x)
        ax.set_xticklabels(GROUP_NAMES)
        ax.set_title(title, loc="left", fontsize=10.2, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", color=GRID)
        ax.axhline(0, color=MUTED, lw=0.7)

    grouped(
        axes[0, 0],
        [cov["balance_inv"]["weights"][g] for g in GROUP_NAMES],
        [con["balance_inv"]["weights"][g] for g in GROUP_NAMES],
        "inv-CE weight (conventional balance)",
    )
    grouped(
        axes[0, 1],
        [abs(cov["typed_corr"]["rho_W"][g]) for g in GROUP_NAMES],
        [abs(con["typed_corr"]["rho_W"][g]) for g in GROUP_NAMES],
        r"$|\mathrm{Corr}(z^{(m)},W)|$  $\equiv\hat c$",
    )
    term_keys = ["video>audio", "video>text", "audio>text"]
    axes[0, 2].bar(
        np.arange(len(term_keys)) - width / 2,
        [cov["typed_corr"]["terms"][k] for k in term_keys],
        width=width,
        color=VIDEO_C,
        alpha=0.95,
        label="covariate hop",
    )
    axes[0, 2].bar(
        np.arange(len(term_keys)) + width / 2,
        [con["typed_corr"]["terms"][k] for k in term_keys],
        width=width,
        color=VIDEO_C,
        alpha=0.40,
        label="concept hop",
    )
    axes[0, 2].set_xticks(np.arange(len(term_keys)))
    axes[0, 2].set_xticklabels(term_keys, fontsize=8.5)
    axes[0, 2].set_title(r"typed $L_{\mathrm{corr}}$ terms (moving$>$quiet)", loc="left", fontsize=10.2, fontweight="bold")
    axes[0, 2].grid(True, axis="y", color=GRID)
    axes[0, 0].legend(frameon=False, fontsize=7.5)

    if amazon:
        loo = amazon["categories_by_loo_erank"]
        labels = [r["label"].replace("_", " ") for r in loo]
        axes[1, 0].barh(np.arange(len(loo))[::-1], [r["delta_erank"] for r in loo][::-1], color="#2C4A6E")
        axes[1, 0].set_yticks(np.arange(len(loo))[::-1])
        axes[1, 0].set_yticklabels(labels[::-1], fontsize=8)
        axes[1, 0].set_title(r"Amazon leave-one-out $\Delta\mathrm{erank}(R)$", loc="left", fontsize=10.2, fontweight="bold")
        axes[1, 0].grid(True, axis="x", color=GRID)
        hops = amazon["hops_by_c"]
        hlab = ["%s→%s" % (h["from"][:6], h["to"][:6]) for h in hops]
        axes[1, 1].barh(np.arange(len(hops))[::-1], [h["c_heat"] for h in hops][::-1], color="#C45C26")
        axes[1, 1].set_yticks(np.arange(len(hops))[::-1])
        axes[1, 1].set_yticklabels(hlab[::-1], fontsize=8)
        axes[1, 1].set_title(r"Amazon hops $1-\cos(\bar x_{t-1},\bar x_t)$", loc="left", fontsize=10.2, fontweight="bold")
        axes[1, 1].grid(True, axis="x", color=GRID)
        axes[1, 2].scatter(
            [h["c_heat"] for h in hops],
            [h["gram_erank"] for h in hops],
            c="#2F6B4F",
            s=42,
            zorder=3,
        )
        for h in hops[:3]:
            axes[1, 2].annotate(
                "%s→%s" % (h["from"][:4], h["to"][:4]),
                (h["c_heat"], h["gram_erank"]),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=7.5,
                color=MUTED,
            )
        axes[1, 2].set_xlabel(r"$1-\cos$")
        axes[1, 2].set_ylabel(r"erank $2\times 2$ Gram")
        axes[1, 2].set_title("hop Gram erank ≡ $1-\\cos$ ranking", loc="left", fontsize=10.2, fontweight="bold")
        axes[1, 2].grid(True, color=GRID)

    fig.suptitle(
        r"Balance is conventional; Corr$(z,W)$ is $\hat c$; typed $L_{\mathrm{corr}}$ is the off-diagonal; erank ranks the heatmap",
        fontsize=12.0,
        fontweight="bold",
        color=INK,
        x=0.04,
        ha="left",
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.91))
    return _save(fig, path)
