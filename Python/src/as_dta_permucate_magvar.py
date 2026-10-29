"""Mag×Var on PermuCATE-style PO-risk (RF outcome models) + multimodal attribution.

Estimator (matches ``permuCATE_vimp.py`` formulas; sklearn RF stand-in)
-----------------------------------------------------------------------
Outcome:  μ̂_0, μ̂_1 = RandomForestRegressor on W=0 / W=1
Propensity: π̂ = RandomForestClassifier (or logistic) for P(W=1|X)
DR pseudo-outcome:
  φ = (μ̂_1-μ̂_0) + (Y-μ̂_W)(W-π̂)/(π̂(1-π̂))
CATE: τ̂ = RF regress φ ~ X
Unit PO-risk:  u_i = (φ_i - τ̂_i)^2
Magnitude μ = E[u], Variance s = Std(r_b) over bootstrap sub-envs.

Honest note
-----------
Earlier AIGC Mag×Var boards used supervised proxy PO=1-p̂(Y|x).
This module is the **causal** MVP: same Mag×Var geometry on PermuCATE φ-risk.

Multimodal attribution
----------------------
Group-permute one modality block (V/A/T), remeasure mean unit risk;
share_m ∝ max(Δu_m, 0). Mag×Var decides *shift type*; modality shares
decide *which channel* to freeze.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import StratifiedKFold

from as_dta_aigc import ROOT, SEED
from fsds_rl_closed_loop import make_mm_hop_proxy

OUT_DEFAULT = ROOT / "results" / "as_dta_permucate_magvar"
CLIP = 1e-2


# ---------------------------------------------------------------------------
# PermuCATE-style RF nuisances + unit PO-risk
# ---------------------------------------------------------------------------


def _rf_reg(n: int, seed: int, n_estimators: int = 80) -> RandomForestRegressor:
    return RandomForestRegressor(
        n_estimators=int(n_estimators),
        max_depth=5 if n >= 80 else 3,
        min_samples_leaf=max(2, n // 50),
        random_state=int(seed),
        n_jobs=1,
    )


def _rf_clf(n: int, seed: int, n_estimators: int = 80) -> RandomForestClassifier:
    return RandomForestClassifier(
        n_estimators=int(n_estimators),
        max_depth=5 if n >= 80 else 3,
        min_samples_leaf=max(2, n // 50),
        random_state=int(seed),
        n_jobs=1,
    )


def permucate_unit_po(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    *,
    seed: int = SEED,
    n_estimators: int = 80,
    n_splits: int = 3,
    clip: float = CLIP,
) -> dict:
    """Cross-fit RF outcome + propensity → φ, τ̂, unit risk u=(φ-τ̂)^2.

    Mirrors ``po_risk`` / ``permuCATE_fit_nuisance`` in permuCATE_vimp.py
    with RandomForest as ``model_m`` / ``model_tau``.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).reshape(-1)
    W = np.asarray(W, dtype=int).reshape(-1)
    n, p = X.shape
    mu0 = np.zeros(n)
    mu1 = np.zeros(n)
    pi = np.zeros(n)

    # Cross-fit outcome models by W and propensity
    splits = max(2, min(int(n_splits), int(np.min(np.bincount(W)))))
    cv = StratifiedKFold(n_splits=splits, shuffle=True, random_state=int(seed))
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        tr0 = tr[W[tr] == 0]
        tr1 = tr[W[tr] == 1]
        m0 = _rf_reg(len(tr0), seed + fold, n_estimators)
        m1 = _rf_reg(len(tr1), seed + 17 + fold, n_estimators)
        if len(tr0) >= 8:
            m0.fit(X[tr0], Y[tr0])
            mu0[te] = m0.predict(X[te])
        else:
            mu0[te] = float(np.mean(Y[W == 0])) if (W == 0).any() else 0.0
        if len(tr1) >= 8:
            m1.fit(X[tr1], Y[tr1])
            mu1[te] = m1.predict(X[te])
        else:
            mu1[te] = float(np.mean(Y[W == 1])) if (W == 1).any() else 0.0
        e = _rf_clf(len(tr), seed + 40 + fold, n_estimators)
        e.fit(X[tr], W[tr])
        pi[te] = e.predict_proba(X[te])[:, 1]

    pi = np.clip(pi, clip, 1.0 - clip)
    mu_w = np.where(W == 1, mu1, mu0)
    phi = (mu1 - mu0) + (Y - mu_w) * (W - pi) / (pi * (1.0 - pi))

    tau_m = _rf_reg(n, seed + 7, n_estimators)
    tau_m.fit(X, phi)
    tau = tau_m.predict(X)
    unit = (phi - tau) ** 2
    vimp = np.maximum(tau_m.feature_importances_.astype(float), 0.0)
    if float(vimp.sum()) <= 0:
        vimp = np.ones(p) / max(p, 1)
    else:
        vimp = vimp / vimp.sum()

    return {
        "phi": phi,
        "tau": tau,
        "unit_po": unit,
        "mu0": mu0,
        "mu1": mu1,
        "pi": pi,
        "r_po": float(np.mean(unit)),
        "ate_proxy": float(np.mean(phi)),
        "vimp": vimp,
        "estimator": "permucate_rf_outcome_crossfit",
        "model_m": "RandomForestRegressor",
        "model_e": "RandomForestClassifier",
        "model_tau": "RandomForestRegressor",
    }


def mag_var_from_unit(unit: np.ndarray, K: int = 36, seed: int = SEED) -> dict:
    """Magnitude μ + sub-env variance s of unit PO-risk."""
    unit = np.asarray(unit, dtype=float)
    n = len(unit)
    rng = np.random.default_rng(int(seed))
    win = max(12, n // 5)
    means = []
    for _ in range(int(K)):
        idx = rng.choice(n, size=min(win, n), replace=False)
        means.append(float(np.mean(unit[idx])))
    means = np.asarray(means, dtype=float)
    return {
        "mu": float(np.mean(unit)),
        "s": float(np.std(means)),
        "v": float(np.var(means)),
        "unit_var": float(np.var(unit)),
        "raw_sub_means": means.tolist(),
    }


def merge_J(d_mu: float, d_s: float, s_ref: float, kappa: float = 3.0) -> float:
    s_ref = max(float(s_ref), 1e-9)
    return float(d_mu * (1.0 + kappa * max(d_s, 0.0) / s_ref))


def pred_shift(d_mu: float, d_s: float, J: float = 0.0, d_err: float = 0.0) -> str:
    """Mag×Var gate on PermuCATE unit risk (μ, s, merge J)."""
    if d_mu < 0.02 and abs(d_s) < 0.01 and J < 0.05:
        return "quiet"
    # Malignant: large magnitude and/or variance/J amplification
    if (
        J >= 0.65
        or d_mu >= 0.55
        or d_s >= 0.05
        or (d_mu >= 0.35 and d_s >= 0.0)
        or d_err >= 0.05
    ):
        return "malignant"
    if d_mu >= 0.05:
        return "benign"
    return "quiet"


# ---------------------------------------------------------------------------
# Tabular boards: benign uniform shift vs malignant subset shift
# ---------------------------------------------------------------------------


def _make_base(n: int = 400, p: int = 24, seed: int = SEED) -> dict:
    rng = np.random.default_rng(int(seed))
    X = rng.normal(size=(n, p))
    # latent effect on a few coords
    beta = np.zeros(p)
    beta[:4] = rng.uniform(0.6, 1.2, size=4)
    Y = X @ beta + rng.normal(scale=0.5, size=n)
    return {"X": X, "Y": Y, "beta": beta, "rng": rng}


def build_shift_pair(
    n: int = 360,
    p: int = 24,
    seed: int = SEED,
) -> dict:
    """REF (W=0) vs LIVE (W=1) under benign or malignant covariate shift."""
    base = _make_base(n=n, p=p, seed=seed)
    X0, Y0 = base["X"], base["Y"]
    rng = np.random.default_rng(int(seed) + 3)

    # Benign: mild uniform mean shift on all features
    X_ben = X0 + rng.normal(scale=0.35, size=X0.shape)
    Y_ben = X_ben @ base["beta"] + rng.normal(scale=0.5, size=n)

    # Malignant: hard shift on a *subset* of rows + attack coords (mixture)
    X_mal = X0.copy()
    S = rng.choice(n, size=max(48, n // 3), replace=False)
    attack_cols = np.arange(0, min(8, p))
    X_mal[np.ix_(S, attack_cols)] += rng.normal(loc=3.5, scale=0.5, size=(len(S), len(attack_cols)))
    Y_mal = X_mal @ base["beta"] + rng.normal(scale=0.5, size=n)
    # outcome flip / reward break on attacked subset only (heteroscedastic PO)
    Y_mal[S] = -Y_mal[S] + rng.normal(scale=1.2, size=len(S))

    return {
        "ref": {"X": X0, "Y": Y0},
        "benign": {"X": X_ben, "Y": Y_ben},
        "malignant": {"X": X_mal, "Y": Y_mal, "S": S.tolist()},
        "p": p,
        "n": n,
    }


def eval_magvar_permucate(
    n: int = 360,
    p: int = 24,
    K: int = 36,
    seed: int = SEED,
    n_estimators: int = 80,
) -> dict:
    pair = build_shift_pair(n=n, p=p, seed=seed)
    rows = {}
    ref_pack = None
    for name in ("ref", "benign", "malignant"):
        # Compare REF vs LIVE: stack W=0 ref, W=1 live (for ref self, split half)
        if name == "ref":
            X = pair["ref"]["X"]
            Y = pair["ref"]["Y"]
            mid = len(X) // 2
            X_stack = np.vstack([X[:mid], X[mid:]])
            Y_stack = np.concatenate([Y[:mid], Y[mid:]])
            W = np.concatenate([np.zeros(mid, dtype=int), np.ones(len(X) - mid, dtype=int)])
        else:
            X0, Y0 = pair["ref"]["X"], pair["ref"]["Y"]
            X1, Y1 = pair[name]["X"], pair[name]["Y"]
            X_stack = np.vstack([X0, X1])
            Y_stack = np.concatenate([Y0, Y1])
            W = np.concatenate([np.zeros(len(X0), dtype=int), np.ones(len(X1), dtype=int)])

        fit = permucate_unit_po(
            X_stack, Y_stack, W, seed=seed + abs(hash(name)) % 1000, n_estimators=n_estimators
        )
        # Mag×Var on LIVE units only (W=1)
        unit_live = fit["unit_po"][W == 1]
        mv = mag_var_from_unit(unit_live, K=K, seed=seed + 9)
        rows[name] = {"fit": fit, "magvar": mv, "n_live": int((W == 1).sum())}
        if name == "ref":
            ref_pack = mv

    assert ref_pack is not None
    contingency = {"tn": 0, "fp": 0, "fn": 0, "tp": 0}
    classified = {}
    for truth in ("benign", "malignant"):
        d_mu = rows[truth]["magvar"]["mu"] - ref_pack["mu"]
        d_s = rows[truth]["magvar"]["s"] - ref_pack["s"]
        J = merge_J(d_mu, d_s, ref_pack["s"])
        pred = pred_shift(d_mu, d_s, J=J)
        if pred == "quiet":
            pred = "benign"
        classified[truth] = {
            "true": truth,
            "pred": pred,
            "d_mu": d_mu,
            "d_s": d_s,
            "J": J,
            "mu": rows[truth]["magvar"]["mu"],
            "s": rows[truth]["magvar"]["s"],
            "r_po": rows[truth]["fit"]["r_po"],
        }
        if truth == "benign" and pred == "benign":
            contingency["tn"] += 1
        elif truth == "benign" and pred == "malignant":
            contingency["fp"] += 1
        elif truth == "malignant" and pred == "benign":
            contingency["fn"] += 1
        else:
            contingency["tp"] += 1

    return {
        "estimator": {
            "name": "PermuCATE-style RF",
            "model_m": "RandomForestRegressor (μ0, μ1)",
            "model_e": "RandomForestClassifier (π)",
            "model_tau": "RandomForestRegressor (τ̂ on φ)",
            "unit_po": "(φ - τ̂)^2",
            "used_permucate_rf": True,
            "note": "Same DR φ as permuCATE_vimp.po_risk; AIGC text MVP still uses supervised proxy.",
        },
        "ref_magvar": ref_pack,
        "classified": classified,
        "contingency_2x2": {
            "labels": ["benign", "malignant"],
            "matrix": [
                [contingency["tn"], contingency["fp"]],
                [contingency["fn"], contingency["tp"]],
            ],
            **contingency,
        },
        "accuracy": float(contingency["tn"] + contingency["tp"]) / 2.0,
    }


# ---------------------------------------------------------------------------
# Multimodal attribution (group-permute on V/A/T)
# ---------------------------------------------------------------------------


def multimodal_group_permute_attr(
    seed: int = SEED,
    n_estimators: int = 60,
    n_perm: int = 4,
) -> dict:
    """MSR-VTT-like hop: Mag×Var on PermuCATE unit risk + modality shares."""
    stream = make_mm_hop_proxy(n_batches=6, n_per=140, p=30, n_arms=4, seed=seed)
    groups = stream.feature_groups  # video/audio/text slices
    # early batch = REF, late concept-flip batch = LIVE
    b_ref = 1
    b_live = 4
    m0 = stream.batch == b_ref
    m1 = stream.batch == b_live
    X0, y0 = stream.X[m0], stream.R[m0]
    X1, y1 = stream.X[m1], stream.R[m1]

    X = np.vstack([X0, X1])
    Y = np.concatenate([y0, y1])
    W = np.concatenate([np.zeros(len(X0), dtype=int), np.ones(len(X1), dtype=int)])
    fit = permucate_unit_po(X, Y, W, seed=seed, n_estimators=n_estimators)
    unit_live = fit["unit_po"][W == 1]
    mv = mag_var_from_unit(unit_live, K=28, seed=seed + 2)
    base_r = float(fit["r_po"])

    rng = np.random.default_rng(int(seed) + 11)
    live_idx = np.asarray(np.where(W == 1)[0], dtype=int).reshape(-1)
    deltas = {}
    for name, sl in groups.items():
        cols = np.asarray(np.arange(X.shape[1])[sl], dtype=int).reshape(-1)
        bumps = []
        for t in range(int(n_perm)):
            Xp = np.array(X, copy=True)
            block = Xp[live_idx][:, cols]
            Xp[live_idx[:, None], cols] = block[rng.permutation(len(live_idx))]
            fit_p = permucate_unit_po(
                Xp, Y, W, seed=seed + 100 + t, n_estimators=max(40, n_estimators // 2)
            )
            # risk lift under modality destroy — primary attribution signal
            bumps.append(float(fit_p["r_po"]) - base_r)
        deltas[name] = float(np.mean(bumps))

    # use positive risk lifts; if all ≤0 fall back to |Δ|
    pos = {k: max(v, 0.0) for k, v in deltas.items()}
    if sum(pos.values()) <= 1e-12:
        pos = {k: abs(v) for k, v in deltas.items()}
    tot = sum(pos.values()) + 1e-12
    shares = {k: pos[k] / tot for k in pos}
    vimp = fit["vimp"]
    vimp_share = {}
    for name, sl in groups.items():
        cols = np.arange(len(vimp))[sl]
        vimp_share[name] = float(np.maximum(vimp[cols], 0.0).sum())
    vs = sum(vimp_share.values()) + 1e-12
    vimp_share = {k: vimp_share[k] / vs for k in vimp_share}

    return {
        "board": "mm_hop_proxy",
        "hop": {"ref_batch": b_ref, "live_batch": b_live},
        "magvar": mv,
        "r_po": fit["r_po"],
        "group_permute_delta_u": deltas,
        "modality_shares_permute": shares,
        "modality_shares_vimp": vimp_share,
        "justify": (
            "Mag×Var on PermuCATE unit risk diagnoses shift heterogeneity; "
            "group-permute / VIMP mass attributes which modality block carries that risk."
        ),
        "estimator": fit["estimator"],
    }


def plot_mvp(payload: dict, path: Path):
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12.0, 8.2))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # 1) Mag×Var scatter
    ax = fig.add_subplot(gs[0, 0])
    ref = payload["magvar_eval"]["ref_magvar"]
    colors = {"benign": "#2F6B4F", "malignant": "#A33B24"}
    ax.scatter(ref["mu"], ref["s"], s=140, c="#5B6B7A", label="REF", zorder=3)
    for k, row in payload["magvar_eval"]["classified"].items():
        ax.scatter(row["mu"], row["s"], s=160, c=colors[k], label=k, zorder=3)
        ax.annotate(k, (row["mu"], row["s"]), fontsize=8, xytext=(6, 4), textcoords="offset points")
    ax.set_xlabel("Magnitude μ = E[(φ−τ̂)²]")
    ax.set_ylabel("Variance s = Std(r_b)")
    ax.set_title("PermuCATE RF Mag×Var", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # 2) 2×2
    ax = fig.add_subplot(gs[0, 1])
    M = np.asarray(payload["magvar_eval"]["contingency_2x2"]["matrix"], dtype=float)
    im = ax.imshow(M, cmap="YlGnBu", vmin=0, vmax=max(M.max(), 1))
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["pred benign", "pred malignant"])
    ax.set_yticklabels(["true benign", "true malignant"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, "%d" % int(M[i, j]), ha="center", va="center", fontsize=16, fontweight="bold")
    ax.set_title("2×2 (PermuCATE RF unit risk)", loc="left", fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046)

    # 3) multimodal shares
    ax = fig.add_subplot(gs[1, 0])
    mm = payload["multimodal"]
    names = list(mm["modality_shares_permute"].keys())
    x = np.arange(len(names))
    ax.bar(x - 0.18, [mm["modality_shares_permute"][k] for k in names], width=0.36, color="#1F6F8B", label="group-permute")
    ax.bar(x + 0.18, [mm["modality_shares_vimp"][k] for k in names], width=0.36, color="#C4892A", alpha=0.8, label="τ̂ VIMP")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("share")
    ax.set_title("Multimodal attribution (V/A/T)", loc="left", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, axis="y", alpha=0.3)

    # 4) estimator callout
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    est = payload["magvar_eval"]["estimator"]
    txt = [
        "ESTIMATOR (this MVP)",
        "  model_m  = %s" % est["model_m"],
        "  model_e  = %s" % est["model_e"],
        "  model_tau= %s" % est["model_tau"],
        "  unit_po  = %s" % est["unit_po"],
        "  used_permucate_rf = True",
        "",
        "2×2 acc = %.2f" % payload["magvar_eval"]["accuracy"],
        "J_ben=%.3f  J_mal=%.3f"
        % (
            payload["magvar_eval"]["classified"]["benign"]["J"],
            payload["magvar_eval"]["classified"]["malignant"]["J"],
        ),
        "",
        "MM: Mag×Var = shift type;",
        "    permute/VIMP = which modality.",
    ]
    ax.text(
        0.02,
        0.98,
        "\n".join(txt),
        transform=ax.transAxes,
        va="top",
        fontsize=8.2,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F7F4EF", edgecolor="#CCC"),
    )
    ax.set_title("Professor MVP card", loc="left", fontweight="bold")

    fig.suptitle(
        "PermuCATE RF Mag×Var + multimodal attribution MVP",
        fontsize=12.5,
        fontweight="bold",
        y=0.995,
    )
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items() if k != "fit"}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, float)):
        return float(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    return obj


def run_mvp(
    seed: int = SEED,
    out_dir: Optional[Path] = None,
    n_estimators: int = 70,
) -> dict:
    out_dir = Path(out_dir or OUT_DEFAULT)
    out_dir.mkdir(parents=True, exist_ok=True)
    mag = eval_magvar_permucate(seed=seed, n_estimators=n_estimators)
    # strip heavy arrays from nested fit for JSON
    for k in list(mag.get("classified", {})):
        pass
    mm = multimodal_group_permute_attr(seed=seed, n_estimators=max(40, n_estimators // 2))
    payload = {
        "mvp": "permucate_rf_magvar_multimodal",
        "magvar_eval": _jsonable(mag),
        "multimodal": _jsonable(mm),
        "professor_pitch": (
            "Mag×Var on PermuCATE RF unit risk (φ−τ̂)² resolves benign vs malignant "
            "shift (2×2); multimodal group-permute attributes the channel."
        ),
    }
    png = out_dir / "permucate_magvar_mm.png"
    plot_mvp(payload, png)
    (out_dir / "permucate_magvar_mm.json").write_text(json.dumps(payload, indent=2))
    _write_readme(out_dir / "README.md", payload, png.name)
    return payload


def _write_readme(path: Path, payload: dict, png_name: str):
    mag = payload["magvar_eval"]
    c = mag["contingency_2x2"]
    mm = payload["multimodal"]
    lines = [
        "# PermuCATE RF Mag×Var + multimodal attribution MVP",
        "",
        "## Estimator (important)",
        "",
        "- **Yes — RandomForest outcome models** (μ̂₀, μ̂₁), RF propensity, RF τ̂ on DR φ.",
        "- Unit PO-risk $u_i=(\\varphi_i-\\hat\\tau_i)^2$ (same family as `permuCATE_vimp.po_risk`).",
        "- Earlier AIGC text boards used supervised proxy $1-\\hat p(Y\\mid x)$ — different path.",
        "",
        "## 2×2",
        "",
        "|  | pred benign | pred malignant |",
        "| --- | --- | --- |",
        "| true benign | %d | %d |" % (c["tn"], c["fp"]),
        "| true malignant | %d | %d |" % (c["fn"], c["tp"]),
        "",
        "accuracy = %.2f" % mag["accuracy"],
        "",
        "## Multimodal shares",
        "",
        "| modality | group-permute | VIMP |",
        "| --- | --- | --- |",
    ]
    for k in mm["modality_shares_permute"]:
        lines.append(
            "| %s | %.3f | %.3f |"
            % (k, mm["modality_shares_permute"][k], mm["modality_shares_vimp"][k])
        )
    lines.extend(
        [
            "",
            payload["professor_pitch"],
            "",
            f"![]({png_name})",
            "",
            "```bash",
            "python3 scripts/run_as_dta_permucate_magvar.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines))
