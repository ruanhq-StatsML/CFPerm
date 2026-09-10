#!/usr/bin/env python3
"""AFFEC FSDS · separate covariate-shift vs concept-drift rankings.

CS: RF Domain Classifier VIMP for W ~ X  → np.argsort(-VIMP)[:k]
CD: RF VIMP on PO pseudo-outcome / within-batch Y~X importance gap → np.argsort(-VIMP)[:k]

Uses cached XYW from results/affec_fsds/affec_fsds_xyw_cache.npz

  python3 scripts/run_affec_fsds_cs_cd_rank.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "affec_fsds"
DOCS = ROOT / "docs" / "method"
CACHE = OUT / "affec_fsds_xyw_cache.npz"

SEED = 2026
TOP_K = 8
ALPHA_CLIP = 0.01

CHANNEL_NAMES = {
    "eye_tracking": {
        0: "fixation_x",
        1: "fixation_y",
        2: "gaze_x_left",
        3: "gaze_y_left",
        4: "gaze_x_right",
        5: "gaze_y_right",
        6: "fixation_duration",
        7: "blink_left",
        8: "blink_right",
        9: "eye_open_left",
        10: "eye_open_right",
        11: "pupil_x_left",
        12: "pupil_y_left",
        13: "pupil_x_right",
        14: "pupil_y_right",
        15: "validity",
    },
    "pupil": {
        0: "pupil_diameter",
        1: "pupil_diameter_raw",
        2: "pupil_diameter_filt",
        3: "pupil_baseline",
        4: "pupil_peak",
        5: "pupil_latency",
        6: "pupil_velocity",
        7: "pupil_acceleration",
        8: "pupil_variability",
        9: "pupil_slope",
        10: "eye_pos_x",
        11: "eye_pos_y",
        12: "eye_pos_x_raw",
        13: "eye_pos_y_raw",
        14: "eye_pos_velocity_x",
        15: "eye_pos_velocity_y",
        16: "eye_pos_dispersion",
        17: "eye_pos_range_x",
        18: "eye_pos_range_y",
        19: "eye_pos_slope",
        20: "pupil_validity",
    },
    "cursor": {
        0: "cursor_x",
        1: "cursor_y",
        2: "cursor_velocity",
        3: "cursor_state",
    },
    "gsr_eda": {
        0: "gsr_raw",
        1: "gsr_filtered",
        2: "gsr_phasic",
        3: "gsr_tonic",
        4: "gsr_scr_count",
        5: "gsr_peaks",
        6: "gsr_slope",
        7: "gsr_mean",
        8: "gsr_std",
        9: "gsr_range",
        10: "body_temp",
        11: "temp_mean",
        12: "temp_std",
        13: "temp_slope",
        14: "acc_x",
        15: "acc_y",
        16: "acc_z",
        17: "acc_magnitude",
        18: "acc_std",
        19: "acc_slope",
        **{i: f"gsr_feat_{i}" for i in range(20, 40)},
    },
    "eeg": {
        **{0: "Fp1", 1: "Fp2", 2: "F3", 3: "F4", 4: "C3", 5: "C4", 6: "P3", 7: "P4", 8: "O1", 9: "O2"},
        **{i: f"EEG_{i}" for i in range(10, 63)},
    },
}


def _name(mod: str, idx: int) -> str:
    return CHANNEL_NAMES.get(mod, {}).get(idx, f"feat_{idx}")


def ranked_by_modality(vimp: np.ndarray, block_slices: dict, mods: list[str], *, k: int = TOP_K):
    out, detail = {}, {}
    for m in mods:
        a, b = block_slices[m]
        vm = np.asarray(vimp[a:b], float)
        order = [int(i) for i in np.argsort(-vm)[: min(k, len(vm))]]
        out[m] = order
        detail[m] = [
            {
                "rank": r + 1,
                "index": idx,
                "global_index": int(a + idx),
                "name": _name(m, idx),
                "vimp": float(vm[idx]),
            }
            for r, idx in enumerate(order)
        ]
    return out, detail


def covariate_shift_vimp(X, W, *, seed: int = SEED):
    """CS ranking signal: RF Domain Classifier · W ~ X."""
    clf = RandomForestClassifier(
        n_estimators=200, max_depth=10, min_samples_leaf=5, random_state=seed, n_jobs=-1
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    rng = np.random.default_rng(seed)
    idx = np.arange(len(W))
    rng.shuffle(idx)
    cut = int(0.8 * len(W))
    tr, te = idx[:cut], idx[cut:]
    clf2 = RandomForestClassifier(
        n_estimators=120, max_depth=10, min_samples_leaf=5, random_state=seed + 1, n_jobs=-1
    )
    clf2.fit(X[tr], W[tr])
    auc = float(roc_auc_score(W[te], clf2.predict_proba(X[te])[:, 1]))
    return vimp, {"rf_domain_auc": auc}


def concept_drift_vimp(X, Y, W, *, seed: int = SEED):
    """CD ranking signal: features that move P(Y|X) across batches.

    1) Cross-fit m(X), e(X); po = (Y-m)(W-e)
    2) RF regressor VIMP for po ~ X  (PO-path drivers)
    3) Blend with |imp_Y|W=1 - imp_Y|W=0| from within-batch Y~X RFs
    """
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    W = np.asarray(W, int)
    n, p = X.shape
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=80, max_depth=8, min_samples_leaf=5, random_state=seed + fold, n_jobs=-1
        )
        e = RandomForestClassifier(
            n_estimators=80, max_depth=8, min_samples_leaf=5, random_state=seed + 50 + fold, n_jobs=-1
        )
        m.fit(X[tr], Y[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1.0 - ALPHA_CLIP)
    po = (Y - m_hat) * (W - e_hat)

    tau = RandomForestRegressor(
        n_estimators=200, max_depth=10, min_samples_leaf=5, random_state=seed + 7, n_jobs=-1
    )
    tau.fit(X, po)
    vimp_po = tau.feature_importances_.astype(float)
    po_risk = float(np.mean(tau.predict(X) ** 2))

    # within-batch Y~X importance gap (concept change of which X drive Y)
    y0 = RandomForestRegressor(
        n_estimators=150, max_depth=10, min_samples_leaf=5, random_state=seed + 11, n_jobs=-1
    )
    y1 = RandomForestRegressor(
        n_estimators=150, max_depth=10, min_samples_leaf=5, random_state=seed + 13, n_jobs=-1
    )
    y0.fit(X[W == 0], Y[W == 0])
    y1.fit(X[W == 1], Y[W == 1])
    vimp_gap = np.abs(y1.feature_importances_ - y0.feature_importances_).astype(float)

    # primary CD score: PO VIMP, light blend with Y~X gap
    vimp = 0.7 * (vimp_po / (vimp_po.sum() + 1e-12)) + 0.3 * (vimp_gap / (vimp_gap.sum() + 1e-12))
    return vimp, {"po_risk": po_risk, "vimp_po": vimp_po, "vimp_gap": vimp_gap}


def write_latex(cs_fi, cs_det, cd_fi, cd_det, meta, path: Path):
    def rows(fi, det):
        lines = []
        for m, idxs in fi.items():
            names = ", ".join(d["name"].replace("_", "\\_") for d in det[m])
            idx_s = "[" + ",".join(str(i) for i in idxs) + "]"
            mod_tex = m.replace("_", "\\_")
            lines.append(f"\\texttt{{{mod_tex}}} & ${idx_s}$ & {names} \\\\")
        return "\n".join(lines)

    k = TOP_K
    auc = meta["rf_domain_auc"]
    por = meta["po_risk"]
    tex = (
        "% AFFEC FSDS · covariate-shift vs concept-drift rankings (RF ordered)\n"
        "\\begin{table}[ht]\n"
        "\\centering\n"
        "\\caption{Covariate-shift ranking on AFFEC. RF Domain Classifier VIMP for $W\\sim X$,\n"
        "ordered by $\\texttt{np.argsort(-VIMP)[:k]}$ with $k="
        + str(k)
        + "$. Batch: early runs $\\{0,1\\}$ vs late runs $\\{2,3\\}$.\n"
        f"Domain AUC ${auc:.3f}$.}}\n"
        "\\label{tab:affec-fsds-cs-rank}\n"
        "\\small\n"
        "\\begin{tabular}{@{}l l p{7.2cm}@{}}\n"
        "\\toprule\n"
        "Modality & Ranked indices & Channel names (same order) \\\\\n"
        "\\midrule\n"
        f"{rows(cs_fi, cs_det)}\n"
        "\\bottomrule\n"
        "\\end{tabular}\n"
        "\\end{table}\n\n"
        "\\begin{table}[ht]\n"
        "\\centering\n"
        "\\caption{Concept-drift ranking on AFFEC. PO-risk path RF VIMP on pseudo-outcome\n"
        "(blended with within-batch $Y\\sim X$ importance gap), ordered by "
        "$\\texttt{np.argsort(-VIMP)[:k]}$ with $k="
        + str(k)
        + f"$. Observed PO-risk ${por:.4f}$.}}\n"
        "\\label{tab:affec-fsds-cd-rank}\n"
        "\\small\n"
        "\\begin{tabular}{@{}l l p{7.2cm}@{}}\n"
        "\\toprule\n"
        "Modality & Ranked indices & Channel names (same order) \\\\\n"
        "\\midrule\n"
        f"{rows(cd_fi, cd_det)}\n"
        "\\bottomrule\n"
        "\\end{tabular}\n"
        "\\end{table}\n"
    )
    path.write_text(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    if not CACHE.exists():
        raise SystemExit(f"missing cache {CACHE}; run scripts/run_affec_fsds.py first")

    z = np.load(CACHE, allow_pickle=True)
    X, Y, W = z["X"], z["Y"], z["W"]
    block_slices = {k: tuple(v) for k, v in z["block_slices"].item().items()}
    mods = list(z["mods"])
    print(f"cache n={len(W)} p={X.shape[1]}", flush=True)

    print("covariate-shift ranking (RF Domain Classifier)…", flush=True)
    cs_vimp, cs_meta = covariate_shift_vimp(X, W, seed=SEED)
    cs_fi, cs_det = ranked_by_modality(cs_vimp, block_slices, mods, k=TOP_K)
    print("CS:", cs_fi, flush=True)

    print("concept-drift ranking (PO-risk RF + Y~X gap)…", flush=True)
    cd_vimp, cd_meta = concept_drift_vimp(X, Y, W, seed=SEED)
    cd_fi, cd_det = ranked_by_modality(cd_vimp, block_slices, mods, k=TOP_K)
    print("CD:", cd_fi, flush=True)

    payload = {
        "batch_def": "W=0: run∈{0,1}; W=1: run∈{2,3}",
        "top_k": TOP_K,
        "n": int(len(W)),
        "p": int(X.shape[1]),
        "covariate_shift": {
            "method": "RF Domain Classifier · W ~ X · argsort(-VIMP)[:k]",
            "rf_domain_auc": cs_meta["rf_domain_auc"],
            "feature_indices": cs_fi,
            "feature_indices_detail": cs_det,
        },
        "concept_drift": {
            "method": "PO-risk RF VIMP on pseudo-outcome + |imp_Y|W=1 - imp_Y|W=0| · argsort(-VIMP)[:k]",
            "po_risk": cd_meta["po_risk"],
            "feature_indices": cd_fi,
            "feature_indices_detail": cd_det,
        },
    }
    dest = OUT / "affec_fsds_cs_cd_ranked_indices.json"
    dest.write_text(json.dumps(payload, indent=2))
    print(f"wrote {dest}", flush=True)

    proto = OUT / "affec_fsds_cs_cd_ranked_indices_prototype.py"
    proto.write_text(
        "# Covariate-shift ranking (RF Domain Classifier, high→low)\n"
        f"feature_indices_covariate_shift = {repr(cs_fi)}\n\n"
        "# Concept-drift ranking (PO-risk path, high→low)\n"
        f"feature_indices_concept_drift = {repr(cd_fi)}\n"
    )
    print(f"wrote {proto}", flush=True)

    md_lines = [
        "# AFFEC FSDS · covariate-shift vs concept-drift rankings",
        "",
        "## Covariate shift (RF Domain Classifier · `argsort(-VIMP)[:k]`)",
        "",
        "```python",
        f"feature_indices_covariate_shift = {cs_fi}",
        "```",
        "",
        f"- domain AUC: {cs_meta['rf_domain_auc']:.3f}",
        "",
        "| modality | ranked indices | names |",
        "|---|---|---|",
    ]
    for m in mods:
        names = ", ".join(d["name"] for d in cs_det[m])
        md_lines.append(f"| {m} | {cs_fi[m]} | {names} |")
    md_lines += [
        "",
        "## Concept drift (PO-risk RF · `argsort(-VIMP)[:k]`)",
        "",
        "```python",
        f"feature_indices_concept_drift = {cd_fi}",
        "```",
        "",
        f"- PO-risk: {cd_meta['po_risk']:.4f}",
        "",
        "| modality | ranked indices | names |",
        "|---|---|---|",
    ]
    for m in mods:
        names = ", ".join(d["name"] for d in cd_det[m])
        md_lines.append(f"| {m} | {cd_fi[m]} | {names} |")
    md = OUT / "affec_fsds_cs_cd_ranked_indices.md"
    md.write_text("\n".join(md_lines) + "\n")
    print(f"wrote {md}", flush=True)

    tex = DOCS / "AFFEC_FSDS_CS_CD_RankedIndices_tables_only.tex"
    write_latex(
        cs_fi,
        cs_det,
        cd_fi,
        cd_det,
        {"rf_domain_auc": cs_meta["rf_domain_auc"], "po_risk": cd_meta["po_risk"]},
        tex,
    )
    print(f"wrote {tex}", flush=True)


if __name__ == "__main__":
    main()
