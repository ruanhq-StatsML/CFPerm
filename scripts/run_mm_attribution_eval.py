#!/usr/bin/env python3
"""Next-stage multimodal attribution evaluation.

Core question: when the batch shift is *known* to live in one modality,
do RF Domain / coord-MMD / PO-risk put VIMP mass on that modality?

Protocol
--------
1. Observational controls (known inductive bias)
   - Indiana Frontal/Lateral  → GT = image
   - COCO keyword outdoor/in  → GT = text   (constructed; optional)
   - COCO early/late image_id → weak/null observational baseline

2. Controlled interventions on COCO early/late (time-order base)
   Keep the same early/late split, then inject a mean shift into the
   *late* batch only:
   - image_inject : X_late[:, :d_img] += α · u_img
   - text_inject  : X_late[:, d_img:] += α · u_txt
   GT modality = the injected block.

3. Metrics (all use ordered scores VIMP)
   - modality_mass_share on image/text
   - selection_AUC treating GT-block coords as positives
   - topk_precision: fraction of global top-k in GT block (k=|GT|=512
     for full-block; also report k=20)
   - domain_auc for RF (shift strength)

  python3 scripts/run_mm_attribution_eval.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor" / "fsds"))
from VIMP_mmd_benchmark import MMD, _subsample_batch  # noqa: E402

DATA = ROOT / "data" / "img_txt"
OUT = ROOT / "results" / "mm_attribution_eval"
DOCS = ROOT / "docs" / "method"

SEED = 2026
D_IMG = 512
D_TXT = 512
TOP_K = 20
MAX_MMD_N = 100
N_SUB = 4000
ALPHA_CLIP = 0.01
INJECT_ALPHA = 0.75  # mean-shift strength in late batch


def metrics_vs_block(scores: np.ndarray, gt: str, *, k_small: int = TOP_K) -> dict:
    scores = np.asarray(scores, float)
    p = len(scores)
    if gt == "image":
        pos = np.arange(0, D_IMG)
    elif gt == "text":
        pos = np.arange(D_IMG, D_IMG + D_TXT)
    else:
        raise ValueError(gt)
    membership = np.zeros(p, int)
    membership[pos] = 1
    auc = float(roc_auc_score(membership, scores))
    # full-block top-|S|
    k_full = len(pos)
    sel_full = set(np.argsort(-scores)[:k_full].tolist())
    tp = len(sel_full & set(pos.tolist()))
    prec_full = tp / k_full
    # small top-k
    sel_k = np.argsort(-scores)[:k_small]
    prec_k = float(np.mean([i in set(pos.tolist()) for i in sel_k]))
    mass = {
        "image": float(scores[:D_IMG].sum()),
        "text": float(scores[D_IMG:].sum()),
    }
    tot = mass["image"] + mass["text"] + 1e-12
    share = {m: mass[m] / tot for m in mass}
    return {
        "gt": gt,
        "selection_AUC": round(auc, 4),
        "topk_precision_fullblock": round(prec_full, 4),
        "topk_precision_k20": round(prec_k, 4),
        "modality_vimp_share": share,
        "mass_on_gt": round(share[gt], 4),
    }


def rf_domain(X0, X1, *, seed=SEED):
    X = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    p = X.shape[1]
    clf = RandomForestClassifier(
        n_estimators=150,
        max_depth=max(2, int(round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(round(np.sqrt(len(X)) // 2))),
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    Xtr, Xte, Wtr, Wte = train_test_split(X, W, test_size=0.25, random_state=seed, stratify=W)
    clf2 = RandomForestClassifier(
        n_estimators=120,
        max_depth=max(2, int(round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(round(np.sqrt(len(X)) // 2))),
        n_jobs=-1,
        random_state=seed + 1,
    )
    clf2.fit(Xtr, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return vimp, auc


def coord_mmd(X0, X1, *, max_n=MAX_MMD_N, seed=SEED):
    rng = np.random.default_rng(seed)
    A = _subsample_batch(X0, max_n, rng)
    B = _subsample_batch(X1, max_n, rng)
    mmd = MMD(compute_kernel="rbf")
    uni = np.zeros(A.shape[1])
    for j in range(A.shape[1]):
        uni[j], _ = mmd(A[:, [j]], B[:, [j]])
    return uni


def po_risk(X, Y, W, *, seed=SEED):
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    Yf = Y.astype(float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=60, max_depth=10, min_samples_leaf=3, random_state=seed + fold, n_jobs=-1
        )
        e = RandomForestClassifier(
            n_estimators=60, max_depth=10, min_samples_leaf=3, random_state=seed + 40 + fold, n_jobs=-1
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1 - ALPHA_CLIP)
    po = (Yf - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=150, max_depth=12, min_samples_leaf=3, random_state=seed + 7, n_jobs=-1
    )
    tau.fit(X, po)
    return tau.feature_importances_.astype(float), float(np.mean(tau.predict(X) ** 2))


def subsample_pair(X0, Y0, X1, Y1, *, n_max=N_SUB, seed=SEED):
    rng = np.random.default_rng(seed)
    n0 = min(len(X0), n_max // 2)
    n1 = min(len(X1), n_max // 2)
    i0 = rng.choice(len(X0), n0, replace=False)
    i1 = rng.choice(len(X1), n1, replace=False)
    return X0[i0], Y0[i0], X1[i1], Y1[i1]


def load_img_txt(name: str):
    d = DATA / name
    if name == "fashion_iq":
        df1 = np.load(d / "df1.npy")
        df2 = np.load(d / "df2.npy")
        X0, Y0 = df1[:, :-1], df1[:, -1]
        X1, Y1 = df2[:, :-1], df2[:, -1]
        return X0, Y0, X1, Y1, "train vs test", None
    img = np.load(d / "img_feats.npy")
    txt = np.load(d / "txt_feats.npy")
    y = np.load(d / "labels.npy").astype(float).reshape(-1)
    meta = pd.read_csv(d / "df_metadata.csv")
    n = min(len(img), len(txt), len(y), len(meta))
    X = np.hstack([img[:n].astype(float), txt[:n].astype(float)])
    y = y[:n]
    meta = meta.iloc[:n]
    if name == "indiana_cxr":
        w = (~meta["projection"].astype(str).str.lower().str.startswith("front")).astype(int).to_numpy()
        batch = "Frontal vs Lateral"
        gt = "image"
    elif name == "coco_outdoor_indoor":
        w = meta["batch"].astype(int).to_numpy()
        batch = "outdoor/indoor keywords"
        gt = "text"
    elif name == "coco_time_order":
        w = meta["batch"].astype(int).to_numpy()
        batch = "early vs late image_id"
        gt = None  # observational weak/null
    else:
        raise ValueError(name)
    return X[w == 0], y[w == 0], X[w == 1], y[w == 1], batch, gt


def inject_shift(X0, X1, *, which: str, alpha: float = INJECT_ALPHA, seed: int = SEED):
    """Add a direction mean-shift to late batch only."""
    rng = np.random.default_rng(seed)
    X1 = X1.copy()
    if which == "image":
        u = rng.normal(size=D_IMG)
        u /= np.linalg.norm(u) + 1e-12
        X1[:, :D_IMG] += alpha * u
        gt = "image"
    elif which == "text":
        u = rng.normal(size=D_TXT)
        u /= np.linalg.norm(u) + 1e-12
        X1[:, D_IMG:] += alpha * u
        gt = "text"
    else:
        raise ValueError(which)
    return X0, X1, gt


def run_scenario(name: str, X0, Y0, X1, Y1, *, batch: str, gt: str | None, tag: str):
    print(f"\n==== {tag} · GT={gt} · {batch} ====", flush=True)
    X0s, Y0s, X1s, Y1s = subsample_pair(X0, Y0, X1, Y1)
    out = {
        "scenario": tag,
        "batch": batch,
        "gt_modality": gt,
        "n0": int(len(X0)),
        "n1": int(len(X1)),
        "n0_sub": int(len(X0s)),
        "n1_sub": int(len(X1s)),
        "methods": {},
    }

    t0 = time.perf_counter()
    rf_vimp, rf_auc = rf_domain(X0s, X1s)
    out["methods"]["rf_domain"] = {
        "domain_auc": round(rf_auc, 4),
        "time_s": round(time.perf_counter() - t0, 2),
    }
    if gt:
        out["methods"]["rf_domain"].update(metrics_vs_block(rf_vimp, gt))
    else:
        share = metrics_vs_block(rf_vimp, "image")["modality_vimp_share"]
        out["methods"]["rf_domain"]["modality_vimp_share"] = share
    print(f"  RF auc={rf_auc:.3f} → {out['methods']['rf_domain']}", flush=True)

    t0 = time.perf_counter()
    mmd_vimp = coord_mmd(X0, X1)
    out["methods"]["coord_mmd"] = {"time_s": round(time.perf_counter() - t0, 2)}
    if gt:
        out["methods"]["coord_mmd"].update(metrics_vs_block(mmd_vimp, gt))
    else:
        out["methods"]["coord_mmd"]["modality_vimp_share"] = metrics_vs_block(mmd_vimp, "image")[
            "modality_vimp_share"
        ]
    print(f"  MMD → {out['methods']['coord_mmd']}", flush=True)

    X = np.vstack([X0s, X1s])
    Y = np.concatenate([Y0s, Y1s])
    W = np.concatenate([np.zeros(len(X0s)), np.ones(len(X1s))]).astype(int)
    t0 = time.perf_counter()
    po_vimp, por = po_risk(X, Y, W)
    out["methods"]["po_risk"] = {
        "po_risk": round(por, 4),
        "time_s": round(time.perf_counter() - t0, 2),
    }
    if gt:
        out["methods"]["po_risk"].update(metrics_vs_block(po_vimp, gt))
    else:
        out["methods"]["po_risk"]["modality_vimp_share"] = metrics_vs_block(po_vimp, "image")[
            "modality_vimp_share"
        ]
    print(f"  PO → {out['methods']['po_risk']}", flush=True)
    return out


def write_artifacts(rows: list[dict]):
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    (OUT / "mm_attribution_eval.json").write_text(json.dumps(rows, indent=2))

    # markdown protocol + results
    md = [
        "# Multimodal attribution evaluation\n\n",
        "## How to evaluate the logic\n\n",
        "1. **Define GT modality from the batch mechanism**, not from the scores.\n",
        "   - view shift (Frontal/Lateral) → GT = image\n",
        "   - caption-keyword split → GT = text (constructed control)\n",
        "   - time-order only → no strong GT (baseline / null)\n",
        "   - synthetic inject on one block → GT = that block\n\n",
        "2. **Score methods** with RF Domain (CS), coord-MMD (CS), PO-risk (CD).\n\n",
        "3. **Judge recovery** via\n",
        "   - `mass_on_gt` (modality VIMP share on GT)\n",
        "   - `selection_AUC` (coords in GT block as positives)\n",
        "   - `topk_precision_k20` / full-block precision\n",
        "   - RF `domain_auc` (was the shift actually detectable?)\n\n",
        "4. **Pass criteria (practical)**\n",
        "   - strong CS controls: `mass_on_gt ≥ 0.6` and `selection_AUC ≥ 0.6`\n",
        "   - inject recoveries: same, with `domain_auc` clearly > 0.5\n",
        "   - time-order baseline: near-balanced mass OK if AUC≈0.5–0.6\n\n",
        "## Results\n\n",
        "| Scenario | GT | Method | domain AUC / PO | mass_on_gt | sel AUC | P@20 |\n",
        "|----------|----|--------|-----------------|------------|---------|------|\n",
    ]
    for r in rows:
        gt = r["gt_modality"] or "—"
        for method, m in r["methods"].items():
            strength = m.get("domain_auc", m.get("po_risk", "—"))
            mass = m.get("mass_on_gt", "—")
            auc = m.get("selection_AUC", "—")
            p20 = m.get("topk_precision_k20", "—")
            md.append(
                f"| {r['scenario']} | {gt} | {method} | {strength} | {mass} | {auc} | {p20} |\n"
            )
    (OUT / "README.md").write_text("".join(md))

    # LaTeX compact recovery table (scenarios with GT only)
    lines = [
        "% Multimodal attribution evaluation · recovery under known GT modality\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Multimodal attribution recovery. Ground-truth modality is fixed by "
        "the batch mechanism (observational control or synthetic late-batch inject). "
        "Metrics: modality VIMP mass on GT, selection AUC over GT coordinates, "
        "and precision of global top-20.}\n",
        "\\label{tab:mm-attr-eval-recovery}\n\\small\n",
        "\\begin{tabular}{@{}l l l c c c@{}}\n\\toprule\n",
        "Scenario & GT & Method & Mass on GT & Sel.\\ AUC & P@20 \\\\\n\\midrule\n",
    ]
    for r in rows:
        if not r["gt_modality"]:
            continue
        for method in ("rf_domain", "coord_mmd", "po_risk"):
            m = r["methods"][method]
            lab = {"rf_domain": "RF Domain", "coord_mmd": "Coord-MMD", "po_risk": "PO-risk"}[method]
            lines.append(
                f"{r['scenario']} & {r['gt_modality']} & {lab} & "
                f"${m['mass_on_gt']:.2f}$ & ${m['selection_AUC']:.2f}$ & "
                f"${m['topk_precision_k20']:.2f}$ \\\\\n"
            )
    lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    # baseline time-order mass table
    lines.append("\\begin{table}[ht]\n\\centering\n")
    lines.append(
        "\\caption{Time-order COCO baseline (early vs late image\\_id): "
        "no strong modality GT. Report modality mass only.}\n"
    )
    lines.append("\\label{tab:mm-attr-eval-time-baseline}\n\\small\n")
    lines.append("\\begin{tabular}{@{}l cc c@{}}\n\\toprule\n")
    lines.append("Method & Image share & Text share & RF AUC / PO \\\\\n\\midrule\n")
    for r in rows:
        if r["scenario"] != "coco_time_order_baseline":
            continue
        for method, lab in [
            ("rf_domain", "RF Domain"),
            ("coord_mmd", "Coord-MMD"),
            ("po_risk", "PO-risk"),
        ]:
            m = r["methods"][method]
            s = m["modality_vimp_share"]
            if "domain_auc" in m:
                strength = f"{m['domain_auc']:.3f}"
            elif "po_risk" in m:
                strength = f"{m['po_risk']:.3f}"
            else:
                strength = "---"
            lines.append(
                f"{lab} & ${s['image']:.2f}$ & ${s['text']:.2f}$ & ${strength}$ \\\\\n"
            )
    lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    tex = "".join(lines)
    (DOCS / "Multimodal_Attribution_Eval_tables_only.tex").write_text(tex)
    (OUT / "Multimodal_Attribution_Eval_tables_only.tex").write_text(tex)
    print(tex)
    print("".join(md))


def main():
    rows = []

    # --- observational controls ---
    for name, tag in [
        ("indiana_cxr", "indiana_view_control"),
        ("coco_outdoor_indoor", "coco_keyword_control"),
        ("coco_time_order", "coco_time_order_baseline"),
    ]:
        X0, Y0, X1, Y1, batch, gt = load_img_txt(name)
        rows.append(run_scenario(name, X0, Y0, X1, Y1, batch=batch, gt=gt, tag=tag))

    # --- controlled inject on time-order COCO ---
    X0, Y0, X1, Y1, batch, _ = load_img_txt("coco_time_order")
    for which in ("image", "text"):
        X0i, X1i, gt = inject_shift(X0, X1, which=which)
        rows.append(
            run_scenario(
                "coco_time_order",
                X0i,
                Y0,
                X1i,
                Y1,
                batch=f"{batch} + late {which} inject α={INJECT_ALPHA}",
                gt=gt,
                tag=f"coco_time_inject_{which}",
            )
        )

    write_artifacts(rows)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
