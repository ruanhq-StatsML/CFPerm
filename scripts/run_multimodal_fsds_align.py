#!/usr/bin/env python3
"""Multimodal FSDS alignment across datasets · RF / MMD / PO-risk.

Datasets (X=[img|txt], last-col Y when built as clip_embedding):
  fashion_iq          — CLIP fashion; W = train vs test
  indiana_cxr         — OpenI CXR+report; W = Frontal vs Lateral
  coco_outdoor_indoor — COCO CLIP; W = outdoor vs indoor caption keywords
  coco_time_order     — COCO CLIP; W = early vs late image_id (order proxy)
  microscopy_clip     — microscopy image/text CLIP; W = long vs short caption

  python3 scripts/run_multimodal_fsds_align.py
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

DATA_ROOT = ROOT / "data" / "img_txt"
OUT = ROOT / "results" / "multimodal_fsds_align"
DOCS = ROOT / "docs" / "method"

SEED = 2026
TOP_K = 20
MAX_MMD_N = 100
ALPHA_CLIP = 0.01
N_SUB_RF = 4000


def load_dataset(name: str):
    """Return df1, df2 as (n, p+1) with last col = Y, plus meta dict."""
    d = DATA_ROOT / name
    if name == "fashion_iq":
        df1 = np.load(d / "df1.npy")
        df2 = np.load(d / "df2.npy")
        return df1, df2, {
            "batch": "train vs test",
            "source": "Fashion-IQ CLIP zip",
            "d_img": 512,
            "d_txt": 512,
        }

    img = np.load(d / "img_feats.npy")
    txt = np.load(d / "txt_feats.npy")
    y = np.load(d / "labels.npy").astype(float).reshape(-1)
    meta = pd.read_csv(d / "df_metadata.csv")
    n = min(len(img), len(txt), len(y), len(meta))
    img, txt, y, meta = img[:n], txt[:n], y[:n], meta.iloc[:n].reset_index(drop=True)
    X = np.hstack([img.astype(float), txt.astype(float)])
    clip = np.hstack([X, y.reshape(-1, 1)])

    if name == "indiana_cxr":
        proj = meta["projection"].astype(str).str.lower()
        w = (~proj.str.startswith("front")).astype(int).to_numpy()
        batch = "Frontal vs Lateral"
        source = "Indiana/OpenI CXR + report CLIP"
    elif name == "coco_outdoor_indoor":
        w = meta["batch"].astype(int).to_numpy()
        batch = "outdoor vs indoor caption keywords"
        source = "COCO CLIP (open-clip ViT-B/32) outdoor/indoor keyword split"
    elif name == "coco_time_order":
        w = meta["batch"].astype(int).to_numpy()
        batch = "early vs late image_id (order)"
        source = "COCO CLIP (open-clip ViT-B/32) early/late image_id windows"
    elif name == "microscopy_clip":
        w = meta["batch"].astype(int).to_numpy()
        batch = "short vs long caption"
        source = "HF kvriza8 microscopy CLIP image+text"
    else:
        raise ValueError(name)

    df1, df2 = clip[w == 0], clip[w == 1]
    return df1, df2, {
        "batch": batch,
        "source": source,
        "d_img": int(img.shape[1]),
        "d_txt": int(txt.shape[1]),
    }


def split_subset(df1, df2):
    p = df1.shape[1] - 1
    return (
        np.asarray(df1[:, :p], float),
        np.asarray(df1[:, p], float),
        np.asarray(df2[:, :p], float),
        np.asarray(df2[:, p], float),
        p,
    )


def rf_domain_vimp(X0, X1, *, seed=SEED):
    X = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    p = X.shape[1]
    clf = RandomForestClassifier(
        n_estimators=150,
        max_depth=max(2, int(np.round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(np.round(np.sqrt(len(X)) // 2))),
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    Xtr, Xte, Wtr, Wte = train_test_split(X, W, test_size=0.25, random_state=seed, stratify=W)
    clf2 = RandomForestClassifier(
        n_estimators=120,
        max_depth=max(2, int(np.round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(np.round(np.sqrt(len(X)) // 2))),
        n_jobs=-1,
        random_state=seed + 1,
    )
    clf2.fit(Xtr, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return vimp, auc


def coord_mmd_vimp(X0, X1, *, max_n=MAX_MMD_N, seed=SEED, d_img=512):
    rng = np.random.default_rng(seed)
    A = _subsample_batch(X0, max_n, rng)
    B = _subsample_batch(X1, max_n, rng)
    mmd = MMD(compute_kernel="rbf")
    p = A.shape[1]
    uni = np.zeros(p)
    for j in range(p):
        uni[j], _ = mmd(A[:, [j]], B[:, [j]])
    full, _ = mmd(A, B)
    no_img, _ = mmd(A[:, d_img:], B[:, d_img:])
    no_txt, _ = mmd(A[:, :d_img], B[:, :d_img])
    block = {
        "mmd_full": float(full),
        "delta_drop_image": float(full - no_img),
        "delta_drop_text": float(full - no_txt),
    }
    return uni, block


def po_risk_vimp(X, Y, W, *, seed=SEED):
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    Yf = Y.astype(float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=80, max_depth=10, min_samples_leaf=3, random_state=seed + fold, n_jobs=-1
        )
        e = RandomForestClassifier(
            n_estimators=80, max_depth=10, min_samples_leaf=3, random_state=seed + 40 + fold, n_jobs=-1
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1 - ALPHA_CLIP)
    po = (Yf - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=200, max_depth=12, min_samples_leaf=3, random_state=seed + 7, n_jobs=-1
    )
    tau.fit(X, po)
    return tau.feature_importances_.astype(float), float(np.mean(tau.predict(X) ** 2))


def rank_modalities(vimp, d_img, d_txt, *, k=TOP_K):
    out, share = {}, {}
    blocks = {"image": (0, d_img), "text": (d_img, d_img + d_txt)}
    for m, (a, b) in blocks.items():
        vm = vimp[a:b]
        out[m] = [int(i) for i in np.argsort(-vm)[: min(k, len(vm))]]
        share[m] = float(vm.sum())
    tot = share["image"] + share["text"] + 1e-12
    share = {m: share[m] / tot for m in share}
    return out, share


def subsample_pair(X0, Y0, X1, Y1, *, n_max, seed):
    rng = np.random.default_rng(seed)
    n0 = min(len(X0), n_max // 2)
    n1 = min(len(X1), n_max // 2)
    i0 = rng.choice(len(X0), n0, replace=False)
    i1 = rng.choice(len(X1), n1, replace=False)
    return X0[i0], Y0[i0], X1[i1], Y1[i1]


def run_one(name: str) -> dict:
    print(f"\n======== {name} ========", flush=True)
    df1, df2, info = load_dataset(name)
    X0, Y0, X1, Y1, p = split_subset(df1, df2)
    d_img, d_txt = info["d_img"], info["d_txt"]
    assert p == d_img + d_txt
    print(
        f"n0={len(X0)} n1={len(X1)} p={p} batch={info['batch']}",
        flush=True,
    )
    X0s, Y0s, X1s, Y1s = subsample_pair(X0, Y0, X1, Y1, n_max=N_SUB_RF, seed=SEED)

    t0 = time.perf_counter()
    rf_vimp, rf_auc = rf_domain_vimp(X0s, X1s, seed=SEED)
    t_rf = time.perf_counter() - t0
    rf_fi, rf_share = rank_modalities(rf_vimp, d_img, d_txt)
    print(f"RF AUC={rf_auc:.3f} share={rf_share} ({t_rf:.1f}s)", flush=True)

    t0 = time.perf_counter()
    mmd_vimp, mmd_block = coord_mmd_vimp(X0, X1, max_n=MAX_MMD_N, seed=SEED, d_img=d_img)
    t_mmd = time.perf_counter() - t0
    mmd_fi, mmd_share = rank_modalities(mmd_vimp, d_img, d_txt)
    print(f"MMD share={mmd_share} block={mmd_block} ({t_mmd:.1f}s)", flush=True)

    X = np.vstack([X0s, X1s])
    Y = np.concatenate([Y0s, Y1s])
    W = np.concatenate([np.zeros(len(X0s)), np.ones(len(X1s))]).astype(int)
    t0 = time.perf_counter()
    po_vimp, po_risk = po_risk_vimp(X, Y, W, seed=SEED + 4)
    t_po = time.perf_counter() - t0
    po_fi, po_share = rank_modalities(po_vimp, d_img, d_txt)
    print(f"PO risk={po_risk:.3f} share={po_share} ({t_po:.1f}s)", flush=True)

    return {
        "dataset": name,
        "source": info["source"],
        "batch": info["batch"],
        "n0": int(len(X0)),
        "n1": int(len(X1)),
        "d_img": d_img,
        "d_txt": d_txt,
        "rf_domain_auc": rf_auc,
        "po_risk": po_risk,
        "time_seconds": {"rf": t_rf, "mmd": t_mmd, "po": t_po},
        "mmd_block": mmd_block,
        "boards": {
            "rf_domain": {"feature_indices": rf_fi, "modality_vimp_share": rf_share},
            "coord_mmd": {"feature_indices": mmd_fi, "modality_vimp_share": mmd_share},
            "po_risk": {"feature_indices": po_fi, "modality_vimp_share": po_share},
        },
    }


def write_artifacts(results: list[dict]):
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    summary = {r["dataset"]: r for r in results}
    (OUT / "multimodal_fsds_align_summary.json").write_text(json.dumps(summary, indent=2))

    proto = ["# Multimodal FSDS alignment · RF / coord-MMD / PO-risk\n"]
    for r in results:
        proto.append(f"\n# --- {r['dataset']} · {r['batch']} ---\n")
        for method, key in [
            ("rf_domain", "rf_domain"),
            ("coord_mmd", "coord_mmd"),
            ("po_risk", "po_risk"),
        ]:
            b = r["boards"][key]
            proto.append(f"feature_indices_{r['dataset']}_{method} = {repr(b['feature_indices'])}\n")
            proto.append(
                f"modality_vimp_share_{r['dataset']}_{method} = {repr(b['modality_vimp_share'])}\n"
            )
    (OUT / "multimodal_fsds_align_prototype.py").write_text("".join(proto))

    # LaTeX: cross-dataset mass share + per-dataset top indices (compact)
    lines = [
        "% Multimodal FSDS alignment · RF / coord-MMD / PO-risk\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Multimodal FSDS modality VIMP mass share across datasets. "
        "CS methods: RF Domain Classifier and coordinate-wise MMD; "
        "CD method: PO-risk RF. Ranking within modality uses "
        "$\\texttt{np.argsort(-VIMP)[:20]}$.}\n",
        "\\label{tab:mm-fsds-align-mass}\n\\small\n",
        "\\begin{tabular}{@{}l l cc cc cc@{}}\n\\toprule\n",
        "Dataset & Batch $W$ & "
        "\\multicolumn{2}{c}{RF Domain} & "
        "\\multicolumn{2}{c}{Coord-MMD} & "
        "\\multicolumn{2}{c}{PO-risk} \\\\\n",
        "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n",
        " &  & Img & Txt & Img & Txt & Img & Txt \\\\\n\\midrule\n",
    ]
    short = {
        "fashion_iq": "Fashion-IQ",
        "indiana_cxr": "Indiana CXR",
        "coco_outdoor_indoor": "COCO out/in (kw)",
        "coco_time_order": "COCO early/late",
        "microscopy_clip": "Microscopy",
    }
    for r in results:
        srf = r["boards"]["rf_domain"]["modality_vimp_share"]
        sm = r["boards"]["coord_mmd"]["modality_vimp_share"]
        sp = r["boards"]["po_risk"]["modality_vimp_share"]
        lines.append(
            f"{short[r['dataset']]} & {r['batch']} & "
            f"${srf['image']:.2f}$ & ${srf['text']:.2f}$ & "
            f"${sm['image']:.2f}$ & ${sm['text']:.2f}$ & "
            f"${sp['image']:.2f}$ & ${sp['text']:.2f}$ \\\\\n"
        )
    lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    lines.append("\\begin{table}[ht]\n\\centering\n")
    lines.append(
        "\\caption{Diagnostics for multimodal FSDS alignment "
        "(RF domain AUC; observed PO-risk).}\n"
    )
    lines.append("\\label{tab:mm-fsds-align-diag}\n\\small\n")
    lines.append("\\begin{tabular}{@{}l r r c c@{}}\n\\toprule\n")
    lines.append("Dataset & $n_0$ & $n_1$ & RF AUC & PO-risk \\\\\n\\midrule\n")
    for r in results:
        lines.append(
            f"{short[r['dataset']]} & {r['n0']} & {r['n1']} & "
            f"${r['rf_domain_auc']:.3f}$ & ${r['po_risk']:.3f}$ \\\\\n"
        )
    lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    # one compact indices table per dataset for RF (primary CS board)
    for r in results:
        fi = r["boards"]["rf_domain"]["feature_indices"]
        lines.append("\\begin{table}[ht]\n\\centering\n")
        lines.append(
            f"\\caption{{RF Domain ranked local indices on {short[r['dataset']]} "
            f"({r['batch']}). Domain AUC ${r['rf_domain_auc']:.3f}$.}}\n"
        )
        lines.append(f"\\label{{tab:mm-fsds-{r['dataset']}-rf}}\n\\small\n")
        lines.append("\\begin{tabular}{@{}l p{11.2cm}@{}}\n\\toprule\n")
        lines.append("Modality & Ranked local indices (high$\\to$low) \\\\\n\\midrule\n")
        for mod in ("image", "text"):
            idx = ",".join(str(i) for i in fi[mod])
            lines.append(f"\\texttt{{{mod}}} & $[{idx}]$ \\\\\n")
        lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    tex = "".join(lines)
    (DOCS / "Multimodal_FSDS_Align_tables_only.tex").write_text(tex)
    (OUT / "Multimodal_FSDS_Align_tables_only.tex").write_text(tex)

    md = [
        "# Multimodal FSDS alignment\n\n",
        "| Dataset | Batch W | RF img/txt | MMD img/txt | PO img/txt | RF AUC |\n",
        "|---------|---------|------------|-------------|------------|--------|\n",
    ]
    for r in results:
        srf = r["boards"]["rf_domain"]["modality_vimp_share"]
        sm = r["boards"]["coord_mmd"]["modality_vimp_share"]
        sp = r["boards"]["po_risk"]["modality_vimp_share"]
        md.append(
            f"| {r['dataset']} | {r['batch']} | "
            f"{srf['image']:.2f}/{srf['text']:.2f} | "
            f"{sm['image']:.2f}/{sm['text']:.2f} | "
            f"{sp['image']:.2f}/{sp['text']:.2f} | "
            f"{r['rf_domain_auc']:.3f} |\n"
        )
    (OUT / "README.md").write_text("".join(md))
    print(tex)
    print("".join(md))


def main():
    datasets = [
        "fashion_iq",
        "indiana_cxr",
        "coco_outdoor_indoor",
        "coco_time_order",
        "microscopy_clip",
    ]
    results = []
    for name in datasets:
        results.append(run_one(name))
    write_artifacts(results)
    print(f"wrote under {OUT}", flush=True)


if __name__ == "__main__":
    main()
