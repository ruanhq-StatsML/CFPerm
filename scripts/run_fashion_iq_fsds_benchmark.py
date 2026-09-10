#!/usr/bin/env python3
"""Fashion-IQ CLIP embedding · FSDS benchmark_feature_selection subset.

Data layout (FSDS / benchmark_whole_feature_selection convention):
  clip_embedding.npy  shape (n, 1025) = [image(512) | text(512) | label]
  df1.npy / df2.npy   train vs test batches (last column = Y)

Methods (from FSDS bakeoff):
  1) RF Domain Classifier VIMP
  2) LOCO-MMD VIMP
  3) PO-risk RF VIMP  (pseudo-outcome risk / PermuCATE-style scoring)

  python3 scripts/run_fashion_iq_fsds_benchmark.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor" / "fsds"))

from VIMP_mmd_benchmark import (  # noqa: E402
    MMD,
    _subsample_batch,
    compute_loco_mmd_batches,
)

DATA = ROOT / "data" / "img_txt" / "fashion_iq"
OUT = ROOT / "results" / "fashion_iq_fsds"
DOCS = ROOT / "docs" / "method"

SEED = 2026
TOP_K = 20
D_IMG = 512
D_TXT = 512
MAX_MMD_N = 100
ALPHA_CLIP = 0.01
N_SUB_RF = 4000  # cap rows for RF / PO paths


def split_subset(df1: np.ndarray, df2: np.ndarray):
    """FSDS split_subset: last column is label."""
    p = df1.shape[1] - 1
    return (
        np.asarray(df1[:, :p], dtype=float),
        np.asarray(df1[:, p], dtype=float),
        np.asarray(df2[:, :p], dtype=float),
        np.asarray(df2[:, p], dtype=float),
        p,
    )


def compute_mmd_feature_vimp(df1_X, df2_X, *, max_n: int = MAX_MMD_N, seed: int = SEED):
    """MMD feature scores for high-d CLIP.

    Full leave-one-coordinate-out MMD is numerically flat when p≈1024
    (ΔMMD ≈ 0). Use coordinate-wise 1-D MMD as the feature board, and
    report modality-block LOCO-MMD as a sanity check.
    """
    rng = np.random.default_rng(seed)
    X0 = _subsample_batch(np.asarray(df1_X, dtype=float), max_n, rng)
    X1 = _subsample_batch(np.asarray(df2_X, dtype=float), max_n, rng)
    mmd = MMD(compute_kernel="rbf")
    p = X0.shape[1]
    uni = np.zeros(p, dtype=float)
    for j in range(p):
        uni[j], _ = mmd(X0[:, [j]], X1[:, [j]])
        if (j + 1) % 256 == 0:
            print(f"  coord-MMD {j+1}/{p}", flush=True)

    # modality-block LOCO (image vs text)
    mmd_full, _ = mmd(X0, X1)
    mmd_no_img, _ = mmd(X0[:, D_IMG:], X1[:, D_IMG:])
    mmd_no_txt, _ = mmd(X0[:, :D_IMG], X1[:, :D_IMG])
    block = {
        "mmd_full": float(mmd_full),
        "delta_drop_image": float(mmd_full - mmd_no_img),
        "delta_drop_text": float(mmd_full - mmd_no_txt),
    }

    # keep classic LOCO call on a tiny probe to document flatness
    loco_probe, _ = compute_loco_mmd_batches(
        X0[:, :32], X1[:, :32], max_n=max_n, seed=seed
    )
    block["loco_probe_absmax_first32"] = float(np.max(np.abs(loco_probe)))
    return uni, block


def compute_rf_domain_vimp(df1_X, df2_X, *, seed: int = SEED):
    """FSDS compute_rf_domain_vimp."""
    X = np.vstack([df1_X, df2_X])
    W = np.concatenate([np.zeros(len(df1_X)), np.ones(len(df2_X))]).astype(int)
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
    Xtr, Xte, Wtr, Wte = train_test_split(
        X, W, test_size=0.25, random_state=seed, stratify=W
    )
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


def po_risk_vimp(X, Y, W, *, seed: int = SEED):
    """PO-risk path: cross-fit m,e → pseudo-outcome → RF VIMP on τ(X)."""
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
    vimp = tau.feature_importances_.astype(float)
    po_risk = float(np.mean(tau.predict(X) ** 2))
    return vimp, po_risk


def metrics_from_scores(scores, feature_ind, p):
    """Top-|S| selection vs ground-truth index set (FSDS metrics_from_scores)."""
    scores = np.asarray(scores, dtype=float)
    feature_ind = np.asarray(feature_ind, dtype=int).ravel()
    k = len(feature_ind)
    selected = np.argsort(-scores)[:k]
    true = set(feature_ind.tolist())
    sel = set(selected.tolist())
    tp = len(sel & true)
    fp = len(sel - true)
    fn = len(true - sel)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    fdr = fp / (tp + fp) if (tp + fp) else 0.0
    out = {
        "TP": int(tp),
        "FP": int(fp),
        "FN": int(fn),
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "F1": round(f1, 4),
        "FDR": round(fdr, 4),
        "selected_top_k": [int(i) for i in selected[: min(20, k)]],
    }
    if 0 < len(true) < p:
        membership = np.zeros(p, dtype=int)
        membership[feature_ind] = 1
        try:
            out["selection_AUC"] = round(float(roc_auc_score(membership, scores)), 4)
        except ValueError:
            out["selection_AUC"] = None
    return out


def rank_modalities(vimp, *, k: int = TOP_K):
    blocks = {"image": (0, D_IMG), "text": (D_IMG, D_IMG + D_TXT)}
    out, detail, mass = {}, {}, {}
    for m, (a, b) in blocks.items():
        vm = vimp[a:b]
        order = [int(i) for i in np.argsort(-vm)[: min(k, len(vm))]]
        out[m] = order
        detail[m] = [
            {"rank": r + 1, "index": idx, "global_index": int(a + idx), "vimp": float(vm[idx])}
            for r, idx in enumerate(order)
        ]
        mass[m] = float(vm.sum())
    total = mass["image"] + mass["text"] + 1e-12
    share = {m: mass[m] / total for m in mass}
    return out, detail, share


def subsample_pair(df1_X, df1_Y, df2_X, df2_Y, *, n_max: int, seed: int):
    rng = np.random.default_rng(seed)
    n1 = min(len(df1_X), n_max // 2)
    n2 = min(len(df2_X), n_max // 2)
    i1 = rng.choice(len(df1_X), n1, replace=False)
    i2 = rng.choice(len(df2_X), n2, replace=False)
    return df1_X[i1], df1_Y[i1], df2_X[i2], df2_Y[i2]


def benchmark_feature_selection(df1, df2, *, seed: int = SEED, max_mmd_n: int = MAX_MMD_N):
    """Subset of FSDS benchmark_whole_feature_selection: RF / MMD / PO-risk."""
    df1_X, df1_Y, df2_X, df2_Y, p = split_subset(df1, df2)
    assert p == D_IMG + D_TXT, f"expected p={D_IMG+D_TXT}, got {p}"

    # Ground-truth modality blocks for selection metrics (image vs text attribution)
    gt_image = np.arange(0, D_IMG)
    gt_text = np.arange(D_IMG, D_IMG + D_TXT)

    df1_Xs, df1_Ys, df2_Xs, df2_Ys = subsample_pair(
        df1_X, df1_Y, df2_X, df2_Y, n_max=N_SUB_RF, seed=seed
    )
    print(
        f"RF/PO subsample n1={len(df1_Xs)} n2={len(df2_Xs)} p={p}",
        flush=True,
    )

    time_dict, scores, metrics = {}, {}, {}

    t0 = time.perf_counter()
    print("RF Domain Classifier…", flush=True)
    rf_vimp, rf_auc = compute_rf_domain_vimp(df1_Xs, df2_Xs, seed=seed)
    time_dict["rf_domain"] = float(time.perf_counter() - t0)
    scores["rf_domain"] = rf_vimp
    metrics["rf_domain"] = {
        "vs_image_block": metrics_from_scores(rf_vimp, gt_image, p),
        "vs_text_block": metrics_from_scores(rf_vimp, gt_text, p),
        "domain_auc": rf_auc,
    }
    print(f"  AUC={rf_auc:.4f} time={time_dict['rf_domain']:.1f}s", flush=True)

    t0 = time.perf_counter()
    print(f"MMD feature VIMP (coord-wise + block LOCO; max_n={max_mmd_n})…", flush=True)
    mmd_vimp, mmd_block = compute_mmd_feature_vimp(
        df1_X, df2_X, max_n=max_mmd_n, seed=seed
    )
    time_dict["loco_mmd"] = float(time.perf_counter() - t0)
    scores["loco_mmd"] = mmd_vimp
    metrics["loco_mmd"] = {
        "vs_image_block": metrics_from_scores(mmd_vimp, gt_image, p),
        "vs_text_block": metrics_from_scores(mmd_vimp, gt_text, p),
        "block_loco": mmd_block,
        "note": (
            "Coordinate-wise 1-D MMD ranking; full leave-one-out MMD is flat "
            "at p=1024 (probe absmax recorded in block_loco)."
        ),
    }
    print(f"  block={mmd_block} time={time_dict['loco_mmd']:.1f}s", flush=True)

    X = np.vstack([df1_Xs, df2_Xs])
    Y = np.concatenate([df1_Ys, df2_Ys])
    W = np.concatenate([np.zeros(len(df1_Xs)), np.ones(len(df2_Xs))]).astype(int)
    t0 = time.perf_counter()
    print("PO-risk RF VIMP…", flush=True)
    po_vimp, po_risk = po_risk_vimp(X, Y, W, seed=seed + 4)
    time_dict["po_risk"] = float(time.perf_counter() - t0)
    scores["po_risk"] = po_vimp
    metrics["po_risk"] = {
        "vs_image_block": metrics_from_scores(po_vimp, gt_image, p),
        "vs_text_block": metrics_from_scores(po_vimp, gt_text, p),
        "po_risk": po_risk,
    }
    print(f"  po_risk={po_risk:.4f} time={time_dict['po_risk']:.1f}s", flush=True)

    boards = {}
    for name, vimp in scores.items():
        fi, det, share = rank_modalities(vimp, k=TOP_K)
        boards[name] = {
            "feature_indices": fi,
            "modality_vimp_share": share,
            "detail": det,
        }
        print(f"  {name} share={share}", flush=True)
        print(f"  {name} feature_indices={fi}", flush=True)

    return {
        "dataset": "fashion_iq_clip",
        "layout": "X=[img(512)|txt(512)], Y=last column; W=train vs test",
        "n1": int(len(df1_X)),
        "n2": int(len(df2_X)),
        "p": int(p),
        "d_img": D_IMG,
        "d_txt": D_TXT,
        "top_k": TOP_K,
        "max_mmd_n": max_mmd_n,
        "time_seconds": time_dict,
        "metrics": metrics,
        "boards": boards,
        "methods": ["rf_domain", "loco_mmd", "po_risk"],
    }


def write_artifacts(payload: dict):
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    (OUT / "fashion_iq_fsds_benchmark.json").write_text(json.dumps(payload, indent=2))

    boards = payload["boards"]
    proto = [
        "# Fashion-IQ CLIP · FSDS benchmark_feature_selection (RF / MMD / PO-risk)\n",
        "# Ranking = np.argsort(-VIMP)[:20] within modality; X=[img|txt], Y=last col\n",
    ]
    for method in payload["methods"]:
        b = boards[method]
        proto.append(f"feature_indices_{method} = {repr(b['feature_indices'])}\n")
        proto.append(f"modality_vimp_share_{method} = {repr(b['modality_vimp_share'])}\n")
    (OUT / "fashion_iq_fsds_benchmark_prototype.py").write_text("".join(proto))

    # LaTeX tables
    lines = [
        "% Fashion-IQ CLIP FSDS · RF / LOCO-MMD / PO-risk ranked indices\n",
        "% Requires booktabs.\n\n",
    ]
    n1, n2 = payload["n1"], payload["n2"]
    auc = payload["metrics"]["rf_domain"]["domain_auc"]
    por = payload["metrics"]["po_risk"]["po_risk"]
    mmd_n = payload["max_mmd_n"]
    s_rf = boards["rf_domain"]["modality_vimp_share"]
    s_mmd = boards["loco_mmd"]["modality_vimp_share"]
    s_po = boards["po_risk"]["modality_vimp_share"]
    captions = {
        "rf_domain": (
            "rf-domain",
            (
                f"RF Domain Classifier on Fashion-IQ CLIP "
                f"($n_1={n1}$, $n_2={n2}$, "
                f"$d_{{\\mathrm{{img}}}}=d_{{\\mathrm{{txt}}}}=512$). "
                f"Batch: train vs test. Domain AUC ${auc:.3f}$; "
                f"modality share image ${s_rf['image']:.3f}$ / text ${s_rf['text']:.3f}$."
            ),
        ),
        "loco_mmd": (
            "coord-mmd",
            (
                f"Coordinate-wise 1-D MMD ranking on Fashion-IQ CLIP "
                f"(max\\_n$={mmd_n}$; full leave-one-out MMD is numerically flat at $p=1024$). "
                f"Modality share image ${s_mmd['image']:.3f}$ / text ${s_mmd['text']:.3f}$."
            ),
        ),
        "po_risk": (
            "po-risk",
            (
                f"PO-risk RF VIMP on Fashion-IQ CLIP. Observed PO-risk ${por:.3f}$; "
                f"modality share image ${s_po['image']:.3f}$ / text ${s_po['text']:.3f}$."
            ),
        ),
    }
    for method, (lab, cap) in captions.items():
        fi = boards[method]["feature_indices"]
        lines.append("\\begin{table}[ht]\n\\centering\n")
        lines.append(f"\\caption{{{cap}}}\n")
        lines.append(f"\\label{{tab:fashion-iq-{lab}}}\n\\small\n")
        lines.append("\\begin{tabular}{@{}l p{11.5cm}@{}}\n\\toprule\n")
        lines.append("Modality & Ranked local indices (high$\\to$low) \\\\\n\\midrule\n")
        for mod in ("image", "text"):
            idx = ",".join(str(i) for i in fi[mod])
            lines.append(f"\\texttt{{{mod}}} & $[{idx}]$ \\\\\n")
        lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n\n")

    # summary mass table
    lines.append("\\begin{table}[ht]\n\\centering\n")
    lines.append(
        "\\caption{Fashion-IQ CLIP modality VIMP mass share across "
        "RF Domain / LOCO-MMD / PO-risk.}\n"
    )
    lines.append("\\label{tab:fashion-iq-modality-mass}\n\\small\n")
    lines.append("\\begin{tabular}{@{}l cc@{}}\n\\toprule\n")
    lines.append("Method & Image share & Text share \\\\\n\\midrule\n")
    for method, label in [
        ("rf_domain", "RF Domain Classifier"),
        ("loco_mmd", "Coord-wise MMD"),
        ("po_risk", "PO-risk RF"),
    ]:
        s = boards[method]["modality_vimp_share"]
        lines.append(f"{label} & ${s['image']:.3f}$ & ${s['text']:.3f}$ \\\\\n")
    lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    # modality-block LOCO-MMD sanity (Δ when dropping image vs text block)
    block = payload["metrics"]["loco_mmd"].get("block_loco", {})
    if block:
        lines.append("\n\\begin{table}[ht]\n\\centering\n")
        lines.append(
            "\\caption{Fashion-IQ CLIP modality-block LOCO-MMD "
            "(drop image block vs drop text block; subsample max\\_n).}\n"
        )
        lines.append("\\label{tab:fashion-iq-block-loco-mmd}\n\\small\n")
        lines.append("\\begin{tabular}{@{}l c@{}}\n\\toprule\n")
        lines.append("Quantity & Value \\\\\n\\midrule\n")
        lines.append(f"MMD full & ${block.get('mmd_full', float('nan')):.6f}$ \\\\\n")
        lines.append(
            f"$\\Delta$ drop image & ${block.get('delta_drop_image', float('nan')):.6f}$ \\\\\n"
        )
        lines.append(
            f"$\\Delta$ drop text & ${block.get('delta_drop_text', float('nan')):.6f}$ \\\\\n"
        )
        lines.append("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    tex = "".join(lines)
    (DOCS / "FashionIQ_FSDS_RF_MMD_POrisk_tables_only.tex").write_text(tex)
    (OUT / "FashionIQ_FSDS_RF_MMD_POrisk_tables_only.tex").write_text(tex)
    print(tex)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df1 = np.load(DATA / "df1.npy")
    df2 = np.load(DATA / "df2.npy")
    print(f"loaded df1={df1.shape} df2={df2.shape} from {DATA}", flush=True)
    payload = benchmark_feature_selection(df1, df2, seed=SEED, max_mmd_n=MAX_MMD_N)
    write_artifacts(payload)
    print(f"wrote under {OUT}", flush=True)


if __name__ == "__main__":
    main()
