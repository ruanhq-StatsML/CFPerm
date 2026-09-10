#!/usr/bin/env python3
"""Image/Text multimodal FSDS · RF Domain Classifier ranked feature indices.

Uses FSDS-style RF Domain Classifier VIMP (not anomaly detection).
Modalities: image = img_feats columns, text = txt_feats columns.
Ground-truth modality blocks are the contiguous index ranges after concat.

Looks for:
  img_feats.npy, txt_feats.npy, labels.npy, df_metadata.csv
under CFPerm_Python / data/img_txt[/dataset] / env CFPERM_PYTHON_DIR.

  python3 scripts/run_img_txt_fsds.py
  python3 scripts/run_img_txt_fsds.py --data-dir data/img_txt/indiana_cxr
  python3 scripts/run_img_txt_fsds.py --all
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

try:
    import pandas as pd
except ImportError:
    import subprocess, sys

    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "pandas"])
    import pandas as pd

from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "img_txt_fsds"
SEED = 2026
TOP_K = 20
N_SUBSAMPLE = 10_000
ALPHA_CLIP = 0.01
NEED = ("img_feats.npy", "txt_feats.npy", "labels.npy", "df_metadata.csv")


def find_data_dirs() -> list[Path]:
    """Return dataset dirs that contain the four required files."""
    env = os.environ.get("CFPERM_PYTHON_DIR")
    bases: list[Path] = []
    if env:
        bases.append(Path(env))
    bases += [
        ROOT / "CFPerm_Python",
        ROOT / "cfperm_python",
        ROOT / "data" / "img_txt",
        ROOT / "data" / "CFPerm_Python",
        Path("/workspace/CFPerm_Python"),
        Path("/workspace/cfperm_python"),
    ]
    found: list[Path] = []
    seen: set[Path] = set()

    def add(p: Path):
        rp = p.resolve()
        if rp in seen:
            return
        if all((p / f).exists() for f in NEED):
            seen.add(rp)
            found.append(p)

    for base in bases:
        if not base or not base.exists():
            continue
        add(base)
        add(base / "data")
        add(base / "data" / "img_txt")
        add(base / "features")
        # nested dataset folders under data/img_txt/<name>/
        if base.name == "img_txt" or (base / "img_txt").is_dir():
            parent = base if base.name == "img_txt" else base / "img_txt"
            for sub in sorted(parent.iterdir()):
                if sub.is_dir():
                    add(sub)
    if not found:
        for p in ROOT.rglob("img_feats.npy"):
            add(p.parent)
    return found


def find_data_dir() -> Path:
    dirs = find_data_dirs()
    if not dirs:
        raise FileNotFoundError(
            "Need img_feats.npy, txt_feats.npy, labels.npy, df_metadata.csv. "
            "Put them under data/img_txt/<dataset>/ or set CFPERM_PYTHON_DIR."
        )
    return dirs[0]


def infer_batch_from_metadata(meta: pd.DataFrame, n: int) -> np.ndarray:
    """Build W∈{0,1} from metadata. Prefer explicit domain/batch/split columns."""
    cols = {c.lower(): c for c in meta.columns}
    # Prefer projection (IU Frontal/Lateral) and train/test split when present.
    preferred = (
        "projection",
        "split",
        "batch",
        "domain",
        "source",
        "shift",
        "group",
        "env",
        "year",
        "time",
    )
    for key in preferred:
        if key not in cols:
            continue
        c = cols[key]
        vals = meta[c].astype(str).to_numpy()
        uniq = pd.unique(vals)
        if len(uniq) < 2:
            continue
        if key == "year":
            try:
                y = pd.to_numeric(meta[c], errors="coerce").to_numpy()
                med = np.nanmedian(y)
                w = (y > med).astype(int)
                if w.sum() > 0 and w.sum() < n:
                    print(f"batch from metadata column '{c}' (year median)", flush=True)
                    return w
            except Exception:
                pass
        # stable binary: first unique → 0, rest → 1
        # for projection prefer Frontal=0 / Lateral=1; for split train=0 / test=1
        if key == "projection":
            w = np.array([0 if v.lower().startswith("front") else 1 for v in vals], dtype=int)
        elif key == "split":
            w = np.array([0 if v.lower().startswith("train") else 1 for v in vals], dtype=int)
        else:
            w = (vals != uniq[0]).astype(int)
        if w.sum() > 0 and w.sum() < n:
            print(f"batch from metadata column '{c}': {uniq.tolist()[:6]}", flush=True)
            return w
    print("batch fallback: first half vs second half of rows", flush=True)
    w = np.zeros(n, int)
    w[n // 2 :] = 1
    return w


def rf_domain_vimp(X, W, *, seed: int = SEED):
    """FSDS covariate-shift board: RF Domain Classifier impurity VIMP."""
    clf = RandomForestClassifier(
        n_estimators=200, max_depth=12, min_samples_leaf=3, random_state=seed, n_jobs=-1
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    Xtr, Xte, Wtr, Wte = train_test_split(X, W, test_size=0.25, random_state=seed, stratify=W)
    clf2 = RandomForestClassifier(
        n_estimators=120, max_depth=12, min_samples_leaf=3, random_state=seed + 1, n_jobs=-1
    )
    clf2.fit(Xtr, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return vimp, auc


def concept_drift_vimp(X, Y, W, *, seed: int = SEED):
    """FSDS concept-drift board: PO-risk RF VIMP + within-batch Y~X gap."""
    n, p = X.shape
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
    vimp_po = tau.feature_importances_.astype(float)
    po_risk = float(np.mean(tau.predict(X) ** 2))

    y0 = RandomForestRegressor(
        n_estimators=120, max_depth=10, min_samples_leaf=3, random_state=seed + 11, n_jobs=-1
    )
    y1 = RandomForestRegressor(
        n_estimators=120, max_depth=10, min_samples_leaf=3, random_state=seed + 13, n_jobs=-1
    )
    y0.fit(X[W == 0], Yf[W == 0])
    y1.fit(X[W == 1], Yf[W == 1])
    gap = np.abs(y1.feature_importances_ - y0.feature_importances_)
    vimp = 0.7 * (vimp_po / (vimp_po.sum() + 1e-12)) + 0.3 * (gap / (gap.sum() + 1e-12))
    return vimp, po_risk


def rank_modalities(vimp, d_img, d_txt, *, k: int = TOP_K):
    """Ground-truth modality blocks: image [0,d_img), text [d_img, d_img+d_txt)."""
    blocks = {"image": (0, d_img), "text": (d_img, d_img + d_txt)}
    out, detail = {}, {}
    for m, (a, b) in blocks.items():
        vm = vimp[a:b]
        order = [int(i) for i in np.argsort(-vm)[: min(k, len(vm))]]
        out[m] = order  # local indices, high→low
        detail[m] = [
            {
                "rank": r + 1,
                "index": idx,
                "global_index": int(a + idx),
                "vimp": float(vm[idx]),
            }
            for r, idx in enumerate(order)
        ]
    mass = {
        "image": float(vimp[:d_img].sum()),
        "text": float(vimp[d_img:].sum()),
    }
    total = mass["image"] + mass["text"] + 1e-12
    mass_share = {m: mass[m] / total for m in mass}
    return out, detail, mass_share, blocks


def subsample_balanced(X, Y, W, *, n_max=N_SUBSAMPLE, seed=SEED):
    rng = np.random.default_rng(seed)
    i0, i1 = np.where(W == 0)[0], np.where(W == 1)[0]
    n_each = min(len(i0), len(i1), n_max // 2)
    take = np.concatenate(
        [rng.choice(i0, n_each, replace=False), rng.choice(i1, n_each, replace=False)]
    )
    rng.shuffle(take)
    return X[take], Y[take], W[take]


def run_one(data_dir: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n==== data_dir={data_dir} ====", flush=True)

    img = np.load(data_dir / "img_feats.npy")
    txt = np.load(data_dir / "txt_feats.npy")
    labels = np.load(data_dir / "labels.npy")
    meta = pd.read_csv(data_dir / "df_metadata.csv")
    print(
        f"img={img.shape} txt={txt.shape} labels={labels.shape} meta={meta.shape} cols={list(meta.columns)[:12]}",
        flush=True,
    )
    n = min(len(img), len(txt), len(labels), len(meta))
    img, txt, labels, meta = img[:n], txt[:n], labels[:n], meta.iloc[:n].reset_index(drop=True)
    d_img, d_txt = int(img.shape[1]), int(txt.shape[1])
    X = np.hstack([img.astype(float), txt.astype(float)])
    Y = labels.reshape(-1)
    W = infer_batch_from_metadata(meta, n)

    X, Y, W = subsample_balanced(X, Y, W)
    print(
        f"subsample n={len(W)} p={X.shape[1]} d_img={d_img} d_txt={d_txt} "
        f"n0={(W == 0).sum()} n1={(W == 1).sum()}",
        flush=True,
    )

    modality_index_ground_truth = {
        "image": list(range(d_img)),
        "text": list(range(d_txt)),
        "image_global": list(range(0, d_img)),
        "text_global": list(range(d_img, d_img + d_txt)),
    }

    print("CS · RF Domain Classifier…", flush=True)
    cs_vimp, cs_auc = rf_domain_vimp(X, W, seed=SEED)
    cs_fi, cs_det, cs_mass, blocks = rank_modalities(cs_vimp, d_img, d_txt, k=TOP_K)
    print("CS feature_indices:", cs_fi, flush=True)
    print("CS modality mass share:", cs_mass, flush=True)

    print("CD · PO-risk RF…", flush=True)
    cd_vimp, po_risk = concept_drift_vimp(X, Y, W, seed=SEED)
    cd_fi, cd_det, cd_mass, _ = rank_modalities(cd_vimp, d_img, d_txt, k=TOP_K)
    print("CD feature_indices:", cd_fi, flush=True)
    print("CD modality mass share:", cd_mass, flush=True)

    payload = {
        "method": "FSDS · RF Domain Classifier (CS) + PO-risk RF (CD)",
        "dataset": data_dir.name,
        "data_dir": str(data_dir),
        "n": int(len(W)),
        "d_img": d_img,
        "d_txt": d_txt,
        "top_k": TOP_K,
        "rf_domain_auc": cs_auc,
        "po_risk": po_risk,
        "modality_index_ground_truth": {
            "image_global": modality_index_ground_truth["image_global"],
            "text_global": modality_index_ground_truth["text_global"],
            "note": "concat X=[img|txt]; local indices are within-modality",
        },
        "covariate_shift": {
            "feature_indices": cs_fi,
            "modality_vimp_share": cs_mass,
            "detail": cs_det,
        },
        "concept_drift": {
            "feature_indices": cd_fi,
            "modality_vimp_share": cd_mass,
            "detail": cd_det,
        },
    }
    stem = data_dir.name
    dest = out_dir / f"{stem}_fsds_feature_indices.json"
    dest.write_text(json.dumps(payload, indent=2))
    print(f"wrote {dest}", flush=True)

    proto = out_dir / f"{stem}_fsds_feature_indices_prototype.py"
    proto.write_text(
        f"# Image/Text FSDS · {stem} · ranked feature indices (order = importance)\n"
        f"feature_indices_covariate_shift = {repr(cs_fi)}\n"
        f"feature_indices_concept_drift = {repr(cd_fi)}\n"
        f"modality_vimp_share_cs = {repr(cs_mass)}\n"
        f"modality_vimp_share_cd = {repr(cd_mass)}\n"
        f"rf_domain_auc = {cs_auc!r}\n"
        f"po_risk = {po_risk!r}\n"
    )
    print(f"wrote {proto}", flush=True)
    print(proto.read_text())
    return payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=None, help="Single dataset directory")
    ap.add_argument("--all", action="store_true", help="Run all discovered img/txt datasets")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    if args.data_dir:
        dirs = [args.data_dir]
    else:
        dirs = find_data_dirs()
        # Prefer nested data/img_txt/<dataset>/ folders when present.
        nested = [p for p in dirs if p.parent.name == "img_txt"]
        dirs = nested or dirs

    if not dirs:
        raise FileNotFoundError("No img/txt datasets found")

    summary = {}
    for d in dirs:
        payload = run_one(d, out)
        summary[d.name] = {
            "feature_indices_covariate_shift": payload["covariate_shift"]["feature_indices"],
            "feature_indices_concept_drift": payload["concept_drift"]["feature_indices"],
            "modality_vimp_share_cs": payload["covariate_shift"]["modality_vimp_share"],
            "modality_vimp_share_cd": payload["concept_drift"]["modality_vimp_share"],
            "rf_domain_auc": payload["rf_domain_auc"],
            "po_risk": payload["po_risk"],
        }

    board = out / "img_txt_fsds_both_datasets_board.py"
    lines = ["# Image/Text FSDS · both datasets · ordered argsort(-VIMP)[:k]\n"]
    for name, s in summary.items():
        safe = name.replace("-", "_")
        lines.append(f"# --- {name} ---\n")
        lines.append(f"feature_indices_covariate_shift_{safe} = {repr(s['feature_indices_covariate_shift'])}\n")
        lines.append(f"feature_indices_concept_drift_{safe} = {repr(s['feature_indices_concept_drift'])}\n")
        lines.append(f"modality_vimp_share_cs_{safe} = {repr(s['modality_vimp_share_cs'])}\n")
        lines.append(f"modality_vimp_share_cd_{safe} = {repr(s['modality_vimp_share_cd'])}\n")
    board.write_text("".join(lines))
    (out / "img_txt_fsds_both_datasets_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {board}", flush=True)
    print(board.read_text())


if __name__ == "__main__":
    main()
