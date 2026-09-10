#!/usr/bin/env python3
"""Hierarchical multimodal attribution board · modality → bbox subgroups.

Board levels (evidence only; business reads the board):
  L1  modality mass share on blocks in X=[image | text | bbox]
  L2  within-block feature ranking  argsort(-VIMP)[:k]
  L2b for bbox block only: subgroup mass
        geo | category_hist | top_boxes

Also runs a late-batch bbox inject recovery check (GT=bbox).

  python3 scripts/run_bbox_attribution_board.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor" / "fsds"))
from VIMP_mmd_benchmark import MMD, _subsample_batch  # noqa: E402

DATA = ROOT / "data" / "img_txt" / "coco_time_order"
OUT = ROOT / "results" / "bbox_attribution"
DOCS = ROOT / "docs" / "method"

SEED = 2026
TOP_K = 20
MAX_MMD_N = 100
N_SUB = 4000
INJECT_ALPHA = 0.75


def load_blocks():
    img = np.load(DATA / "img_feats.npy").astype(float)
    txt = np.load(DATA / "txt_feats.npy").astype(float)
    bbox = np.load(DATA / "bbox_feats.npy").astype(float)
    meta = pd.read_csv(DATA / "df_metadata.csv")
    schema = json.loads((DATA / "bbox_feature_schema.json").read_text())
    w = meta["batch"].astype(int).to_numpy()
    blocks = {
        "image": (0, img.shape[1]),
        "text": (img.shape[1], img.shape[1] + txt.shape[1]),
        "bbox": (
            img.shape[1] + txt.shape[1],
            img.shape[1] + txt.shape[1] + bbox.shape[1],
        ),
    }
    X = np.hstack([img, txt, bbox])
    return X[w == 0], X[w == 1], blocks, schema, meta


def rf_domain_vimp(X0, X1, *, seed=SEED):
    X = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    p = X.shape[1]
    clf = RandomForestClassifier(
        n_estimators=200,
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


def coord_mmd_vimp(X0, X1, *, max_n=MAX_MMD_N, seed=SEED):
    rng = np.random.default_rng(seed)
    A = _subsample_batch(X0, max_n, rng)
    B = _subsample_batch(X1, max_n, rng)
    mmd = MMD(compute_kernel="rbf")
    uni = np.zeros(A.shape[1])
    for j in range(A.shape[1]):
        uni[j], _ = mmd(A[:, [j]], B[:, [j]])
    return uni


def subsample(X0, X1, *, n_max=N_SUB, seed=SEED):
    rng = np.random.default_rng(seed)
    n0 = min(len(X0), n_max // 2)
    n1 = min(len(X1), n_max // 2)
    return X0[rng.choice(len(X0), n0, replace=False)], X1[rng.choice(len(X1), n1, replace=False)]


def modality_board(vimp, blocks, *, k=TOP_K):
    mass = {m: float(vimp[a:b].sum()) for m, (a, b) in blocks.items()}
    tot = sum(mass.values()) + 1e-12
    share = {m: mass[m] / tot for m in mass}
    ranks = {}
    for m, (a, b) in blocks.items():
        vm = vimp[a:b]
        order = [int(i) for i in np.argsort(-vm)[: min(k, len(vm))]]
        ranks[m] = [
            {"rank": r + 1, "local_index": idx, "global_index": int(a + idx), "vimp": float(vm[idx])}
            for r, idx in enumerate(order)
        ]
    return share, ranks


def bbox_subgroup_board(vimp, blocks, schema):
    """L2b: geo / category_hist / top_boxes mass inside bbox block."""
    a, b = blocks["bbox"]
    vb = vimp[a:b]
    g0, g1 = schema["geo_slice"]
    c0, c1 = schema["category_hist_slice"]
    t0, t1 = schema["top_boxes_slice"]
    parts = {
        "geo": float(vb[g0:g1].sum()),
        "category_hist": float(vb[c0:c1].sum()),
        "top_boxes": float(vb[t0:t1].sum()),
    }
    tot = sum(parts.values()) + 1e-12
    share = {k: parts[k] / tot for k in parts}
    # named geo ranking
    geo_names = schema["geo_names"]
    geo_v = vb[g0:g1]
    geo_rank = [
        {"rank": r + 1, "name": geo_names[i], "local_index": int(i), "vimp": float(geo_v[i])}
        for r, i in enumerate(np.argsort(-geo_v))
    ]
    # category hist top
    cat_ids = schema["category_ids"]
    cat_v = vb[c0:c1]
    cat_rank = [
        {
            "rank": r + 1,
            "category_id": int(cat_ids[i]),
            "local_index": int(i),
            "vimp": float(cat_v[i]),
        }
        for r, i in enumerate(np.argsort(-cat_v)[:TOP_K])
    ]
    return share, geo_rank, cat_rank


def metrics_vs_bbox(scores, blocks):
    a, b = blocks["bbox"]
    pos = set(range(a, b))
    membership = np.zeros(len(scores), int)
    membership[a:b] = 1
    auc = float(roc_auc_score(membership, scores))
    top20 = np.argsort(-scores)[:TOP_K]
    p20 = float(np.mean([i in pos for i in top20]))
    mass = float(scores[a:b].sum()) / (float(scores.sum()) + 1e-12)
    return {"selection_AUC": round(auc, 4), "P@20": round(p20, 4), "mass_on_bbox": round(mass, 4)}


def inject_bbox(X0, X1, blocks, *, alpha=INJECT_ALPHA, seed=SEED):
    rng = np.random.default_rng(seed)
    X1 = X1.copy()
    a, b = blocks["bbox"]
    d = b - a
    u = rng.normal(size=d)
    u /= np.linalg.norm(u) + 1e-12
    X1[:, a:b] += alpha * u
    return X0, X1


def run_one(tag, X0, X1, blocks, schema):
    X0s, X1s = subsample(X0, X1)
    t0 = time.perf_counter()
    rf_vimp, rf_auc = rf_domain_vimp(X0s, X1s)
    t_rf = time.perf_counter() - t0
    rf_share, rf_ranks = modality_board(rf_vimp, blocks)
    rf_bbox_share, rf_geo, rf_cat = bbox_subgroup_board(rf_vimp, blocks, schema)

    t0 = time.perf_counter()
    mmd_vimp = coord_mmd_vimp(X0, X1)
    t_mmd = time.perf_counter() - t0
    mmd_share, mmd_ranks = modality_board(mmd_vimp, blocks)
    mmd_bbox_share, mmd_geo, mmd_cat = bbox_subgroup_board(mmd_vimp, blocks, schema)

    out = {
        "scenario": tag,
        "n0": int(len(X0)),
        "n1": int(len(X1)),
        "block_dims": {m: b - a for m, (a, b) in blocks.items()},
        "rf_domain": {
            "domain_auc": round(rf_auc, 4),
            "time_s": round(t_rf, 2),
            "modality_mass_share": rf_share,
            "feature_ranking": {
                m: [r["local_index"] for r in rf_ranks[m]] for m in rf_ranks
            },
            "bbox_subgroup_mass_share": rf_bbox_share,
            "bbox_geo_ranking": rf_geo,
            "bbox_category_ranking_top20": rf_cat,
        },
        "coord_mmd": {
            "time_s": round(t_mmd, 2),
            "modality_mass_share": mmd_share,
            "feature_ranking": {
                m: [r["local_index"] for r in mmd_ranks[m]] for m in mmd_ranks
            },
            "bbox_subgroup_mass_share": mmd_bbox_share,
            "bbox_geo_ranking": mmd_geo,
            "bbox_category_ranking_top20": mmd_cat,
        },
    }
    if "inject" in tag:
        out["rf_domain"]["bbox_recovery"] = metrics_vs_bbox(rf_vimp, blocks)
        out["coord_mmd"]["bbox_recovery"] = metrics_vs_bbox(mmd_vimp, blocks)
    return out


def write_artifacts(baseline, inject):
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {"baseline_time_order": baseline, "bbox_inject": inject}
    (OUT / "bbox_attribution_board.json").write_text(json.dumps(payload, indent=2))

    # prototype board (numbers only)
    proto = [
        "# Hierarchical attribution board · image / text / bbox\n",
        "# L1 = modality mass share; L2 = within-block ranking; L2b = bbox subgroups\n",
        f"modality_mass_share_rf = {repr(baseline['rf_domain']['modality_mass_share'])}\n",
        f"feature_ranking_rf = {repr(baseline['rf_domain']['feature_ranking'])}\n",
        f"bbox_subgroup_mass_share_rf = {repr(baseline['rf_domain']['bbox_subgroup_mass_share'])}\n",
        f"bbox_geo_ranking_rf = {repr([x['name'] for x in baseline['rf_domain']['bbox_geo_ranking']])}\n",
        f"modality_mass_share_mmd = {repr(baseline['coord_mmd']['modality_mass_share'])}\n",
        f"bbox_inject_recovery_rf = {repr(inject['rf_domain']['bbox_recovery'])}\n",
        f"bbox_inject_recovery_mmd = {repr(inject['coord_mmd']['bbox_recovery'])}\n",
    ]
    (OUT / "bbox_attribution_board_prototype.py").write_text("".join(proto))

    br = baseline["rf_domain"]
    bm = baseline["coord_mmd"]
    ir = inject["rf_domain"]
    im = inject["coord_mmd"]

    md = [
        "# Bounding-box hierarchical attribution board\n\n",
        "Evidence board only. Business decides what the rankings mean.\n\n",
        "## L1 · modality mass (COCO early/late + bbox block)\n\n",
        "| Method | image | text | bbox | RF AUC |\n|--------|-------|------|------|--------|\n",
        f"| RF | {br['modality_mass_share']['image']:.3f} | {br['modality_mass_share']['text']:.3f} | "
        f"{br['modality_mass_share']['bbox']:.3f} | {br['domain_auc']:.3f} |\n",
        f"| MMD | {bm['modality_mass_share']['image']:.3f} | {bm['modality_mass_share']['text']:.3f} | "
        f"{bm['modality_mass_share']['bbox']:.3f} | — |\n\n",
        "## L2b · bbox subgroup mass (RF)\n\n",
        f"| geo | category_hist | top_boxes |\n|-----|---------------|----------|\n"
        f"| {br['bbox_subgroup_mass_share']['geo']:.3f} | "
        f"{br['bbox_subgroup_mass_share']['category_hist']:.3f} | "
        f"{br['bbox_subgroup_mass_share']['top_boxes']:.3f} |\n\n",
        "## BBox inject recovery (GT=bbox)\n\n",
        f"- RF: mass_on_bbox={ir['bbox_recovery']['mass_on_bbox']}, "
        f"selAUC={ir['bbox_recovery']['selection_AUC']}, P@20={ir['bbox_recovery']['P@20']}\n"
        f"- MMD: mass_on_bbox={im['bbox_recovery']['mass_on_bbox']}, "
        f"selAUC={im['bbox_recovery']['selection_AUC']}, P@20={im['bbox_recovery']['P@20']}\n",
    ]
    (OUT / "README.md").write_text("".join(md))

    # LaTeX
    lines = [
        "% Hierarchical attribution board with bounding-box module\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Hierarchical attribution board on COCO early/late "
        "(time-order) with an added bounding-box feature block. "
        "L1 reports modality VIMP mass on $\\texttt{image}$/$\\texttt{text}$/$\\texttt{bbox}$.}\n",
        "\\label{tab:bbox-attr-l1-mass}\n\\small\n",
        "\\begin{tabular}{@{}l ccc c@{}}\n\\toprule\n",
        "Method & Image & Text & BBox & RF AUC \\\\\n\\midrule\n",
        f"RF Domain & ${br['modality_mass_share']['image']:.3f}$ & "
        f"${br['modality_mass_share']['text']:.3f}$ & "
        f"${br['modality_mass_share']['bbox']:.3f}$ & ${br['domain_auc']:.3f}$ \\\\\n",
        f"Coord-MMD & ${bm['modality_mass_share']['image']:.3f}$ & "
        f"${bm['modality_mass_share']['text']:.3f}$ & "
        f"${bm['modality_mass_share']['bbox']:.3f}$ & --- \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{BBox-block subgroup mass (L2b) under RF Domain on the same board. "
        "Subgroups: geometry summaries, COCO category histogram, top-$M$ box descriptors.}\n",
        "\\label{tab:bbox-attr-l2b-subgroups}\n\\small\n",
        "\\begin{tabular}{@{}l ccc@{}}\n\\toprule\n",
        "Method & Geo & Category hist & Top boxes \\\\\n\\midrule\n",
        f"RF Domain & ${br['bbox_subgroup_mass_share']['geo']:.3f}$ & "
        f"${br['bbox_subgroup_mass_share']['category_hist']:.3f}$ & "
        f"${br['bbox_subgroup_mass_share']['top_boxes']:.3f}$ \\\\\n",
        f"Coord-MMD & ${bm['bbox_subgroup_mass_share']['geo']:.3f}$ & "
        f"${bm['bbox_subgroup_mass_share']['category_hist']:.3f}$ & "
        f"${bm['bbox_subgroup_mass_share']['top_boxes']:.3f}$ \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{RF Domain named geo ranking inside the bbox block "
        "(high$\\to$low VIMP).}\n",
        "\\label{tab:bbox-attr-geo-rank}\n\\small\n",
        "\\begin{tabular}{@{}r l c@{}}\n\\toprule\n",
        "Rank & Geo feature & VIMP \\\\\n\\midrule\n",
    ]
    for g in br["bbox_geo_ranking"]:
        lines.append(f"{g['rank']} & \\texttt{{{g['name']}}} & ${g['vimp']:.4f}$ \\\\\n")
    lines += [
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{BBox inject recovery on COCO early/late "
        "(GT modality = bbox). Late batch only receives a mean-shift in the bbox block.}\n",
        "\\label{tab:bbox-attr-inject-recovery}\n\\small\n",
        "\\begin{tabular}{@{}l ccc@{}}\n\\toprule\n",
        "Method & Mass on bbox & Sel.\\ AUC & P@20 \\\\\n\\midrule\n",
        f"RF Domain & ${ir['bbox_recovery']['mass_on_bbox']:.3f}$ & "
        f"${ir['bbox_recovery']['selection_AUC']:.3f}$ & "
        f"${ir['bbox_recovery']['P@20']:.3f}$ \\\\\n",
        f"Coord-MMD & ${im['bbox_recovery']['mass_on_bbox']:.3f}$ & "
        f"${im['bbox_recovery']['selection_AUC']:.3f}$ & "
        f"${im['bbox_recovery']['P@20']:.3f}$ \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n",
    ]
    tex = "".join(lines)
    (DOCS / "BBox_Hierarchical_Attribution_tables_only.tex").write_text(tex)
    (OUT / "BBox_Hierarchical_Attribution_tables_only.tex").write_text(tex)
    print(tex)
    print("".join(md))


def main():
    X0, X1, blocks, schema, meta = load_blocks()
    print(
        f"loaded n0={len(X0)} n1={len(X1)} dims={ {k: v[1]-v[0] for k,v in blocks.items()} }",
        flush=True,
    )
    baseline = run_one("coco_time_order_with_bbox", X0, X1, blocks, schema)
    print("baseline L1", baseline["rf_domain"]["modality_mass_share"], flush=True)
    print("baseline L2b", baseline["rf_domain"]["bbox_subgroup_mass_share"], flush=True)

    X0i, X1i = inject_bbox(X0, X1, blocks)
    inject = run_one("coco_time_order_bbox_inject", X0i, X1i, blocks, schema)
    print("inject recovery RF", inject["rf_domain"]["bbox_recovery"], flush=True)

    write_artifacts(baseline, inject)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
