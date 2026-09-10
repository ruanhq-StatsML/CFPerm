#!/usr/bin/env python3
"""Interpretable bbox attribution board · handcrafted named features only.

Business rule: do NOT rank opaque CLIP coordinates (e.g. index 10/768).
BBox evidence must be human-readable feature names.

Feature groups (all named):
  counts      — n_objects, n_people, supercategory counts, ...
  geometry    — max/mean/total area frac, center, spread
  presence    — binary flags for key business classes
  largest     — category / size / position of the largest box

Board:
  L1  modality mass: image_clip | text_clip | bbox_named
      (CLIP blocks only for modality-level mass, never as claimed indices)
  L2  ranking inside bbox_named by feature NAME + VIMP

  python3 scripts/run_bbox_named_attribution.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "img_txt" / "coco_time_order"
ANN = ROOT / "data" / "raw" / "coco" / "instances_train2017.json"
OUT = ROOT / "results" / "bbox_attribution"
DOCS = ROOT / "docs" / "method"

SEED = 2026
N_SUB = 4000
INJECT_ALPHA = 0.75

# Business-facing class flags (interpretable presence features)
KEY_CLASSES = [
    "person",
    "car",
    "truck",
    "bus",
    "bicycle",
    "motorcycle",
    "dog",
    "cat",
    "chair",
    "couch",
    "dining table",
    "tv",
    "laptop",
    "cell phone",
    "bottle",
    "cup",
    "bowl",
    "book",
]

SUPER_GROUPS = [
    "person",
    "vehicle",
    "animal",
    "outdoor",
    "accessory",
    "sports",
    "kitchen",
    "food",
    "furniture",
    "electronic",
    "appliance",
    "indoor",
]


def build_named_bbox_features(meta: pd.DataFrame) -> tuple[np.ndarray, list[str], dict]:
    inst = json.loads(ANN.read_text())
    id2hw = {im["id"]: (im["width"], im["height"]) for im in inst["images"]}
    by = defaultdict(list)
    for a in inst["annotations"]:
        by[a["image_id"]].append(a)
    cat_id2name = {c["id"]: c["name"] for c in inst["categories"]}
    cat_id2super = {c["id"]: c["supercategory"] for c in inst["categories"]}
    name2id = {c["name"]: c["id"] for c in inst["categories"]}

    # Feature name order is the contract for the board.
    names: list[str] = [
        "n_objects",
        "n_people",
        "n_unique_categories",
        "total_area_frac",
        "mean_area_frac",
        "max_area_frac",
        "largest_cx",
        "largest_cy",
        "largest_w_frac",
        "largest_h_frac",
        "mean_cx",
        "mean_cy",
        "std_cx",
        "std_cy",
        "coverage_union_proxy",  # sum area capped — overlaps inflate
    ]
    for g in SUPER_GROUPS:
        names.append(f"count_{g}")
        names.append(f"area_frac_{g}")
    for cls in KEY_CLASSES:
        names.append(f"has_{cls.replace(' ', '_')}")
        names.append(f"count_{cls.replace(' ', '_')}")
        names.append(f"area_frac_{cls.replace(' ', '_')}")
    names += [
        "largest_is_person",
        "largest_is_vehicle",
        "largest_is_animal",
        "largest_is_furniture",
        "largest_is_food",
    ]

    X = np.zeros((len(meta), len(names)), dtype=np.float32)
    name2j = {n: i for i, n in enumerate(names)}

    for i, row in meta.reset_index(drop=True).iterrows():
        iid = int(row["coco_image_id"])
        W, H = id2hw.get(iid, (1, 1))
        anns = by.get(iid, [])
        if not anns:
            continue

        areas, cxs, cys, ws, hs = [], [], [], [], []
        per_cat_count = defaultdict(int)
        per_cat_area = defaultdict(float)
        per_super_count = defaultdict(int)
        per_super_area = defaultdict(float)
        best = None  # (area, cx, cy, wfrac, hfrac, cid)

        for a in anns:
            x, y, w, h = a["bbox"]
            area = (w * h) / max(W * H, 1.0)
            cx = (x + 0.5 * w) / max(W, 1.0)
            cy = (y + 0.5 * h) / max(H, 1.0)
            wfrac = w / max(W, 1.0)
            hfrac = h / max(H, 1.0)
            cid = int(a["category_id"])
            cname = cat_id2name.get(cid, "unknown")
            super_c = cat_id2super.get(cid, "indoor")
            areas.append(area)
            cxs.append(cx)
            cys.append(cy)
            ws.append(wfrac)
            hs.append(hfrac)
            per_cat_count[cname] += 1
            per_cat_area[cname] += area
            per_super_count[super_c] += 1
            per_super_area[super_c] += area
            if best is None or area > best[0]:
                best = (area, cx, cy, wfrac, hfrac, cid)

        areas_a = np.asarray(areas, float)
        cxs_a = np.asarray(cxs, float)
        cys_a = np.asarray(cys, float)

        def setv(n, v):
            X[i, name2j[n]] = float(v)

        setv("n_objects", len(anns))
        setv("n_people", per_cat_count.get("person", 0))
        setv("n_unique_categories", len(per_cat_count))
        setv("total_area_frac", min(float(areas_a.sum()), 5.0))
        setv("mean_area_frac", float(areas_a.mean()))
        setv("max_area_frac", float(areas_a.max()))
        setv("largest_cx", best[1])
        setv("largest_cy", best[2])
        setv("largest_w_frac", best[3])
        setv("largest_h_frac", best[4])
        setv("mean_cx", float(cxs_a.mean()))
        setv("mean_cy", float(cys_a.mean()))
        setv("std_cx", float(cxs_a.std()) if len(cxs_a) > 1 else 0.0)
        setv("std_cy", float(cys_a.std()) if len(cys_a) > 1 else 0.0)
        setv("coverage_union_proxy", min(float(areas_a.sum()), 1.0))

        for g in SUPER_GROUPS:
            setv(f"count_{g}", per_super_count.get(g, 0))
            setv(f"area_frac_{g}", min(per_super_area.get(g, 0.0), 5.0))

        for cls in KEY_CLASSES:
            key = cls.replace(" ", "_")
            setv(f"has_{key}", 1.0 if per_cat_count.get(cls, 0) > 0 else 0.0)
            setv(f"count_{key}", per_cat_count.get(cls, 0))
            setv(f"area_frac_{key}", min(per_cat_area.get(cls, 0.0), 5.0))

        largest_super = cat_id2super.get(best[5], "")
        setv("largest_is_person", 1.0 if largest_super == "person" else 0.0)
        setv("largest_is_vehicle", 1.0 if largest_super == "vehicle" else 0.0)
        setv("largest_is_animal", 1.0 if largest_super == "animal" else 0.0)
        setv("largest_is_furniture", 1.0 if largest_super == "furniture" else 0.0)
        setv("largest_is_food", 1.0 if largest_super == "food" else 0.0)

    groups = {
        "counts_geometry": names[:15],
        "supercategory": [n for n in names if n.startswith("count_") and n.split("count_")[1] in SUPER_GROUPS]
        + [n for n in names if n.startswith("area_frac_") and n.split("area_frac_")[1] in SUPER_GROUPS],
        "key_class": [n for n in names if any(n.endswith(c.replace(" ", "_")) or n.endswith(c.replace(" ", "_")) for c in [])],
    }
    # cleaner group slices by index ranges we built
    schema = {
        "feature_names": names,
        "n_features": len(names),
        "groups": {
            "geometry_counts": list(range(0, 15)),
            "supercategory": list(range(15, 15 + 2 * len(SUPER_GROUPS))),
            "key_class": list(
                range(
                    15 + 2 * len(SUPER_GROUPS),
                    15 + 2 * len(SUPER_GROUPS) + 3 * len(KEY_CLASSES),
                )
            ),
            "largest_flags": list(
                range(
                    15 + 2 * len(SUPER_GROUPS) + 3 * len(KEY_CLASSES),
                    len(names),
                )
            ),
        },
        "note": "Handcrafted interpretable bbox features only. No CLIP coordinate ranking.",
    }
    return X, names, schema


def rf_domain(X0, X1, *, seed=SEED):
    X = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    p = X.shape[1]
    clf = RandomForestClassifier(
        n_estimators=300,
        max_depth=max(3, int(round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(round(np.sqrt(len(X)) // 2))),
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(X, W)
    vimp = clf.feature_importances_.astype(float)
    Xtr, Xte, Wtr, Wte = train_test_split(X, W, test_size=0.25, random_state=seed, stratify=W)
    clf2 = RandomForestClassifier(
        n_estimators=150,
        max_depth=max(3, int(round(np.sqrt(p)))),
        min_samples_leaf=max(1, int(round(np.sqrt(len(X)) // 2))),
        n_jobs=-1,
        random_state=seed + 1,
    )
    clf2.fit(Xtr, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return vimp, auc


def subsample(X0, X1, *, n_max=N_SUB, seed=SEED):
    rng = np.random.default_rng(seed)
    n0 = min(len(X0), n_max // 2)
    n1 = min(len(X1), n_max // 2)
    return (
        X0[rng.choice(len(X0), n0, replace=False)],
        X1[rng.choice(len(X1), n1, replace=False)],
    )


def named_ranking(vimp: np.ndarray, names: list[str], *, k: int = 25):
    order = np.argsort(-vimp)
    rows = []
    for r, j in enumerate(order[:k]):
        rows.append(
            {
                "rank": r + 1,
                "feature": names[int(j)],
                "vimp": float(vimp[int(j)]),
            }
        )
    return rows


def group_mass(vimp, schema):
    out = {}
    for g, idxs in schema["groups"].items():
        out[g] = float(vimp[np.array(idxs, int)].sum())
    tot = sum(out.values()) + 1e-12
    return {g: out[g] / tot for g in out}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(DATA / "df_metadata.csv")
    assert "coco_image_id" in meta.columns
    bbox_X, names, schema = build_named_bbox_features(meta)
    np.save(DATA / "bbox_named_feats.npy", bbox_X)
    (DATA / "bbox_named_feature_names.json").write_text(
        json.dumps({"names": names, "schema": schema}, indent=2)
    )
    print(f"bbox named features: n={len(meta)} d={len(names)}", flush=True)

    img = np.load(DATA / "img_feats.npy").astype(float)
    txt = np.load(DATA / "txt_feats.npy").astype(float)
    w = meta["batch"].astype(int).to_numpy()

    # L1 uses CLIP only as modality blocks; bbox block is named features.
    X = np.hstack([img, txt, bbox_X.astype(float)])
    d_img, d_txt, d_bbox = img.shape[1], txt.shape[1], bbox_X.shape[1]
    blocks = {
        "image_clip": (0, d_img),
        "text_clip": (d_img, d_img + d_txt),
        "bbox_named": (d_img + d_txt, d_img + d_txt + d_bbox),
    }

    X0, X1 = X[w == 0], X[w == 1]
    X0s, X1s = subsample(X0, X1)
    vimp, auc = rf_domain(X0s, X1s)
    mass = {
        k: float(vimp[a:b].sum())
        for k, (a, b) in blocks.items()
    }
    tot = sum(mass.values()) + 1e-12
    mass = {k: mass[k] / tot for k in mass}

    # Interpretable ranking: ONLY inside bbox_named
    a, b = blocks["bbox_named"]
    bbox_vimp = vimp[a:b]
    ranking = named_ranking(bbox_vimp, names, k=25)
    gmass = group_mass(bbox_vimp, schema)

    print("L1 mass", {k: round(v, 4) for k, v in mass.items()}, "AUC", round(auc, 4), flush=True)
    print("bbox group mass", {k: round(v, 4) for k, v in gmass.items()}, flush=True)
    print("top named features:", [(r["feature"], round(r["vimp"], 5)) for r in ranking[:10]], flush=True)

    # Inject only on bbox_named → recovery should land on named bbox block
    rng = np.random.default_rng(SEED)
    X1i = X1.copy()
    u = rng.normal(size=d_bbox)
    u /= np.linalg.norm(u) + 1e-12
    X1i[:, a:b] += INJECT_ALPHA * u
    X0s, X1is = subsample(X0, X1i, seed=SEED + 1)
    vimp_i, auc_i = rf_domain(X0s, X1is, seed=SEED + 5)
    mass_i = {k: float(vimp_i[aa:bb].sum()) for k, (aa, bb) in blocks.items()}
    tot_i = sum(mass_i.values()) + 1e-12
    mass_i = {k: mass_i[k] / tot_i for k in mass_i}
    top20 = np.argsort(-vimp_i)[:20]
    p20 = float(np.mean([(a <= j < b) for j in top20]))
    ranking_i = named_ranking(vimp_i[a:b], names, k=10)

    payload = {
        "principle": "BBox board ranks handcrafted named features only; never CLIP coordinate indices.",
        "n0": int((w == 0).sum()),
        "n1": int((w == 1).sum()),
        "d_bbox_named": d_bbox,
        "feature_names": names,
        "baseline": {
            "rf_domain_auc": round(auc, 4),
            "modality_mass_share": mass,
            "bbox_group_mass_share": gmass,
            "named_feature_ranking": ranking,
        },
        "bbox_named_inject": {
            "rf_domain_auc": round(auc_i, 4),
            "modality_mass_share": mass_i,
            "P@20_on_bbox_named": round(p20, 4),
            "named_feature_ranking": ranking_i,
        },
    }
    (OUT / "bbox_named_attribution_board.json").write_text(json.dumps(payload, indent=2))
    (OUT / "bbox_named_attribution_prototype.py").write_text(
        "# Interpretable bbox attribution · named features only\n"
        f"modality_mass_share = {repr(mass)}\n"
        f"bbox_group_mass_share = {repr(gmass)}\n"
        f"named_feature_ranking_top25 = {repr([(r['feature'], r['vimp']) for r in ranking])}\n"
        f"inject_mass_bbox_named = {mass_i['bbox_named']!r}\n"
    )

    md = [
        "# Interpretable bbox attribution (named features)\n\n",
        "Do **not** claim CLIP index importance. BBox evidence is handcrafted names.\n\n",
        f"## L1 modality mass (RF AUC={auc:.3f})\n\n",
        f"| image_clip | text_clip | bbox_named |\n|---:|---:|---:|\n"
        f"| {mass['image_clip']:.3f} | {mass['text_clip']:.3f} | {mass['bbox_named']:.3f} |\n\n",
        "## BBox group mass\n\n",
        "| group | mass |\n|---|---:|\n",
    ]
    for g, v in gmass.items():
        md.append(f"| {g} | {v:.3f} |\n")
    md += [
        "\n## Named feature ranking (top 15)\n\n",
        "| rank | feature | vimp |\n|---:|---|---:|\n",
    ]
    for r in ranking[:15]:
        md.append(f"| {r['rank']} | `{r['feature']}` | {r['vimp']:.5f} |\n")
    md += [
        "\n## Inject recovery (GT=bbox_named)\n\n",
        f"- mass_on_bbox_named={mass_i['bbox_named']:.3f}, AUC={auc_i:.3f}, P@20={p20:.3f}\n",
    ]
    (OUT / "README.md").write_text("".join(md))

    # LaTeX
    lines = [
        "% Interpretable bbox attribution — handcrafted named features only\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Hierarchical board with an interpretable bbox module. "
        "CLIP blocks are used only for modality-level mass; "
        "bbox evidence is restricted to handcrafted named features "
        "(no CLIP coordinate claims).}\n",
        "\\label{tab:bbox-named-l1}\n\\small\n",
        "\\begin{tabular}{@{}l ccc c@{}}\n\\toprule\n",
        "Method & Image CLIP & Text CLIP & BBox named & RF AUC \\\\\n\\midrule\n",
        f"RF Domain & ${mass['image_clip']:.3f}$ & ${mass['text_clip']:.3f}$ & "
        f"${mass['bbox_named']:.3f}$ & ${auc:.3f}$ \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{BBox named-feature group mass.}\n",
        "\\label{tab:bbox-named-groups}\n\\small\n",
        "\\begin{tabular}{@{}l c@{}}\n\\toprule\n",
        "Group & Mass share \\\\\n\\midrule\n",
    ]
    for g, v in gmass.items():
        lines.append(f"\\texttt{{{g}}} & ${v:.3f}$ \\\\\n")
    lines += [
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Top named bbox features by RF Domain VIMP "
        "(human-readable names only).}\n",
        "\\label{tab:bbox-named-ranking}\n\\small\n",
        "\\begin{tabular}{@{}r l c@{}}\n\\toprule\n",
        "Rank & Feature & VIMP \\\\\n\\midrule\n",
    ]
    for r in ranking[:15]:
        feat = r["feature"].replace("_", "\\_")
        lines.append(f"{r['rank']} & \\texttt{{{feat}}} & ${r['vimp']:.5f}$ \\\\\n")
    lines += [
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Named-bbox inject recovery (GT = bbox\\_named block).}\n",
        "\\label{tab:bbox-named-inject}\n\\small\n",
        "\\begin{tabular}{@{}l ccc@{}}\n\\toprule\n",
        "Method & Mass on bbox named & RF AUC & P@20 \\\\\n\\midrule\n",
        f"RF Domain & ${mass_i['bbox_named']:.3f}$ & ${auc_i:.3f}$ & ${p20:.3f}$ \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n",
    ]
    tex = "".join(lines)
    (DOCS / "BBox_Named_Attribution_tables_only.tex").write_text(tex)
    (OUT / "BBox_Named_Attribution_tables_only.tex").write_text(tex)
    print(tex)
    print("".join(md))


if __name__ == "__main__":
    main()
