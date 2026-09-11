#!/usr/bin/env python3
"""Multi-layer COCO bbox features — beyond conventional mean/std.

Design (business dashboard hierarchy):
  L1  modality mass: image / text / bbox_named
  L2  within-bbox family mass: IoU / crowding / area / semantic / geometry
  L3  named feature ranks inside each family

Non-conventional emphasis (user ask):
  - IoU family: pairwise mean/max, overlap pair frac, person–person IoU
  - Crowding: coverage, object density, NN spacing, iscrowd ratio, congest score
  - Area: total/largest/person/vehicle fracs, overlap-adjusted free area

Geometry mean/std kept only as a thin L0 baseline block.
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
RAW = ROOT / "data" / "raw" / "coco"
DATA = ROOT / "data" / "img_txt" / "coco_time_order"
OUT = ROOT / "results" / "bbox_attribution"
DOCS = ROOT / "docs" / "method"

PERSON = {"person"}
VEHICLE = {"bicycle", "car", "motorcycle", "airplane", "bus", "train", "truck", "boat"}


def iou_xywh(a, b) -> float:
    ax, ay, aw, ah = a
    bx, by, bw, bh = b
    x1, y1 = max(ax, bx), max(ay, by)
    x2, y2 = min(ax + aw, bx + bw), min(ay + ah, by + bh)
    inter = max(0.0, x2 - x1) * max(0.0, y2 - y1)
    if inter <= 0:
        return 0.0
    union = aw * ah + bw * bh - inter
    return float(inter / max(union, 1e-12))


def mean_std(arr):
    a = np.asarray(arr, float)
    if len(a) == 0:
        return 0.0, 0.0
    if len(a) == 1:
        return float(a[0]), 0.0
    return float(a.mean()), float(a.std())


def feature_names() -> list[str]:
    # L0 thin geometry (normalized)
    geo = [
        "cx_mean",
        "cy_mean",
        "cx_std",
        "cy_std",
        "w_mean",
        "h_mean",
        "aspect_mean",
        "aspect_std",
    ]
    # L1 IoU family
    iou = [
        "pairwise_iou_mean",
        "pairwise_iou_max",
        "overlap_pair_frac",  # frac of pairs with IoU>0
        "heavy_overlap_pair_frac",  # IoU>0.3
        "person_person_iou_mean",
        "max_box_overlap_count",  # how many others overlap the busiest box
    ]
    # L2 crowding / congestion
    crowd = [
        "iscrowd_ratio",
        "coverage_area",  # sum(area)/WH clipped? we use sum of normalized areas
        "n_objects_norm",  # n / 20 soft scale
        "obj_density",  # n / max(coverage,eps) — crowded packing
        "nn_center_dist_mean",
        "nn_center_dist_min",
        "congest_score",  # coverage * (1 - nn_mean) * (1 + iou_mean)
        "densest_quad_area_prop",
        "periphery_vs_center_area",
    ]
    # L3 area semantics
    area = [
        "area_mean",
        "area_std",
        "largest_area_frac",
        "person_area_frac",
        "vehicle_area_frac",
        "free_area_approx",  # 1 - coverage*(1-iou_mean) rough non-overlap free
        "area_gini",  # inequality of box areas
    ]
    # L4 light semantic composition
    sem = [
        "person_count_norm",
        "vehicle_count_norm",
        "n_categories_present",
        "top1_cat_share",
    ]
    return geo + iou + crowd + area + sem


def family_of(name: str) -> str:
    if name.startswith(("cx_", "cy_", "w_", "h_", "aspect_")):
        return "geometry"
    if "iou" in name or name.startswith(("overlap_", "heavy_", "max_box_overlap")):
        return "iou"
    if name in {
        "iscrowd_ratio",
        "coverage_area",
        "n_objects_norm",
        "obj_density",
        "nn_center_dist_mean",
        "nn_center_dist_min",
        "congest_score",
        "densest_quad_area_prop",
        "periphery_vs_center_area",
    }:
        return "crowding"
    if "area" in name or name.startswith("free_"):
        return "area"
    return "semantic"


def build_matrix(meta: pd.DataFrame) -> tuple[np.ndarray, list[str], dict]:
    inst = json.loads((RAW / "instances_train2017.json").read_text())
    id2hw = {im["id"]: (im["width"], im["height"]) for im in inst["images"]}
    id2name = {c["id"]: c["name"] for c in inst["categories"]}
    by = defaultdict(list)
    for a in inst["annotations"]:
        by[a["image_id"]].append(a)

    names = feature_names()
    X = np.zeros((len(meta), len(names)), dtype=np.float32)
    j = {n: i for i, n in enumerate(names)}

    for i, row in meta.reset_index(drop=True).iterrows():
        iid = int(row["coco_image_id"])
        W, H = id2hw.get(iid, (1, 1))
        anns = by.get(iid, [])
        if not anns:
            continue

        cxs, cys, ws, hs, aspects, areas, crowds = [], [], [], [], [], [], []
        boxes, names_box = [], []
        quad = np.zeros(4, float)
        cat_counts: dict[str, int] = defaultdict(int)
        person_area = vehicle_area = 0.0

        for a in anns:
            x, y, w, h = a["bbox"]
            cx = (x + 0.5 * w) / max(W, 1)
            cy = (y + 0.5 * h) / max(H, 1)
            wf, hf = w / max(W, 1), h / max(H, 1)
            area = (w * h) / max(W * H, 1)
            ar = max(wf / max(hf, 1e-6), hf / max(wf, 1e-6))
            cname = id2name.get(int(a["category_id"]), "unk")
            cxs.append(cx)
            cys.append(cy)
            ws.append(wf)
            hs.append(hf)
            aspects.append(ar)
            areas.append(area)
            crowds.append(float(a.get("iscrowd", 0)))
            boxes.append((x, y, w, h))
            names_box.append(cname)
            cat_counts[cname] += 1
            if cname in PERSON:
                person_area += area
            if cname in VEHICLE:
                vehicle_area += area
            q = (0 if cy < 0.5 else 2) + (0 if cx < 0.5 else 1)
            quad[q] += area

        # geometry
        for key, arr in [
            ("cx", cxs),
            ("cy", cys),
            ("w", ws),
            ("h", hs),
            ("aspect", aspects),
        ]:
            m, s = mean_std(arr)
            if f"{key}_mean" in j:
                X[i, j[f"{key}_mean"]] = m
            if f"{key}_std" in j:
                X[i, j[f"{key}_std"]] = s

        # IoU family
        ious = []
        person_ious = []
        overlap_hits = 0
        heavy_hits = 0
        n = len(boxes)
        overlap_count = np.zeros(n, int)
        for u in range(n):
            for v in range(u + 1, n):
                val = iou_xywh(boxes[u], boxes[v])
                ious.append(val)
                if val > 0:
                    overlap_hits += 1
                    overlap_count[u] += 1
                    overlap_count[v] += 1
                if val > 0.3:
                    heavy_hits += 1
                if names_box[u] in PERSON and names_box[v] in PERSON:
                    person_ious.append(val)
        n_pairs = max(n * (n - 1) // 2, 1)
        iou_mean = float(np.mean(ious)) if ious else 0.0
        X[i, j["pairwise_iou_mean"]] = iou_mean
        X[i, j["pairwise_iou_max"]] = float(np.max(ious)) if ious else 0.0
        X[i, j["overlap_pair_frac"]] = overlap_hits / n_pairs
        X[i, j["heavy_overlap_pair_frac"]] = heavy_hits / n_pairs
        X[i, j["person_person_iou_mean"]] = (
            float(np.mean(person_ious)) if person_ious else 0.0
        )
        X[i, j["max_box_overlap_count"]] = float(overlap_count.max()) if n else 0.0

        # crowding
        coverage = float(np.sum(areas))
        nn_mean = nn_min = 0.0
        pts = np.stack([cxs, cys], 1) if cxs else np.zeros((0, 2))
        if len(pts) >= 2:
            dmat = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
            np.fill_diagonal(dmat, np.inf)
            nn = dmat.min(1)
            nn_mean = float(nn.mean())
            nn_min = float(nn.min())
        dens_q = int(np.argmax(quad)) if quad.sum() > 0 else 0
        center_a = float(
            np.sum(
                [
                    a
                    for (cx, cy, a) in zip(cxs, cys, areas)
                    if 0.25 <= cx <= 0.75 and 0.25 <= cy <= 0.75
                ]
            )
        )
        peri_a = max(coverage - center_a, 0.0)
        X[i, j["iscrowd_ratio"]] = float(np.mean(crowds)) if crowds else 0.0
        X[i, j["coverage_area"]] = coverage
        X[i, j["n_objects_norm"]] = n / 20.0
        X[i, j["obj_density"]] = n / max(coverage, 1e-6)
        X[i, j["nn_center_dist_mean"]] = nn_mean
        X[i, j["nn_center_dist_min"]] = nn_min
        X[i, j["congest_score"]] = coverage * (1.0 - nn_mean) * (1.0 + iou_mean)
        X[i, j["densest_quad_area_prop"]] = float(quad[dens_q] / max(coverage, 1e-6))
        X[i, j["periphery_vs_center_area"]] = peri_a / max(center_a, 1e-6)

        # area
        a_m, a_s = mean_std(areas)
        X[i, j["area_mean"]] = a_m
        X[i, j["area_std"]] = a_s
        X[i, j["largest_area_frac"]] = float(np.max(areas) / max(coverage, 1e-6)) if areas else 0.0
        X[i, j["person_area_frac"]] = person_area / max(coverage, 1e-6)
        X[i, j["vehicle_area_frac"]] = vehicle_area / max(coverage, 1e-6)
        X[i, j["free_area_approx"]] = max(0.0, 1.0 - coverage * (1.0 - 0.5 * iou_mean))
        # Gini of areas
        if len(areas) >= 2:
            s = np.sort(np.asarray(areas, float))
            k = np.arange(1, len(s) + 1)
            gini = float((2 * (k * s).sum()) / (len(s) * s.sum()) - (len(s) + 1) / len(s))
            X[i, j["area_gini"]] = max(gini, 0.0)
        else:
            X[i, j["area_gini"]] = 0.0

        # semantic
        X[i, j["person_count_norm"]] = cat_counts.get("person", 0) / 20.0
        X[i, j["vehicle_count_norm"]] = sum(cat_counts[c] for c in VEHICLE) / 20.0
        X[i, j["n_categories_present"]] = len(cat_counts) / 20.0
        if cat_counts:
            X[i, j["top1_cat_share"]] = max(cat_counts.values()) / max(n, 1)

    families = {
        "geometry": [n for n in names if family_of(n) == "geometry"],
        "iou": [n for n in names if family_of(n) == "iou"],
        "crowding": [n for n in names if family_of(n) == "crowding"],
        "area": [n for n in names if family_of(n) == "area"],
        "semantic": [n for n in names if family_of(n) == "semantic"],
    }
    schema = {
        "layers": {
            "L1_modality": ["image_clip", "text_clip", "bbox_layered"],
            "L2_bbox_family": list(families.keys()),
            "L3_named_feature": names,
        },
        "families": families,
        "n_features": len(names),
        "emphasis": ["iou", "crowding", "area"],
        "note": "All geometric quantities normalized by image W/H; IoU/crowd/area are primary signal blocks.",
    }
    return X, names, schema


SEED = 0
N_SUB = 4000


def rf_domain(img, txt, bbox, y, seed=SEED):
    """Match business-schema protocol: subsample + held-out AUC."""
    X = np.hstack([img.astype(float), txt.astype(float), bbox.astype(float)])
    d_img, d_txt = img.shape[1], txt.shape[1]
    a0 = d_img + d_txt
    X0, X1 = X[y == 0], X[y == 1]
    rng = np.random.default_rng(seed)
    n0 = min(len(X0), N_SUB // 2)
    n1 = min(len(X1), N_SUB // 2)
    X0 = X0[rng.choice(len(X0), n0, replace=False)]
    X1 = X1[rng.choice(len(X1), n1, replace=False)]
    Xtr = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    depth = max(3, int(round(np.sqrt(Xtr.shape[1]))))
    leaf = max(1, int(round(np.sqrt(len(Xtr)) // 2)))
    clf = RandomForestClassifier(
        n_estimators=300,
        max_depth=depth,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(Xtr, W)
    vimp = clf.feature_importances_.astype(float)
    Xtr2, Xte, Wtr, Wte = train_test_split(
        Xtr, W, test_size=0.25, random_state=seed, stratify=W
    )
    clf2 = RandomForestClassifier(
        n_estimators=120,
        max_depth=depth,
        min_samples_leaf=leaf,
        n_jobs=-1,
        random_state=seed + 1,
    )
    clf2.fit(Xtr2, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return vimp, auc, a0


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(DATA / "df_metadata.csv")
    img = np.load(DATA / "img_feats.npy")
    txt = np.load(DATA / "txt_feats.npy")
    y = meta["batch"].astype(int).to_numpy()  # early/late domain label

    Xb, names, schema = build_matrix(meta)
    np.save(DATA / "bbox_layered_feats.npy", Xb)
    (DATA / "bbox_layered_feature_schema.json").write_text(json.dumps(schema, indent=2))

    vimp, auc, a0 = rf_domain(img, txt, Xb, y)
    d_img, d_txt = img.shape[1], txt.shape[1]
    mass = {
        "image_clip": float(vimp[:d_img].sum()),
        "text_clip": float(vimp[d_img:a0].sum()),
        "bbox_layered": float(vimp[a0:].sum()),
    }
    s = sum(mass.values()) or 1.0
    mass_share = {k: v / s for k, v in mass.items()}

    vb = vimp[a0:]
    # L2 family mass inside bbox
    fam_mass = {}
    for fam, flist in schema["families"].items():
        idx = [names.index(f) for f in flist]
        fam_mass[fam] = float(vb[idx].sum())
    fs = sum(fam_mass.values()) or 1.0
    fam_share = {k: v / fs for k, v in fam_mass.items()}

    order = np.argsort(-vb)
    ranking = [
        {
            "rank": r + 1,
            "feature": names[int(j)],
            "family": family_of(names[int(j)]),
            "vimp": float(vb[int(j)]),
        }
        for r, j in enumerate(order[:20])
    ]

    # within-family top-5 for emphasis families
    family_tops = {}
    for fam in ["iou", "crowding", "area"]:
        flist = schema["families"][fam]
        idx = [names.index(f) for f in flist]
        local = sorted(idx, key=lambda t: -vb[t])
        family_tops[fam] = [
            {"feature": names[t], "vimp": float(vb[t])} for t in local[:5]
        ]

    payload = {
        "dataset": "coco_time_order",
        "rf_domain_auc": auc,
        "modality_mass_share": mass_share,
        "bbox_family_mass_share": fam_share,
        "named_feature_ranking_top20": ranking,
        "emphasis_family_tops": family_tops,
        "schema": schema,
    }
    (OUT / "bbox_layered_coco_board.json").write_text(json.dumps(payload, indent=2))

    # markdown
    md = [
        "# COCO multi-layer bbox features (IoU / crowding / area)\n\n",
        f"RF Domain AUC = {auc:.3f}\n\n",
        "## L1 modality mass\n\n",
        "| block | share |\n|---|---:|\n",
    ]
    for k, v in mass_share.items():
        md.append(f"| {k} | {v:.3f} |\n")
    md += [
        "\n## L2 bbox family mass\n\n",
        "| family | share |\n|---|---:|\n",
    ]
    for k, v in sorted(fam_share.items(), key=lambda kv: -kv[1]):
        md.append(f"| {k} | {v:.3f} |\n")
    md += [
        "\n## L3 top-20 named\n\n",
        "| rank | family | feature | vimp |\n|---:|---|---|---:|\n",
    ]
    for r in ranking:
        md.append(
            f"| {r['rank']} | {r['family']} | `{r['feature']}` | {r['vimp']:.5f} |\n"
        )
    (OUT / "README_layered_coco.md").write_text("".join(md))

    # LaTeX
    def esc(s: str) -> str:
        return s.replace("_", "\\_")

    tex = [
        "% COCO multi-layer bbox features: IoU / crowding / area\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{COCO early/late multi-layer attribution. "
        "L1 modality mass on concat $X=(\\mathrm{image},\\mathrm{text},\\mathrm{bbox\\_layered})$; "
        "L2 family mass inside the bbox block (IoU / crowding / area emphasized).}\n",
        "\\label{tab:coco-layered-mass}\n\\small\n",
        "\\begin{tabular}{@{}l c c@{}}\n\\toprule\n",
        "Level & Block & Mass share \\\\\n\\midrule\n",
        f"L1 modality & Image CLIP & ${mass_share['image_clip']:.3f}$ \\\\\n",
        f"L1 modality & Text CLIP & ${mass_share['text_clip']:.3f}$ \\\\\n",
        f"L1 modality & BBox layered & ${mass_share['bbox_layered']:.3f}$ \\\\\n",
        "\\midrule\n",
    ]
    for fam, v in sorted(fam_share.items(), key=lambda kv: -kv[1]):
        tex.append(f"L2 bbox family & {esc(fam)} & ${v:.3f}$ \\\\\n")
    tex += [
        "\\midrule\n",
        f"RF Domain AUC & \\multicolumn{{2}}{{c}}{{${auc:.3f}$}} \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Emphasis families on COCO: top features inside IoU, crowding, and area "
        "(RF Domain VIMP within the bbox block).}\n",
        "\\label{tab:coco-layered-emphasis}\n\\small\n",
        "\\begin{tabular}{@{}l l c@{}}\n\\toprule\n",
        "Family & Feature & VIMP \\\\\n\\midrule\n",
    ]
    for fam in ["iou", "crowding", "area"]:
        for t in family_tops[fam]:
            tex.append(
                f"{esc(fam)} & \\texttt{{{esc(t['feature'])}}} & ${t['vimp']:.5f}$ \\\\\n"
            )
        tex.append("\\midrule\n")
    # drop last midrule
    tex[-1] = "\\bottomrule\n"
    tex += ["\\end{tabular}\n\\end{table}\n\n"]

    tex += [
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Top-20 named layered bbox features on COCO early/late.}\n",
        "\\label{tab:coco-layered-top20}\n\\small\n",
        "\\begin{tabular}{@{}r l l c@{}}\n\\toprule\n",
        "Rank & Family & Feature & VIMP \\\\\n\\midrule\n",
    ]
    for r in ranking:
        tex.append(
            f"{r['rank']} & {esc(r['family'])} & "
            f"\\texttt{{{esc(r['feature'])}}} & ${r['vimp']:.5f}$ \\\\\n"
        )
    tex += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n"]

    text = "".join(tex)
    (OUT / "COCO_Layered_BBox_tables_only.tex").write_text(text)
    (DOCS / "COCO_Layered_BBox_tables_only.tex").write_text(text)

    print("auc", round(auc, 4))
    print("modality", {k: round(v, 3) for k, v in mass_share.items()})
    print("family", {k: round(v, 3) for k, v in fam_share.items()})
    print("top5", [(r["family"], r["feature"], round(r["vimp"], 5)) for r in ranking[:5]])
    print("wrote", OUT / "COCO_Layered_BBox_tables_only.tex")


if __name__ == "__main__":
    main()
