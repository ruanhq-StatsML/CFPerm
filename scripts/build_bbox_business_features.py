#!/usr/bin/env python3
"""Per-sample interpretable bbox feature table (business schema).

Schema (per sample; NO bbox-count; NO RoI-512 ranking):
  category histogram          80
  position cx,cy mean/std      4
  size w,h mean/std            4
  area mean/std                2
  aspect ratio mean/std        2
  detection score mean/std     2   (GT COCO → score=1 placeholder)
  iscrowd ratio                1
  densest region + NN dist     2
  pairwise IoU mean            1
  -------------------------
  total named dims           ~98

Highly interpretable: every coordinate has a name.
RoI embedding(512) intentionally omitted from business ranking.

  python3 scripts/build_bbox_business_features.py
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
SEED = 2026
N_SUB = 4000


def iou_xywh(a, b):
    ax, ay, aw, ah = a
    bx, by, bw, bh = b
    ax2, ay2 = ax + aw, ay + ah
    bx2, by2 = bx + bw, by + bh
    ix1, iy1 = max(ax, bx), max(ay, by)
    ix2, iy2 = min(ax2, bx2), min(ay2, by2)
    iw, ih = max(0.0, ix2 - ix1), max(0.0, iy2 - iy1)
    inter = iw * ih
    union = aw * ah + bw * bh - inter + 1e-12
    return inter / union


def build_feature_names(cat_ids: list[int]) -> list[str]:
    names = [f"cat_hist_{cid}" for cid in cat_ids]
    names += [
        "cx_mean",
        "cy_mean",
        "cx_std",
        "cy_std",
        "w_mean",
        "h_mean",
        "w_std",
        "h_std",
        "area_mean",
        "area_std",
        "aspect_mean",
        "aspect_std",
        "score_mean",  # placeholder for GT
        "score_std",
        "iscrowd_ratio",
        "densest_quad",  # 0=TL,1=TR,2=BL,3=BR
        "nn_center_dist_mean",
        "pairwise_iou_mean",
    ]
    # also keep a few human aliases for densest quad one-hot (more interpretable than int code)
    names += [
        "densest_is_TL",
        "densest_is_TR",
        "densest_is_BL",
        "densest_is_BR",
    ]
    return names


def build_matrix(meta: pd.DataFrame) -> tuple[np.ndarray, list[str], dict]:
    inst = json.loads((RAW / "instances_train2017.json").read_text())
    id2hw = {im["id"]: (im["width"], im["height"]) for im in inst["images"]}
    by = defaultdict(list)
    for a in inst["annotations"]:
        by[a["image_id"]].append(a)
    cats = sorted({c["id"] for c in inst["categories"]})
    cat2i = {c: i for i, c in enumerate(cats)}
    names = build_feature_names(cats)
    X = np.zeros((len(meta), len(names)), dtype=np.float32)
    name2j = {n: i for i, n in enumerate(names)}

    for i, row in meta.reset_index(drop=True).iterrows():
        iid = int(row["coco_image_id"])
        W, H = id2hw.get(iid, (1, 1))
        anns = by.get(iid, [])
        if not anns:
            continue

        hist = np.zeros(len(cats), float)
        cxs, cys, ws, hs, areas, aspects, scores, crowds = [], [], [], [], [], [], [], []
        boxes = []
        quad_a = np.zeros(4, float)  # TL TR BL BR

        for a in anns:
            x, y, w, h = a["bbox"]
            cid = int(a["category_id"])
            if cid in cat2i:
                hist[cat2i[cid]] += 1
            cx = (x + 0.5 * w) / max(W, 1)
            cy = (y + 0.5 * h) / max(H, 1)
            wf, hf = w / max(W, 1), h / max(H, 1)
            area = (w * h) / max(W * H, 1)
            ar = max(wf / max(hf, 1e-6), hf / max(wf, 1e-6))
            cxs.append(cx)
            cys.append(cy)
            ws.append(wf)
            hs.append(hf)
            areas.append(area)
            aspects.append(ar)
            scores.append(float(a.get("score", 1.0)))  # GT → 1.0
            crowds.append(float(a.get("iscrowd", 0)))
            boxes.append((x, y, w, h))
            q = (0 if cy < 0.5 else 2) + (0 if cx < 0.5 else 1)
            quad_a[q] += area

        if hist.sum() > 0:
            hist = hist / hist.sum()
        for cid in cats:
            X[i, name2j[f"cat_hist_{cid}"]] = hist[cat2i[cid]]

        def st(arr):
            a = np.asarray(arr, float)
            return float(a.mean()), float(a.std()) if len(a) > 1 else 0.0

        cx_m, cx_s = st(cxs)
        cy_m, cy_s = st(cys)
        w_m, w_s = st(ws)
        h_m, h_s = st(hs)
        a_m, a_s = st(areas)
        ar_m, ar_s = st(aspects)
        sc_m, sc_s = st(scores)
        X[i, name2j["cx_mean"]] = cx_m
        X[i, name2j["cy_mean"]] = cy_m
        X[i, name2j["cx_std"]] = cx_s
        X[i, name2j["cy_std"]] = cy_s
        X[i, name2j["w_mean"]] = w_m
        X[i, name2j["h_mean"]] = h_m
        X[i, name2j["w_std"]] = w_s
        X[i, name2j["h_std"]] = h_s
        X[i, name2j["area_mean"]] = a_m
        X[i, name2j["area_std"]] = a_s
        X[i, name2j["aspect_mean"]] = ar_m
        X[i, name2j["aspect_std"]] = ar_s
        X[i, name2j["score_mean"]] = sc_m
        X[i, name2j["score_std"]] = sc_s
        X[i, name2j["iscrowd_ratio"]] = float(np.mean(crowds))

        dens = int(np.argmax(quad_a)) if quad_a.sum() > 0 else 0
        X[i, name2j["densest_quad"]] = dens
        for k, lab in enumerate(["TL", "TR", "BL", "BR"]):
            X[i, name2j[f"densest_is_{lab}"]] = 1.0 if dens == k else 0.0

        # mean nearest-neighbor center distance (normalized)
        pts = np.stack([cxs, cys], 1)
        if len(pts) >= 2:
            dmat = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
            np.fill_diagonal(dmat, np.inf)
            nn = dmat.min(1)
            X[i, name2j["nn_center_dist_mean"]] = float(nn.mean())
        else:
            X[i, name2j["nn_center_dist_mean"]] = 0.0

        # pairwise IoU mean
        if len(boxes) >= 2:
            ious = []
            for u in range(len(boxes)):
                for v in range(u + 1, len(boxes)):
                    ious.append(iou_xywh(boxes[u], boxes[v]))
            X[i, name2j["pairwise_iou_mean"]] = float(np.mean(ious))
        else:
            X[i, name2j["pairwise_iou_mean"]] = 0.0

    schema = {
        "feature_names": names,
        "n_features": len(names),
        "blocks": {
            "category_histogram": [f"cat_hist_{c}" for c in cats],
            "position": ["cx_mean", "cy_mean", "cx_std", "cy_std"],
            "size": ["w_mean", "h_mean", "w_std", "h_std"],
            "area": ["area_mean", "area_std"],
            "aspect": ["aspect_mean", "aspect_std"],
            "score": ["score_mean", "score_std"],
            "occlusion": ["iscrowd_ratio"],
            "spatial_distribution": [
                "densest_quad",
                "nn_center_dist_mean",
                "densest_is_TL",
                "densest_is_TR",
                "densest_is_BL",
                "densest_is_BR",
            ],
            "spatial_relation": ["pairwise_iou_mean"],
        },
        "omitted": {
            "bbox_count": "skipped by design",
            "roi_embedding_512": "omitted from business ranking (not interpretable)",
        },
        "notes": {
            "score": "COCO GT has no detector score; filled with 1.0 (std=0). Replace when using detector outputs.",
            "interpretable": True,
        },
    }
    return X, names, schema


def rf_rank(img, txt, bbox, names, w):
    X = np.hstack([img.astype(float), txt.astype(float), bbox.astype(float)])
    d_img, d_txt, d_b = img.shape[1], txt.shape[1], bbox.shape[1]
    a0 = d_img + d_txt
    X0, X1 = X[w == 0], X[w == 1]
    rng = np.random.default_rng(SEED)
    n0 = min(len(X0), N_SUB // 2)
    n1 = min(len(X1), N_SUB // 2)
    X0 = X0[rng.choice(len(X0), n0, replace=False)]
    X1 = X1[rng.choice(len(X1), n1, replace=False)]
    Xtr = np.vstack([X0, X1])
    W = np.concatenate([np.zeros(len(X0)), np.ones(len(X1))]).astype(int)
    clf = RandomForestClassifier(
        n_estimators=300,
        max_depth=max(3, int(round(np.sqrt(Xtr.shape[1])))),
        min_samples_leaf=max(1, int(round(np.sqrt(len(Xtr)) // 2))),
        n_jobs=-1,
        random_state=SEED,
    )
    clf.fit(Xtr, W)
    vimp = clf.feature_importances_.astype(float)
    mass = {
        "image_clip": float(vimp[:d_img].sum()),
        "text_clip": float(vimp[d_img:a0].sum()),
        "bbox_named": float(vimp[a0:].sum()),
    }
    tot = sum(mass.values()) + 1e-12
    mass = {k: mass[k] / tot for k in mass}
    vb = vimp[a0:]
    order = np.argsort(-vb)
    ranking = [
        {"rank": r + 1, "feature": names[int(j)], "vimp": float(vb[int(j)])}
        for r, j in enumerate(order[:20])
    ]
    Xtr2, Xte, Wtr, Wte = train_test_split(Xtr, W, test_size=0.25, random_state=SEED, stratify=W)
    clf2 = RandomForestClassifier(
        n_estimators=120,
        max_depth=max(3, int(round(np.sqrt(Xtr.shape[1])))),
        min_samples_leaf=max(1, int(round(np.sqrt(len(Xtr)) // 2))),
        n_jobs=-1,
        random_state=SEED + 1,
    )
    clf2.fit(Xtr2, Wtr)
    auc = float(roc_auc_score(Wte, clf2.predict_proba(Xte)[:, 1]))
    return mass, ranking, auc


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(DATA / "df_metadata.csv")
    assert "coco_image_id" in meta.columns
    Xb, names, schema = build_matrix(meta)
    np.save(DATA / "bbox_business_feats.npy", Xb)
    (DATA / "bbox_business_feature_schema.json").write_text(json.dumps(schema, indent=2))

    img = np.load(DATA / "img_feats.npy")
    txt = np.load(DATA / "txt_feats.npy")
    w = meta["batch"].astype(int).to_numpy()
    mass, ranking, auc = rf_rank(img, txt, Xb, names, w)

    payload = {
        "business_value": (
            "Per-sample bbox summary stats are human-readable "
            "(histogram/position/size/area/aspect/crowd/spatial). "
            "VIMP ranks NAMES, not embedding indices."
        ),
        "schema": schema,
        "n_features": len(names),
        "rf_domain_auc": round(auc, 4),
        "modality_mass_share": mass,
        "named_feature_ranking_top20": ranking,
    }
    (OUT / "bbox_business_schema_board.json").write_text(json.dumps(payload, indent=2))

    md = [
        "# Business bbox feature schema (highly interpretable)\n\n",
        "| 类别 | 特征 | 维度 |\n|---|---|---:|\n",
        "| 类别 | category histogram | 80 |\n",
        "| 位置 | cx, cy mean/std | 4 |\n",
        "| 尺寸 | w, h mean/std | 4 |\n",
        "| 面积 | area mean/std | 2 |\n",
        "| 宽高比 | aspect mean/std | 2 |\n",
        "| 置信度 | score mean/std (GT placeholder) | 2 |\n",
        "| 遮挡 | iscrowd 比例 | 1 |\n",
        "| 空间分布 | densest quad + NN dist (+one-hot) | 6 |\n",
        "| 空间关系 | pairwise IoU mean | 1 |\n",
        "| ~~数量~~ | ~~count~~ | skipped |\n",
        "| ~~RoI embedding~~ | ~~512~~ | omitted from business ranking |\n\n",
        f"Total named dims: **{len(names)}**. RF AUC={auc:.3f}. "
        f"mass image/text/bbox={mass['image_clip']:.3f}/{mass['text_clip']:.3f}/{mass['bbox_named']:.3f}\n\n",
        "## Top named ranks\n\n| rank | feature | vimp |\n|---:|---|---:|\n",
    ]
    for r in ranking:
        md.append(f"| {r['rank']} | `{r['feature']}` | {r['vimp']:.5f} |\n")
    (OUT / "README_business_schema.md").write_text("".join(md))

    lines = [
        "% Business interpretable bbox feature schema\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Per-sample interpretable bounding-box feature schema "
        "(highly readable names). BBox count and RoI-512 embeddings are omitted "
        "from business ranking.}\n",
        "\\label{tab:bbox-business-schema}\n\\small\n",
        "\\begin{tabular}{@{}l l r@{}}\n\\toprule\n",
        "Category & Features & Dims \\\\\n\\midrule\n",
        "Category & histogram over COCO classes & 80 \\\\\n",
        "Position & $c_x,c_y$ mean/std & 4 \\\\\n",
        "Size & $w,h$ mean/std & 4 \\\\\n",
        "Area & area mean/std & 2 \\\\\n",
        "Aspect & aspect mean/std & 2 \\\\\n",
        "Score & detection score mean/std (GT placeholder) & 2 \\\\\n",
        "Occlusion & iscrowd ratio & 1 \\\\\n",
        "Spatial distribution & densest quadrant + NN distance (+one-hot) & 6 \\\\\n",
        "Spatial relation & pairwise IoU mean & 1 \\\\\n",
        "\\midrule\n",
        f"Total named & & {len(names)} \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        (
            "\\caption{Top named business bbox features by RF Domain VIMP on COCO "
            "early/late (AUC $"
            + format(auc, ".3f")
            + "$).}\n"
        ),
        "\\label{tab:bbox-business-rank}\n\\small\n",
        "\\begin{tabular}{@{}r l c@{}}\n\\toprule\n",
        "Rank & Feature & VIMP \\\\\n\\midrule\n",
    ]
    for r in ranking[:15]:
        lines.append(
            f"{r['rank']} & \\texttt{{{r['feature'].replace('_', '\\_')}}} & ${r['vimp']:.5f}$ \\\\\n"
        )
    lines += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n"]
    tex = "".join(lines)
    (DOCS / "BBox_Business_Schema_tables_only.tex").write_text(tex)
    (OUT / "BBox_Business_Schema_tables_only.tex").write_text(tex)
    print(tex)
    print("".join(md))


if __name__ == "__main__":
    main()
