#!/usr/bin/env python3
"""Named bbox features (region proportions) · concat · cross-dataset VIMP consensus.

Business goal
-------------
Handcraft interpretable bbox features (including region-proportion logic),
concatenate with image/text CLIP, rank by RF VIMP, then check whether a
known feature set can be *selected back* across multiple dataset splits.

If the same named features are recovered on several boards, that is evidence
the feature design has value (not opaque CLIP indices).

Datasets (all COCO CLIP + same named bbox schema):
  1) coco_time_order      — early vs late image_id
  2) coco_outdoor_indoor  — outdoor vs indoor caption keywords
  3) coco_center_split    — center-heavy vs periphery-heavy images (from bbox itself)

Protocol
--------
  A. Baseline VIMP ranking on concat X=[img|txt|bbox_named]
  B. Targeted inject on a known named set S*
  C. Recovery: does S* land in top-k named ranks on each dataset?
  D. Consensus: features recovered on ≥2 datasets

  python3 scripts/run_bbox_named_consensus.py
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
DATA_ROOT = ROOT / "data" / "img_txt"
OUT = ROOT / "results" / "bbox_attribution"
DOCS = ROOT / "docs" / "method"

SEED = 2026
N_SUB = 4000
TOP_K = 20
INJECT_ALPHA = 1.25

KEY_CLASSES = [
    "person",
    "car",
    "truck",
    "bus",
    "bicycle",
    "dog",
    "cat",
    "chair",
    "couch",
    "dining table",
    "tv",
    "laptop",
    "bottle",
    "cup",
    "book",
]

SUPER_GROUPS = [
    "person",
    "vehicle",
    "animal",
    "outdoor",
    "sports",
    "kitchen",
    "food",
    "furniture",
    "electronic",
    "appliance",
    "indoor",
    "accessory",
]

# Known valuable named set we try to recover across datasets
S_STAR = [
    "n_objects",
    "n_people",
    "area_frac_person",
    "area_frac_vehicle",
    "region_center_area_prop",
    "region_periphery_area_prop",
    "quad_TL_area_prop",
    "quad_TR_area_prop",
    "quad_BL_area_prop",
    "quad_BR_area_prop",
]


def _load_coco():
    inst = json.loads((RAW / "instances_train2017.json").read_text())
    caps = json.loads((RAW / "captions_train2017.json").read_text())
    id2hw = {im["id"]: (im["width"], im["height"]) for im in inst["images"]}
    by = defaultdict(list)
    for a in inst["annotations"]:
        by[a["image_id"]].append(a)
    cat_id2name = {c["id"]: c["name"] for c in inst["categories"]}
    cat_id2super = {c["id"]: c["supercategory"] for c in inst["categories"]}
    cap2id = {}
    for a in caps["annotations"]:
        t = " ".join(a["caption"].lower().split())
        cap2id[t] = a["image_id"]
    return id2hw, by, cat_id2name, cat_id2super, cap2id


def ensure_coco_image_id(meta: pd.DataFrame, cap2id: dict) -> pd.DataFrame:
    meta = meta.copy()
    if "coco_image_id" not in meta.columns or meta["coco_image_id"].isna().any():
        keyed = meta["caption"].astype(str).map(lambda s: " ".join(s.lower().split()))
        meta["coco_image_id"] = keyed.map(cap2id)
    meta = meta.dropna(subset=["coco_image_id"]).reset_index(drop=True)
    meta["coco_image_id"] = meta["coco_image_id"].astype(int)
    return meta


def feature_name_list() -> list[str]:
    names = [
        # counts / size
        "n_objects",
        "n_people",
        "n_unique_categories",
        "total_area_frac",
        "mean_area_frac",
        "max_area_frac",
        "coverage_union_proxy",
        # largest box
        "largest_cx",
        "largest_cy",
        "largest_w_frac",
        "largest_h_frac",
        "largest_aspect_ratio",
        # overall layout
        "mean_cx",
        "mean_cy",
        "std_cx",
        "std_cy",
        # region proportion logic (area shares sum≈1 over objects)
        "region_center_area_prop",
        "region_periphery_area_prop",
        "region_center_count_prop",
        "region_periphery_count_prop",
        "quad_TL_area_prop",
        "quad_TR_area_prop",
        "quad_BL_area_prop",
        "quad_BR_area_prop",
        "quad_TL_count_prop",
        "quad_TR_count_prop",
        "quad_BL_count_prop",
        "quad_BR_count_prop",
        # relative proportions
        "person_to_total_area_prop",
        "vehicle_to_total_area_prop",
        "animal_to_total_area_prop",
        "furniture_to_total_area_prop",
        "largest_to_total_area_prop",
    ]
    for g in SUPER_GROUPS:
        names += [f"count_{g}", f"area_frac_{g}"]
    for cls in KEY_CLASSES:
        key = cls.replace(" ", "_")
        names += [f"has_{key}", f"count_{key}", f"area_frac_{key}"]
    names += [
        "largest_is_person",
        "largest_is_vehicle",
        "largest_is_animal",
        "largest_is_furniture",
        "largest_is_food",
    ]
    return names


def build_named_bbox_matrix(meta: pd.DataFrame, coco) -> tuple[np.ndarray, list[str]]:
    id2hw, by, cat_id2name, cat_id2super, _ = coco
    names = feature_name_list()
    name2j = {n: i for i, n in enumerate(names)}
    X = np.zeros((len(meta), len(names)), dtype=np.float32)

    def setv(i, n, v):
        X[i, name2j[n]] = float(v)

    for i, row in meta.reset_index(drop=True).iterrows():
        iid = int(row["coco_image_id"])
        W, H = id2hw.get(iid, (1, 1))
        anns = by.get(iid, [])
        if not anns:
            continue

        areas, cxs, cys = [], [], []
        per_cat_c = defaultdict(int)
        per_cat_a = defaultdict(float)
        per_super_c = defaultdict(int)
        per_super_a = defaultdict(float)
        best = None

        # region accumulators
        center_a = peri_a = 0.0
        center_n = peri_n = 0
        quad_a = {"TL": 0.0, "TR": 0.0, "BL": 0.0, "BR": 0.0}
        quad_n = {"TL": 0, "TR": 0, "BL": 0, "BR": 0}

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
            per_cat_c[cname] += 1
            per_cat_a[cname] += area
            per_super_c[super_c] += 1
            per_super_a[super_c] += area
            if best is None or area > best[0]:
                best = (area, cx, cy, wfrac, hfrac, cid, max(wfrac / max(hfrac, 1e-6), hfrac / max(wfrac, 1e-6)))

            # center = middle 50% box in normalized coords
            if 0.25 <= cx <= 0.75 and 0.25 <= cy <= 0.75:
                center_a += area
                center_n += 1
            else:
                peri_a += area
                peri_n += 1
            q = ("T" if cy < 0.5 else "B") + ("L" if cx < 0.5 else "R")
            quad_a[q] += area
            quad_n[q] += 1

        areas_a = np.asarray(areas, float)
        cxs_a = np.asarray(cxs, float)
        cys_a = np.asarray(cys, float)
        tot_a = float(areas_a.sum()) + 1e-12
        n_obj = float(len(anns))

        setv(i, "n_objects", n_obj)
        setv(i, "n_people", per_cat_c.get("person", 0))
        setv(i, "n_unique_categories", len(per_cat_c))
        setv(i, "total_area_frac", min(tot_a, 5.0))
        setv(i, "mean_area_frac", float(areas_a.mean()))
        setv(i, "max_area_frac", float(areas_a.max()))
        setv(i, "coverage_union_proxy", min(tot_a, 1.0))
        setv(i, "largest_cx", best[1])
        setv(i, "largest_cy", best[2])
        setv(i, "largest_w_frac", best[3])
        setv(i, "largest_h_frac", best[4])
        setv(i, "largest_aspect_ratio", best[6])
        setv(i, "mean_cx", float(cxs_a.mean()))
        setv(i, "mean_cy", float(cys_a.mean()))
        setv(i, "std_cx", float(cxs_a.std()) if len(cxs_a) > 1 else 0.0)
        setv(i, "std_cy", float(cys_a.std()) if len(cys_a) > 1 else 0.0)

        setv(i, "region_center_area_prop", center_a / tot_a)
        setv(i, "region_periphery_area_prop", peri_a / tot_a)
        setv(i, "region_center_count_prop", center_n / n_obj)
        setv(i, "region_periphery_count_prop", peri_n / n_obj)
        for q in ("TL", "TR", "BL", "BR"):
            setv(i, f"quad_{q}_area_prop", quad_a[q] / tot_a)
            setv(i, f"quad_{q}_count_prop", quad_n[q] / n_obj)

        setv(i, "person_to_total_area_prop", per_super_a.get("person", 0.0) / tot_a)
        setv(i, "vehicle_to_total_area_prop", per_super_a.get("vehicle", 0.0) / tot_a)
        setv(i, "animal_to_total_area_prop", per_super_a.get("animal", 0.0) / tot_a)
        setv(i, "furniture_to_total_area_prop", per_super_a.get("furniture", 0.0) / tot_a)
        setv(i, "largest_to_total_area_prop", best[0] / tot_a)

        for g in SUPER_GROUPS:
            setv(i, f"count_{g}", per_super_c.get(g, 0))
            setv(i, f"area_frac_{g}", min(per_super_a.get(g, 0.0), 5.0))
        for cls in KEY_CLASSES:
            key = cls.replace(" ", "_")
            setv(i, f"has_{key}", 1.0 if per_cat_c.get(cls, 0) > 0 else 0.0)
            setv(i, f"count_{key}", per_cat_c.get(cls, 0))
            setv(i, f"area_frac_{key}", min(per_cat_a.get(cls, 0.0), 5.0))

        largest_super = cat_id2super.get(best[5], "")
        for flag, super_name in [
            ("largest_is_person", "person"),
            ("largest_is_vehicle", "vehicle"),
            ("largest_is_animal", "animal"),
            ("largest_is_furniture", "furniture"),
            ("largest_is_food", "food"),
        ]:
            setv(i, flag, 1.0 if largest_super == super_name else 0.0)

    return X, names


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


def named_rank(vimp_bbox, names, *, k=TOP_K):
    order = np.argsort(-vimp_bbox)
    return [
        {"rank": r + 1, "feature": names[int(j)], "vimp": float(vimp_bbox[int(j)])}
        for r, j in enumerate(order[:k])
    ]


def prepare_center_split(meta: pd.DataFrame, bbox_X: np.ndarray, names: list[str]):
    """Third dataset: batch by center vs periphery area proportion."""
    j = names.index("region_center_area_prop")
    score = bbox_X[:, j]
    # drop empties (no boxes → 0/0 already 0)
    med = np.median(score)
    batch = (score > med).astype(int)  # 1 = center-heavy
    meta = meta.copy()
    meta["batch"] = batch
    meta["domain"] = np.where(batch == 0, "periphery_heavy", "center_heavy")
    return meta


def run_board(tag, img, txt, bbox_X, names, w, *, inject_s_star=False):
    X = np.hstack([img.astype(float), txt.astype(float), bbox_X.astype(float)])
    d_img, d_txt, d_bbox = img.shape[1], txt.shape[1], bbox_X.shape[1]
    a0, a1 = d_img + d_txt, d_img + d_txt + d_bbox
    name2j = {n: i for i, n in enumerate(names)}

    X0, X1 = X[w == 0], X[w == 1]
    if inject_s_star:
        rng = np.random.default_rng(SEED + 17)
        X1 = X1.copy()
        # shift only S* coordinates inside bbox block
        for n in S_STAR:
            j = a0 + name2j[n]
            # signed mean-shift on that named feature
            X1[:, j] += INJECT_ALPHA * (0.5 + rng.random())
        mode = "inject_S_star"
    else:
        mode = "baseline"

    X0s, X1s = subsample(X0, X1, seed=SEED + (3 if inject_s_star else 0))
    vimp, auc = rf_domain(X0s, X1s, seed=SEED + (9 if inject_s_star else 0))
    mass = {
        "image_clip": float(vimp[:d_img].sum()),
        "text_clip": float(vimp[d_img:a0].sum()),
        "bbox_named": float(vimp[a0:a1].sum()),
    }
    tot = sum(mass.values()) + 1e-12
    mass = {k: mass[k] / tot for k in mass}
    ranking = named_rank(vimp[a0:a1], names, k=TOP_K)
    top_feats = [r["feature"] for r in ranking]
    recovered = [f for f in S_STAR if f in top_feats]
    # also rank positions of each S* feature
    full_order = list(np.argsort(-vimp[a0:a1]))
    s_ranks = {f: int(full_order.index(name2j[f]) + 1) for f in S_STAR}

    out = {
        "dataset": tag,
        "mode": mode,
        "n0": int((w == 0).sum()),
        "n1": int((w == 1).sum()),
        "rf_domain_auc": round(auc, 4),
        "modality_mass_share": mass,
        "named_feature_ranking_top20": ranking,
        "S_star_recovered_in_top20": recovered,
        "S_star_recovery_rate": round(len(recovered) / len(S_STAR), 4),
        "S_star_rank_positions": s_ranks,
    }
    print(
        f"[{tag}/{mode}] AUC={auc:.3f} mass_bbox={mass['bbox_named']:.3f} "
        f"recover={len(recovered)}/{len(S_STAR)} {recovered}",
        flush=True,
    )
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    coco = _load_coco()
    _, _, _, _, cap2id = coco

    boards = []

    # --- dataset 1: time order ---
    d1 = DATA_ROOT / "coco_time_order"
    meta1 = ensure_coco_image_id(pd.read_csv(d1 / "df_metadata.csv"), cap2id)
    # align arrays to filtered meta rows: original file order assumed same length
    # re-read raw and filter by matched ids present in meta1 index from original
    meta1_full = pd.read_csv(d1 / "df_metadata.csv")
    meta1_full = ensure_coco_image_id(meta1_full, cap2id)
    img1 = np.load(d1 / "img_feats.npy")[: len(meta1_full)]
    txt1 = np.load(d1 / "txt_feats.npy")[: len(meta1_full)]
    bbox1, names = build_named_bbox_matrix(meta1_full, coco)
    np.save(d1 / "bbox_named_feats.npy", bbox1)
    (d1 / "bbox_named_feature_names.json").write_text(
        json.dumps({"names": names, "S_star": S_STAR}, indent=2)
    )
    w1 = meta1_full["batch"].astype(int).to_numpy()
    boards.append(run_board("coco_time_order", img1, txt1, bbox1, names, w1, inject_s_star=False))
    boards.append(run_board("coco_time_order", img1, txt1, bbox1, names, w1, inject_s_star=True))

    # --- dataset 2: outdoor/indoor ---
    d2 = DATA_ROOT / "coco_outdoor_indoor"
    meta2 = ensure_coco_image_id(pd.read_csv(d2 / "df_metadata.csv"), cap2id)
    # outdoor_indoor arrays may be subset already aligned with csv rows
    img2 = np.load(d2 / "img_feats.npy")
    txt2 = np.load(d2 / "txt_feats.npy")
    n2 = min(len(meta2), len(img2), len(txt2))
    meta2 = meta2.iloc[:n2].reset_index(drop=True)
    img2, txt2 = img2[:n2], txt2[:n2]
    bbox2, names2 = build_named_bbox_matrix(meta2, coco)
    assert names2 == names
    np.save(d2 / "bbox_named_feats.npy", bbox2)
    w2 = meta2["batch"].astype(int).to_numpy()
    boards.append(run_board("coco_outdoor_indoor", img2, txt2, bbox2, names, w2, inject_s_star=False))
    boards.append(run_board("coco_outdoor_indoor", img2, txt2, bbox2, names, w2, inject_s_star=True))

    # --- dataset 3: center vs periphery (built from time-order features) ---
    meta3 = prepare_center_split(meta1_full, bbox1, names)
    out3 = DATA_ROOT / "coco_center_split"
    out3.mkdir(parents=True, exist_ok=True)
    np.save(out3 / "img_feats.npy", img1.astype(np.float32))
    np.save(out3 / "txt_feats.npy", txt1.astype(np.float32))
    np.save(out3 / "bbox_named_feats.npy", bbox1)
    meta3.to_csv(out3 / "df_metadata.csv", index=False)
    (out3 / "README.md").write_text(
        "# COCO center vs periphery split\n\n"
        "Batch from `region_center_area_prop` median split. "
        "Same named bbox schema as other COCO boards.\n"
    )
    w3 = meta3["batch"].astype(int).to_numpy()
    boards.append(run_board("coco_center_split", img1, txt1, bbox1, names, w3, inject_s_star=False))
    boards.append(run_board("coco_center_split", img1, txt1, bbox1, names, w3, inject_s_star=True))

    # --- consensus on inject recoveries ---
    inject_boards = [b for b in boards if b["mode"] == "inject_S_star"]
    baseline_boards = [b for b in boards if b["mode"] == "baseline"]
    # feature -> datasets where recovered in top20 under inject
    consensus = {}
    for f in S_STAR:
        hits = [b["dataset"] for b in inject_boards if f in b["S_star_recovered_in_top20"]]
        consensus[f] = {
            "n_datasets_recovered": len(hits),
            "datasets": hits,
            "consensus_ok": len(hits) >= 2,
        }
    consensus_features = [f for f, v in consensus.items() if v["consensus_ok"]]

    # also: natural top20 intersection across baseline boards (optional signal)
    base_sets = [set(r["feature"] for r in b["named_feature_ranking_top20"]) for b in baseline_boards]
    natural_inter = set.intersection(*base_sets) if base_sets else set()

    payload = {
        "principle": "Named bbox features + concat VIMP; value = cross-dataset recovery of S*.",
        "S_star": S_STAR,
        "boards": boards,
        "inject_consensus": consensus,
        "consensus_features_recovered_on_ge_2_datasets": consensus_features,
        "natural_top20_intersection_baseline": sorted(natural_inter),
    }
    (OUT / "bbox_named_consensus_board.json").write_text(json.dumps(payload, indent=2))

    # markdown + latex
    md = [
        "# Named bbox concat VIMP · cross-dataset consensus\n\n",
        f"S* = `{S_STAR}`\n\n",
        "## Inject recovery by dataset\n\n",
        "| dataset | AUC | mass_bbox | recover | recovered features |\n|---|---:|---:|---:|---|\n",
    ]
    for b in inject_boards:
        md.append(
            f"| {b['dataset']} | {b['rf_domain_auc']:.3f} | {b['modality_mass_share']['bbox_named']:.3f} | "
            f"{b['S_star_recovery_rate']:.2f} | {', '.join(b['S_star_recovered_in_top20']) or '—'} |\n"
        )
    md += [
        "\n## Consensus (recovered on ≥2 datasets under inject)\n\n",
        f"**{len(consensus_features)}/{len(S_STAR)}**: {', '.join(consensus_features) or '—'}\n\n",
        "## Baseline top-20 named ranks (for reference)\n\n",
    ]
    for b in baseline_boards:
        top = ", ".join(r["feature"] for r in b["named_feature_ranking_top20"][:8])
        md.append(
            f"- `{b['dataset']}` AUC={b['rf_domain_auc']:.3f} "
            f"mass_bbox={b['modality_mass_share']['bbox_named']:.3f}: {top}\n"
        )
    (OUT / "README.md").write_text("".join(md))

    lines = [
        "% Named bbox features with region proportions · cross-dataset VIMP consensus\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Cross-dataset recovery of handcrafted bbox features $S^{\\star}$ "
        "(region proportions included). Features are concatenated with image/text CLIP; "
        "RF Domain VIMP ranks \\emph{named} bbox features only. "
        "A feature counts as recovered if it appears in the top-20 named ranks "
        "after a targeted late-batch inject on $S^{\\star}$.}\n",
        "\\label{tab:bbox-named-consensus-recovery}\n\\small\n",
        "\\begin{tabular}{@{}l c c c p{6.5cm}@{}}\n\\toprule\n",
        "Dataset & RF AUC & BBox mass & Recover & Recovered $S^{\\star}$ features \\\\\n\\midrule\n",
    ]
    for b in inject_boards:
        rec = ", ".join(f.replace("_", "\\_") for f in b["S_star_recovered_in_top20"]) or "---"
        lines.append(
            f"{b['dataset'].replace('_', '\\_')} & ${b['rf_domain_auc']:.3f}$ & "
            f"${b['modality_mass_share']['bbox_named']:.3f}$ & "
            f"${b['S_star_recovery_rate']:.2f}$ & {rec} \\\\\n"
        )
    lines += [
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Consensus features: members of $S^{\\star}$ recovered on at least two datasets.}\n",
        "\\label{tab:bbox-named-consensus-features}\n\\small\n",
        "\\begin{tabular}{@{}l c l@{}}\n\\toprule\n",
        "Feature & \\# datasets & Datasets \\\\\n\\midrule\n",
    ]
    for f, v in consensus.items():
        ds = ", ".join(d.replace("_", "\\_") for d in v["datasets"]) or "---"
        lines.append(f"\\texttt{{{f.replace('_', '\\_')}}} & {v['n_datasets_recovered']} & {ds} \\\\\n")
    lines += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n\n"]

    # baseline mass table
    lines += [
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Baseline modality mass on concat $X=(\\mathrm{img},\\mathrm{txt},\\mathrm{bbox\\_named})$ "
        "before inject.}\n",
        "\\label{tab:bbox-named-baseline-mass}\n\\small\n",
        "\\begin{tabular}{@{}l cccc@{}}\n\\toprule\n",
        "Dataset & Image CLIP & Text CLIP & BBox named & RF AUC \\\\\n\\midrule\n",
    ]
    for b in baseline_boards:
        m = b["modality_mass_share"]
        lines.append(
            f"{b['dataset'].replace('_', '\\_')} & ${m['image_clip']:.3f}$ & ${m['text_clip']:.3f}$ & "
            f"${m['bbox_named']:.3f}$ & ${b['rf_domain_auc']:.3f}$ \\\\\n"
        )
    lines += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n"]

    tex = "".join(lines)
    (DOCS / "BBox_Named_Consensus_tables_only.tex").write_text(tex)
    (OUT / "BBox_Named_Consensus_tables_only.tex").write_text(tex)
    (OUT / "bbox_named_consensus_prototype.py").write_text(
        "# Named bbox concat VIMP · cross-dataset consensus\n"
        f"S_star = {repr(S_STAR)}\n"
        f"consensus_features = {repr(consensus_features)}\n"
        f"inject_recovery = {repr({b['dataset']: b['S_star_recovered_in_top20'] for b in inject_boards})}\n"
    )
    print(tex)
    print("".join(md))
    print("consensus_features:", consensus_features, flush=True)


if __name__ == "__main__":
    main()
