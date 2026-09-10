#!/usr/bin/env python3
"""Post-hoc leave-one-bounding-box-out (LOBO) attribution on COCO.

Decision chain (token-attribution analogue):
  1) Format each sample as {coco_image_id, caption, bbx_text: {bbx_k: {box, text, category}}}
  2) Image-level named features via recalculate_bbox_features(boxes)
  3) Select mechanisms with RF Domain + PO-risk VIMP
  4) Locate high-shift images under selected f*
  5) For each image: leave-one-box-out → recompute features → Δ PO contribution

  python3 scripts/run_coco_lobo_bbox_attribution.py
"""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "coco"
DATA = ROOT / "data" / "img_txt" / "coco_time_order"
OUT = ROOT / "results" / "bbox_attribution"
DOCS = ROOT / "docs" / "method"
FORM = DATA / "coco_bbx_text_form"

SEED = 2026
ALPHA_CLIP = 0.05
N_FIT = 2400
N_LOCATE = 40
MAX_BOXES = 8  # per image for LOBO cost
MIN_BOX_AREA = 32 * 32

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
    return float(inter / max(aw * ah + bw * bh - inter, 1e-12))


def mean_std(arr):
    a = np.asarray(arr, float)
    if len(a) == 0:
        return 0.0, 0.0
    if len(a) == 1:
        return float(a[0]), 0.0
    return float(a.mean()), float(a.std())


FEATURE_NAMES = [
    # geometry
    "cx_mean",
    "cy_mean",
    "cx_std",
    "cy_std",
    "w_mean",
    "h_mean",
    "aspect_mean",
    "aspect_std",
    # iou
    "pairwise_iou_mean",
    "pairwise_iou_max",
    "overlap_pair_frac",
    "heavy_overlap_pair_frac",
    "person_person_iou_mean",
    "max_box_overlap_count",
    # crowding
    "iscrowd_ratio",
    "coverage_area",
    "n_objects_norm",
    "obj_density",
    "nn_center_dist_mean",
    "nn_center_dist_min",
    "congest_score",
    "densest_quad_area_prop",
    "periphery_vs_center_area",
    # area
    "area_mean",
    "area_std",
    "largest_area_frac",
    "person_area_frac",
    "vehicle_area_frac",
    "free_area_approx",
    "area_gini",
    # semantic / text×bbox light
    "person_count_norm",
    "vehicle_count_norm",
    "n_categories_present",
    "top1_cat_share",
    "mention_hit_rate",
    "orphan_box_frac",
]


def recalculate_bbox_features(boxes: list[dict], image_hw: tuple[int, int], caption: str = "") -> dict:
    """Recompute image-level named features from the *remaining* box list.

    boxes: [{box:[x,y,w,h], category:str, text:str, iscrowd:float}, ...]
    """
    H, W = image_hw
    out = {n: 0.0 for n in FEATURE_NAMES}
    if not boxes or W <= 0 or H <= 0:
        return out

    cxs, cys, ws, hs, aspects, areas, crowds = [], [], [], [], [], [], []
    xywhs, cats = [], []
    quad = np.zeros(4, float)
    person_area = vehicle_area = 0.0
    cat_counts: Counter = Counter()
    cap = (caption or "").lower()

    for b in boxes:
        x, y, w, h = b["box"]
        cx = (x + 0.5 * w) / W
        cy = (y + 0.5 * h) / H
        wf, hf = w / W, h / H
        area = (w * h) / (W * H)
        ar = max(wf / max(hf, 1e-6), hf / max(wf, 1e-6))
        cname = str(b.get("category", "unk"))
        cxs.append(cx)
        cys.append(cy)
        ws.append(wf)
        hs.append(hf)
        aspects.append(ar)
        areas.append(area)
        crowds.append(float(b.get("iscrowd", 0.0)))
        xywhs.append((x, y, w, h))
        cats.append(cname)
        cat_counts[cname] += 1
        if cname in PERSON:
            person_area += area
        if cname in VEHICLE:
            vehicle_area += area
        q = (0 if cy < 0.5 else 2) + (0 if cx < 0.5 else 1)
        quad[q] += area

    n = len(boxes)
    for key, arr in [
        ("cx", cxs),
        ("cy", cys),
        ("w", ws),
        ("h", hs),
        ("aspect", aspects),
    ]:
        m, s = mean_std(arr)
        if f"{key}_mean" in out:
            out[f"{key}_mean"] = m
        if f"{key}_std" in out:
            out[f"{key}_std"] = s

    ious, person_ious = [], []
    overlap_hits = heavy_hits = 0
    overlap_count = np.zeros(n, int)
    for u in range(n):
        for v in range(u + 1, n):
            val = iou_xywh(xywhs[u], xywhs[v])
            ious.append(val)
            if val > 0:
                overlap_hits += 1
                overlap_count[u] += 1
                overlap_count[v] += 1
            if val > 0.3:
                heavy_hits += 1
            if cats[u] in PERSON and cats[v] in PERSON:
                person_ious.append(val)
    n_pairs = max(n * (n - 1) // 2, 1)
    iou_mean = float(np.mean(ious)) if ious else 0.0
    out["pairwise_iou_mean"] = iou_mean
    out["pairwise_iou_max"] = float(np.max(ious)) if ious else 0.0
    out["overlap_pair_frac"] = overlap_hits / n_pairs
    out["heavy_overlap_pair_frac"] = heavy_hits / n_pairs
    out["person_person_iou_mean"] = float(np.mean(person_ious)) if person_ious else 0.0
    out["max_box_overlap_count"] = float(overlap_count.max()) if n else 0.0

    coverage = float(np.sum(areas))
    nn_mean = nn_min = 0.0
    pts = np.stack([cxs, cys], 1)
    if len(pts) >= 2:
        dmat = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
        np.fill_diagonal(dmat, np.inf)
        nn = dmat.min(1)
        nn_mean = float(nn.mean())
        nn_min = float(nn.min())
    dens_q = int(np.argmax(quad)) if quad.sum() > 0 else 0
    center_a = float(
        sum(
            a
            for cx, cy, a in zip(cxs, cys, areas)
            if 0.25 <= cx <= 0.75 and 0.25 <= cy <= 0.75
        )
    )
    peri_a = max(coverage - center_a, 0.0)
    out["iscrowd_ratio"] = float(np.mean(crowds))
    out["coverage_area"] = coverage
    out["n_objects_norm"] = n / 20.0
    out["obj_density"] = n / max(coverage, 1e-6)
    out["nn_center_dist_mean"] = nn_mean
    out["nn_center_dist_min"] = nn_min
    out["congest_score"] = coverage * (1.0 - nn_mean) * (1.0 + iou_mean)
    out["densest_quad_area_prop"] = float(quad[dens_q] / max(coverage, 1e-6))
    out["periphery_vs_center_area"] = peri_a / max(center_a, 1e-6)

    a_m, a_s = mean_std(areas)
    out["area_mean"] = a_m
    out["area_std"] = a_s
    out["largest_area_frac"] = float(np.max(areas) / max(coverage, 1e-6))
    out["person_area_frac"] = person_area / max(coverage, 1e-6)
    out["vehicle_area_frac"] = vehicle_area / max(coverage, 1e-6)
    out["free_area_approx"] = max(0.0, 1.0 - coverage * (1.0 - 0.5 * iou_mean))
    if len(areas) >= 2:
        s = np.sort(np.asarray(areas, float))
        k = np.arange(1, len(s) + 1)
        gini = float((2 * (k * s).sum()) / (len(s) * s.sum()) - (len(s) + 1) / len(s))
        out["area_gini"] = max(gini, 0.0)

    out["person_count_norm"] = cat_counts.get("person", 0) / 20.0
    out["vehicle_count_norm"] = sum(cat_counts[c] for c in VEHICLE) / 20.0
    out["n_categories_present"] = len(cat_counts) / 20.0
    out["top1_cat_share"] = max(cat_counts.values()) / max(n, 1)

    # light bbox×text: category string appears in caption?
    hits = sum(1 for c in cats if re.search(rf"\b{re.escape(c)}\b", cap))
    out["mention_hit_rate"] = hits / max(n, 1)
    out["orphan_box_frac"] = 1.0 - out["mention_hit_rate"]
    return out


def feat_vector(d: dict) -> np.ndarray:
    return np.array([d[n] for n in FEATURE_NAMES], dtype=np.float64)


def box_text_description(category: str, caption: str) -> str:
    """Per-box text: category + caption snippet if category mentioned, else category label."""
    cap = caption or ""
    if category and re.search(rf"\b{re.escape(category)}\b", cap, flags=re.I):
        return f"{category}: {cap}"
    return f"{category}" if category else (cap[:120] if cap else "object")


def build_bbx_text_samples(meta: pd.DataFrame) -> list[dict]:
    """COCO → list of {image meta, bbx_text} records."""
    inst = json.loads((RAW / "instances_train2017.json").read_text())
    caps = json.loads((RAW / "captions_train2017.json").read_text())
    id2hw = {im["id"]: (im["height"], im["width"]) for im in inst["images"]}
    id2name = {c["id"]: c["name"] for c in inst["categories"]}
    by = defaultdict(list)
    for a in inst["annotations"]:
        by[a["image_id"]].append(a)
    # prefer metadata caption; fallback first coco caption
    cap_by = defaultdict(list)
    for a in caps["annotations"]:
        cap_by[a["image_id"]].append(a["caption"])

    samples = []
    for i, row in meta.reset_index(drop=True).iterrows():
        iid = int(row["coco_image_id"])
        H, W = id2hw.get(iid, (1, 1))
        caption = str(row.get("caption") or (cap_by[iid][0] if cap_by[iid] else ""))
        anns = by.get(iid, [])
        # keep largest boxes up to MAX_BOXES
        ranked = sorted(anns, key=lambda a: a["bbox"][2] * a["bbox"][3], reverse=True)
        boxes = []
        bbx_text = {}
        for k, a in enumerate(ranked[:MAX_BOXES]):
            x, y, w, h = a["bbox"]
            if w * h < MIN_BOX_AREA:
                continue
            cat = id2name.get(int(a["category_id"]), "unk")
            text = box_text_description(cat, caption)
            bname = f"bbx_{k+1}"
            rec = {
                "bbx_id": bname,
                "box": [float(x), float(y), float(w), float(h)],
                "category": cat,
                "text": text,
                "iscrowd": float(a.get("iscrowd", 0)),
                "area": float(w * h),
            }
            boxes.append(rec)
            bbx_text[bname] = {
                "box": rec["box"],
                "text_description": text,
                "category": cat,
            }
        if len(boxes) < 2:
            continue
        samples.append(
            {
                "row_index": int(i),
                "coco_image_id": iid,
                "batch": int(row["batch"]),
                "domain": str(row["domain"]),
                "caption": caption,
                "image_hw": [int(H), int(W)],
                "bbx_text": bbx_text,
                "boxes": boxes,  # working list for LOBO
            }
        )
    return samples


def fit_po_tau(X, Y, W, *, seed=SEED):
    """Cross-fit m,e then fit τ on pseudo-outcome; return τ model + po_risk scalar."""
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    Yf = Y.astype(float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=80,
            max_depth=10,
            min_samples_leaf=3,
            random_state=seed + fold,
            n_jobs=-1,
        )
        e = RandomForestClassifier(
            n_estimators=80,
            max_depth=10,
            min_samples_leaf=3,
            random_state=seed + 40 + fold,
            n_jobs=-1,
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, ALPHA_CLIP, 1 - ALPHA_CLIP)
    po = (Yf - m_hat) * (W - e_hat)
    tau = RandomForestRegressor(
        n_estimators=160,
        max_depth=12,
        min_samples_leaf=3,
        random_state=seed + 7,
        n_jobs=-1,
    )
    tau.fit(X, po)
    pred = tau.predict(X)
    return tau, tau.feature_importances_.astype(float), float(np.mean(pred**2)), pred


def rf_domain_vimp(X, W, *, seed=SEED):
    clf = RandomForestClassifier(
        n_estimators=300,
        max_depth=max(3, int(round(np.sqrt(X.shape[1])))),
        min_samples_leaf=max(1, int(round(np.sqrt(len(X)) // 2))),
        n_jobs=-1,
        random_state=seed,
    )
    clf.fit(X, W)
    return clf.feature_importances_.astype(float)


def lobo_one_image(sample: dict, tau, f_star: list[str]) -> dict:
    """Leave-one-box-out: remove bbx_k (+text) → recalculate → Δ τ² and Δ f*."""
    boxes = sample["boxes"]
    hw = tuple(sample["image_hw"])
    cap = sample["caption"]
    base_feat = recalculate_bbox_features(boxes, hw, cap)
    base_x = feat_vector(base_feat).reshape(1, -1)
    base_tau = float(tau.predict(base_x)[0])
    base_score = base_tau**2

    rows = []
    for k, b in enumerate(boxes):
        remain = boxes[:k] + boxes[k + 1 :]
        feat_k = recalculate_bbox_features(remain, hw, cap)
        xk = feat_vector(feat_k).reshape(1, -1)
        tau_k = float(tau.predict(xk)[0])
        score_k = tau_k**2
        delta_po = base_score - score_k
        # relative drop for readable business scale
        delta_po_rel = delta_po / (base_score + 1e-12)
        delta_f = {f: float(base_feat[f] - feat_k[f]) for f in f_star}
        rows.append(
            {
                "bbx_id": b["bbx_id"],
                "box": b["box"],
                "category": b["category"],
                "text_description": b["text"],
                "delta_po_contrib": delta_po,
                "delta_po_rel": delta_po_rel,
                "delta_f_star": delta_f,
                "abs_delta_po": abs(delta_po),
                "abs_delta_po_rel": abs(delta_po_rel),
            }
        )
    rows.sort(key=lambda r: -r["abs_delta_po_rel"])
    top = rows[0] if rows else None
    return {
        "coco_image_id": sample["coco_image_id"],
        "batch": sample["batch"],
        "caption": sample["caption"],
        "n_boxes": len(boxes),
        "base_po_contrib": base_score,
        "base_f_star": {f: float(base_feat[f]) for f in f_star},
        "top_box": top,
        "ranking": rows,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    FORM.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(DATA / "df_metadata.csv")
    print("building bbx_text form…", flush=True)
    samples = build_bbx_text_samples(meta)
    # persist compact form (no working boxes duplication beyond bbx_text)
    form_path = FORM / "samples_bbx_text.jsonl"
    with form_path.open("w") as f:
        for s in samples:
            rec = {
                "row_index": s["row_index"],
                "coco_image_id": s["coco_image_id"],
                "batch": s["batch"],
                "domain": s["domain"],
                "caption": s["caption"],
                "image_hw": s["image_hw"],
                "bbx_text": s["bbx_text"],
            }
            f.write(json.dumps(rec) + "\n")
    (FORM / "README.md").write_text(
        "# COCO bbx_text form\n\n"
        "Each line: `{coco_image_id, caption, image_hw, bbx_text: {bbx_k: {box, text_description, category}}}`\n"
        f"n_samples={len(samples)}, max_boxes={MAX_BOXES}\n"
    )
    print(f"wrote {form_path} n={len(samples)}", flush=True)

    # feature matrix
    X_all = np.vstack(
        [
            feat_vector(
                recalculate_bbox_features(s["boxes"], tuple(s["image_hw"]), s["caption"])
            )
            for s in samples
        ]
    )
    W_all = np.array([s["batch"] for s in samples], int)
    # outcome Y: labels.npy (row-aligned) or n_boxes; normalize to [0,1] via rank quantile
    labels = np.load(DATA / "labels.npy").astype(float)
    Y_raw = np.array(
        [
            float(labels[s["row_index"]]) if s["row_index"] < len(labels) else float(len(s["boxes"]))
            for s in samples
        ]
    )
    # quantile map → [0, 1] (average ranks for ties)
    order = np.argsort(Y_raw, kind="mergesort")
    ranks = np.empty(len(Y_raw), dtype=float)
    ranks[order] = np.arange(len(Y_raw), dtype=float)
    # tie-average
    _, inv, counts = np.unique(Y_raw, return_inverse=True, return_counts=True)
    sum_ranks = np.zeros(len(counts), dtype=float)
    np.add.at(sum_ranks, inv, ranks)
    mean_ranks = sum_ranks / counts
    ranks = mean_ranks[inv]
    Y_all = ranks / max(len(Y_raw) - 1, 1)

    rng = np.random.default_rng(SEED)
    # balanced subsample for fitting
    i0 = np.where(W_all == 0)[0]
    i1 = np.where(W_all == 1)[0]
    n_each = min(len(i0), len(i1), N_FIT // 2)
    fit_idx = np.concatenate(
        [rng.choice(i0, n_each, replace=False), rng.choice(i1, n_each, replace=False)]
    )
    rng.shuffle(fit_idx)
    X, Y, W = X_all[fit_idx], Y_all[fit_idx], W_all[fit_idx]

    print("RF Domain VIMP…", flush=True)
    rf_vimp = rf_domain_vimp(X, W)
    print("PO-risk τ…", flush=True)
    tau, po_vimp, por, po_pred = fit_po_tau(X, Y, W)

    # select f*: union of top RF and top PO named features
    top_rf = [FEATURE_NAMES[j] for j in np.argsort(-rf_vimp)[:8]]
    top_po = [FEATURE_NAMES[j] for j in np.argsort(-po_vimp)[:8]]
    # prefer crowding/area/iou mechanisms
    pref = [
        n
        for n in (top_po + top_rf)
        if any(
            k in n
            for k in (
                "density",
                "congest",
                "coverage",
                "gini",
                "iou",
                "nn_",
                "largest",
                "person_area",
                "mention",
                "orphan",
            )
        )
    ]
    f_star = []
    for n in pref + top_po + top_rf:
        if n not in f_star:
            f_star.append(n)
        if len(f_star) >= 5:
            break

    rf_rank = [
        {"rank": r + 1, "feature": FEATURE_NAMES[j], "vimp": float(rf_vimp[j])}
        for r, j in enumerate(np.argsort(-rf_vimp)[:15])
    ]
    po_rank = [
        {"rank": r + 1, "feature": FEATURE_NAMES[j], "vimp": float(po_vimp[j])}
        for r, j in enumerate(np.argsort(-po_vimp)[:15])
    ]

    # locate images: high |τ| on full sample matrix using fitted τ
    tau_all = tau.predict(X_all)
    score_all = tau_all**2
    # also mechanism score = sum |zscore(f*)|
    Z = X_all.copy()
    mu, sg = Z.mean(0), Z.std(0) + 1e-12
    Z = (Z - mu) / sg
    f_idx = [FEATURE_NAMES.index(f) for f in f_star]
    mech = np.abs(Z[:, f_idx]).sum(1)
    locate_score = 0.5 * (score_all / (score_all.max() + 1e-12)) + 0.5 * (
        mech / (mech.max() + 1e-12)
    )
    # take from both batches
    loc = []
    for b in (0, 1):
        idx_b = np.where(W_all == b)[0]
        take = idx_b[np.argsort(-locate_score[idx_b])[: N_LOCATE // 2]]
        loc.extend(int(i) for i in take)

    print(f"f_star={f_star}", flush=True)
    print(f"LOBO on {len(loc)} images…", flush=True)
    attributions = []
    for i in loc:
        attributions.append(lobo_one_image(samples[i], tau, f_star))

    # category counters among top boxes
    cat_counter = Counter(a["top_box"]["category"] for a in attributions if a.get("top_box"))
    # aggregate mean |Δ| by category
    cat_delta = defaultdict(list)
    for a in attributions:
        if a.get("top_box"):
            cat_delta[a["top_box"]["category"]].append(a["top_box"]["abs_delta_po_rel"])
    cat_summary = [
        {
            "category": c,
            "n_top": cat_counter[c],
            "mean_abs_delta_po_rel": float(np.mean(cat_delta[c])),
        }
        for c in sorted(cat_counter, key=lambda x: -cat_counter[x])
    ]

    # example table rows
    examples = []
    for a in attributions[:12]:
        t = a["top_box"]
        if not t:
            continue
        examples.append(
            {
                "coco_image_id": a["coco_image_id"],
                "batch": a["batch"],
                "bbx_id": t["bbx_id"],
                "category": t["category"],
                "delta_po_rel": round(t["delta_po_rel"], 4),
                "text": t["text_description"][:80],
            }
        )

    payload = {
        "principle": (
            "Post-hoc LOBO bbox attribution: select named f* with RF/PO VIMP, "
            "locate images, remove one box(+text), recalculate_bbox_features, Δ PO contrib."
        ),
        "form": str(form_path.relative_to(ROOT)),
        "y_normalization": "rank_quantile_to_[0,1]",
        "y_raw_summary": {
            "min": float(Y_raw.min()),
            "max": float(Y_raw.max()),
            "mean": float(Y_raw.mean()),
        },
        "n_samples_form": len(samples),
        "n_fit": int(len(fit_idx)),
        "po_risk": por,
        "f_star": f_star,
        "rf_named_ranking_top15": rf_rank,
        "po_named_ranking_top15": po_rank,
        "n_lobo_images": len(attributions),
        "top_box_category_counter": cat_summary,
        "examples": examples,
        "attributions": [
            {
                "coco_image_id": a["coco_image_id"],
                "batch": a["batch"],
                "caption": a["caption"][:120],
                "n_boxes": a["n_boxes"],
                "base_po_contrib": a["base_po_contrib"],
                "base_f_star": a["base_f_star"],
                "top_box": a["top_box"],
                "ranking_top3": a["ranking"][:3],
            }
            for a in attributions
        ],
    }
    (OUT / "coco_lobo_bbox_attribution_board.json").write_text(json.dumps(payload, indent=2))

    # markdown
    md = [
        "# COCO post-hoc LOBO bounding-box attribution\n\n",
        f"- form: `{form_path.relative_to(ROOT)}`\n",
        f"- PO-risk (fit) = {por:.4f}\n",
        f"- f* = {f_star}\n\n",
        "## Named selection (PO VIMP top)\n\n| rank | feature | vimp |\n|---:|---|---:|\n",
    ]
    for r in po_rank[:10]:
        md.append(f"| {r['rank']} | `{r['feature']}` | {r['vimp']:.5f} |\n")
    md += [
        "\n## Top-box category counter (among LOBO images)\n\n",
        "| category | n_top | mean |Δ PO| |\n|---|---:|---:|\n",
    ]
    for c in cat_summary[:10]:
        md.append(
            f"| {c['category']} | {c['n_top']} | {c['mean_abs_delta_po_rel']:.4f} |\n"
        )
    md += [
        "\n## Examples\n\n",
        "| image | batch | bbx | category | Δ PO rel | text |\n|---:|---:|---|---|---:|---|\n",
    ]
    for e in examples:
        md.append(
            f"| {e['coco_image_id']} | {e['batch']} | {e['bbx_id']} | {e['category']} | "
            f"{e['delta_po_rel']:.4f} | {e['text'].replace('|', '/')} |\n"
        )
    (OUT / "README_lobo_bbox.md").write_text("".join(md))

    # LaTeX
    def esc(s: str) -> str:
        return str(s).replace("_", "\\_").replace("&", "\\&")

    tex = [
        "% COCO post-hoc leave-one-bounding-box-out attribution\n",
        "% Requires booktabs.\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Post-hoc LOBO bbox attribution on COCO early/late. "
        "Named mechanisms $f^{\\star}$ selected by PO-risk VIMP; "
        "images located by $\\widehat\\tau^2$ and $f^{\\star}$; "
        "each box is removed with its text, features recalculated, "
        "and $\\Delta$ PO contribution ranked.}\n",
        "\\label{tab:coco-lobo-pipeline}\n\\small\n",
        "\\begin{tabular}{@{}l l@{}}\n\\toprule\n",
        "Step & Operation \\\\\n\\midrule\n",
        "Form & $\\{$caption, bbx\\_text: box + text\\_description$\\}$ \\\\\n",
        "Select $f^{\\star}$ & PO / RF named VIMP on image-level bbox features \\\\\n",
        "Locate images & high $\\widehat\\tau(X)^2$ and $|f^{\\star}|$ \\\\\n",
        "LOBO & remain $=$ boxes$\\setminus\\{k\\}$; recalculate\\_bbox\\_features \\\\\n",
        "Score & $\\Delta_k^{\\mathrm{rel}}=(\\widehat\\tau(X)^2-\\widehat\\tau(X^{(-k)})^2)/\\widehat\\tau(X)^2$ \\\\\n",
        "Y & rank-quantile mapped to $[0,1]$ \\\\\n",
        f"PO-risk (fit) & ${por:.4f}$ \\\\\n",
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Selected mechanisms $f^{\\star}$ and top PO named ranks.}\n",
        "\\label{tab:coco-lobo-fstar}\n\\small\n",
        "\\begin{tabular}{@{}r l c@{}}\n\\toprule\n",
        "Rank & Feature & PO VIMP \\\\\n\\midrule\n",
    ]
    for r in po_rank[:10]:
        mark = " $*$" if r["feature"] in f_star else ""
        tex.append(
            f"{r['rank']} & \\texttt{{{esc(r['feature'])}}}{mark} & ${r['vimp']:.5f}$ \\\\\n"
        )
    tex += [
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Category counter of LOBO top boxes "
        f"(n={len(attributions)} located images)." + "}\n",
        "\\label{tab:coco-lobo-cat}\n\\small\n",
        "\\begin{tabular}{@{}l c c@{}}\n\\toprule\n",
        "Category & \\# as top box & mean $|\\Delta$ PO rel$|$ \\\\\n\\midrule\n",
    ]
    for c in cat_summary[:8]:
        tex.append(
            f"{esc(c['category'])} & {c['n_top']} & ${c['mean_abs_delta_po_rel']:.4f}$ \\\\\n"
        )
    tex += [
        "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n",
        "\\begin{table}[ht]\n\\centering\n",
        "\\caption{Example images: which bounding box has largest "
        "leave-one-out relative $\\Delta$ PO contribution "
        "$(\\widehat\\tau^2-\\widehat\\tau_{(-k)}^2)/\\widehat\\tau^2$.}\n",
        "\\label{tab:coco-lobo-examples}\n\\small\n",
        "\\begin{tabular}{@{}r c l l c@{}}\n\\toprule\n",
        "Image & Batch & Box & Category & $\\Delta$ PO rel \\\\\n\\midrule\n",
    ]
    for e in examples[:10]:
        tex.append(
            f"{e['coco_image_id']} & {e['batch']} & \\texttt{{{esc(e['bbx_id'])}}} & "
            f"{esc(e['category'])} & ${e['delta_po_rel']:.4f}$ \\\\\n"
        )
    tex += ["\\bottomrule\n\\end{tabular}\n\\end{table}\n"]
    text = "".join(tex)
    (OUT / "COCO_LOBO_BBox_Attribution_tables_only.tex").write_text(text)
    (DOCS / "COCO_LOBO_BBox_Attribution_tables_only.tex").write_text(text)

    print("category_counter", cat_summary[:5], flush=True)
    print("examples", examples[:3], flush=True)
    print("wrote", OUT / "COCO_LOBO_BBox_Attribution_tables_only.tex", flush=True)


if __name__ == "__main__":
    main()
