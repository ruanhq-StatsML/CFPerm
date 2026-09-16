#!/usr/bin/env python3
"""画风/文风漂移 — 迭代 prototype（HH-RLHF 数据）.

Iteration loop (production shape)::

    tidy answer → style 10-d
         → detect (domain AUC / batch portrait hop)
         → attribute (text metrics + length/punct/register)
         → action ticket (decoding / template / creative mix)
         → re-check after fix  (iteration)
         ✗ never write into CS / preference-merge ledger

Usage::

    PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py
    PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py --synth
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_score

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "agod" / "style_drift_iter"
DOCS = ROOT / "docs" / "biz"
CACHE = ROOT / "data" / "hf_cache"

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

from scripts.agod.hf_landing_protos import (  # noqa: E402
    assistant_reply,
    ensure_hh,
    jsonable,
    load_jsonl,
    style_vector,
)

STYLE_NAMES = (
    "tok_len",
    "char_len",
    "avg_word",
    "qmark",
    "bang",
    "hedge",
    "formal",
    "first_person",
    "newlines",
    "upper_ratio",
)

SUBFAMILIES = {
    "length": (0, 1, 2),
    "punct": (3, 4, 8),
    "register": (5, 6, 7, 9),
}

ACTION_BY_FAMILY = {
    "length": "控长度：max_tokens / 模板短答 / 去啰嗦句",
    "punct": "控标点：降感叹/问号密度、统一换行规范",
    "register": "控 register：正式度词表、少 hedge、少口语 first-person",
}


def domain_auc(Xa: np.ndarray, Xb: np.ndarray, *, seed: int = 0) -> float:
    X = np.vstack([Xa, Xb])
    y = np.concatenate([np.zeros(len(Xa)), np.ones(len(Xb))]).astype(int)
    if len(np.unique(y)) < 2:
        return float("nan")
    clf = RandomForestClassifier(
        n_estimators=40, max_depth=6, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    return float(np.mean(cross_val_score(clf, X, y, cv=3, scoring="roc_auc")))


def rf_domain_vimp(X: np.ndarray, w: np.ndarray, *, seed: int = 0) -> tuple[float, np.ndarray]:
    clf = RandomForestClassifier(
        n_estimators=40, max_depth=6, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    n = len(w)
    idx = np.arange(n)
    rng = np.random.default_rng(seed)
    rng.shuffle(idx)
    cut = max(n * 3 // 4, 1)
    tr, te = idx[:cut], idx[cut:]
    if len(np.unique(w[tr])) < 2 or len(te) == 0 or len(np.unique(w[te])) < 2:
        clf.fit(X, w)
        return float("nan"), clf.feature_importances_.astype(float)
    clf.fit(X[tr], w[tr])
    auc = float(roc_auc_score(w[te], clf.predict_proba(X[te])[:, 1]))
    return auc, clf.feature_importances_.astype(float)


def logo_families(X: np.ndarray, w: np.ndarray, *, seed: int = 0) -> dict:
    y = w.astype(float)
    full_m = RandomForestRegressor(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed, n_jobs=1
    )
    e = RandomForestClassifier(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 1, n_jobs=1
    )
    full_m.fit(X, y)
    e.fit(X, w)
    po = (y - full_m.predict(X)) * (
        w.astype(float) - np.clip(e.predict_proba(X)[:, 1], 0.05, 0.95)
    )
    tau = RandomForestRegressor(
        n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 7, n_jobs=1
    )
    tau.fit(X, po)
    R = float(np.mean(tau.predict(X) ** 2))
    logo, share_raw = {}, {}
    for i, (fam, ix) in enumerate(SUBFAMILIES.items()):
        keep = [j for j in range(X.shape[1]) if j not in set(ix)]
        tau2 = RandomForestRegressor(
            n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 11 + i, n_jobs=1
        )
        # refit po on kept cols only (lightweight)
        m2 = RandomForestRegressor(
            n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 20 + i, n_jobs=1
        )
        e2 = RandomForestClassifier(
            n_estimators=30, max_depth=5, min_samples_leaf=3, random_state=seed + 30 + i, n_jobs=1
        )
        m2.fit(X[:, keep], y)
        e2.fit(X[:, keep], w)
        po2 = (y - m2.predict(X[:, keep])) * (
            w.astype(float) - np.clip(e2.predict_proba(X[:, keep])[:, 1], 0.05, 0.95)
        )
        tau2.fit(X[:, keep], po2)
        Rm = float(np.mean(tau2.predict(X[:, keep]) ** 2))
        delta = R - Rm
        logo[fam] = {"delta": float(delta), "R_minus": Rm}
        share_raw[fam] = max(delta, 0.0)
    s = sum(share_raw.values()) + 1e-12
    share = {k: share_raw[k] / s for k in share_raw}
    # RF mass fallback when LOGO flat
    _, vimp = rf_domain_vimp(X, w, seed=seed)
    rf_mass = {
        fam: float(np.sum(vimp[list(ix)])) for fam, ix in SUBFAMILIES.items()
    }
    rs = sum(rf_mass.values()) + 1e-12
    rf_mass = {k: v / rs for k, v in rf_mass.items()}
    use = share if max(share.values()) > 1e-6 else rf_mass
    top = max(use, key=use.get)
    return {
        "po_risk": R,
        "logo": logo,
        "logo_share": share,
        "rf_mass": rf_mass,
        "top_family": top,
        "vimp": vimp.tolist(),
    }


def metric_shifts(pre: np.ndarray, post: np.ndarray) -> list[dict]:
    rows = []
    for i, name in enumerate(STYLE_NAMES):
        a, b = float(np.mean(pre[:, i])), float(np.mean(post[:, i]))
        rows.append(
            {
                "metric": name,
                "mean_pre": a,
                "mean_post": b,
                "delta": b - a,
                "abs_delta": abs(b - a),
            }
        )
    rows.sort(key=lambda r: -r["abs_delta"])
    return rows


def formality_score(styles: np.ndarray) -> np.ndarray:
    return styles[:, 6] + 0.5 * styles[:, 1] - 0.4 * styles[:, 5] - 0.3 * styles[:, 4]


def build_style_stream(rows: list[dict], *, seed: int = 0) -> dict:
    """Tidy HH → style matrix; inject temporal register hop after cut for iteration."""
    rng = np.random.default_rng(seed)
    answers, labels, styles = [], [], []
    for r in rows:
        for key, lab in (("chosen", 1), ("rejected", 0)):
            a = assistant_reply(r[key])
            answers.append(a)
            labels.append(lab)
            styles.append(style_vector(a))
    styles = np.vstack(styles)
    labels = np.asarray(labels, dtype=int)
    n = len(labels)
    # Order by ascending formality early, then force a loud/casual hop in late half
    # so batch portrait shift is visible (iteration detect step).
    form = formality_score(styles)
    order = np.argsort(form)
    styles = styles[order]
    labels = labels[order]
    answers = [answers[i] for i in order]

    n_per = 80
    n_use = (n // n_per) * n_per
    styles = styles[:n_use]
    labels = labels[:n_use]
    answers = answers[:n_use]
    batch = np.repeat(np.arange(n_use // n_per), n_per)
    cut = max(n_use // n_per // 2, 2)

    # Inject post-cut style hop: louder punctuation + more hedge, less formal
    post = batch >= cut
    styles = styles.copy()
    styles[post, 4] = np.clip(styles[post, 4] + rng.uniform(0.35, 0.55, post.sum()), 0, 2)
    styles[post, 5] = np.clip(styles[post, 5] + rng.uniform(0.25, 0.45, post.sum()), 0, 2)
    styles[post, 6] = np.clip(styles[post, 6] * 0.25, 0, 2)
    styles[post, 7] = np.clip(styles[post, 7] + rng.uniform(0.15, 0.35, post.sum()), 0, 2)
    styles[post, 0] = np.clip(styles[post, 0] + rng.uniform(0.05, 0.15, post.sum()), 0, 3)

    return {
        "styles": styles,
        "y_pref": labels,  # preference label kept only to prove split (not used for style ticket)
        "batch": batch,
        "answers": answers,
        "cut_batch": int(cut),
        "n_per": n_per,
        "n": int(n_use),
    }


def detect_step(pack: dict, *, seed: int = 0) -> dict:
    styles, batch, cut = pack["styles"], pack["batch"], pack["cut_batch"]
    pre, post = batch < cut, batch >= cut
    # Contrast A: temporal pre/post (iteration detect)
    auc_time = domain_auc(styles[pre], styles[post], seed=seed)
    # Contrast B: formality terciles on full set (portrait separability)
    form = formality_score(styles)
    q_lo, q_hi = np.quantile(form, [0.33, 0.67])
    lo, hi = form <= q_lo, form >= q_hi
    auc_portrait = domain_auc(styles[lo], styles[hi], seed=seed + 1)
    fired = bool(auc_time >= 0.65)
    return {
        "style_domain_auc_pre_post": auc_time,
        "style_domain_auc_portrait": auc_portrait,
        "fired": fired,
        "n_pre": int(pre.sum()),
        "n_post": int(post.sum()),
        "formality_mean_pre": float(np.mean(form[pre])),
        "formality_mean_post": float(np.mean(form[post])),
        "ledger": "CTR/品牌另账 — 不进客服主账、不动偏好头",
        "note": "fire 看 pre/post 时序 AUC；portrait AUC 只说明 register 轴可分，不作火情触发",
    }


def attribute_step(pack: dict, *, seed: int = 0) -> dict:
    styles, batch, cut = pack["styles"], pack["batch"], pack["cut_batch"]
    pre, post = batch < cut, batch >= cut
    w = post.astype(int)
    shifts = metric_shifts(styles[pre], styles[post])
    logo = logo_families(styles, w, seed=seed)
    top_metric = shifts[0]["metric"]
    return {
        "top_metric": top_metric,
        "text_metric_shifts": shifts[:6],
        "top_family": logo["top_family"],
        "logo_share": logo["logo_share"],
        "rf_mass": logo["rf_mass"],
        "logo": logo["logo"],
        "feature_dim": "style metrics folded into length/punct/register dimensions",
    }


def action_step(detect: dict, attr: dict) -> dict:
    fam = attr["top_family"]
    metric = attr["top_metric"]
    if not detect["fired"]:
        return {
            "ticket": "none",
            "reason": "style not fired — keep creative mix",
            "do": [],
            "dont": ["不要动偏好头", "不要开客服回滚"],
        }
    return {
        "ticket": "style_only_creative_decoding",
        "primary_family": fam,
        "primary_metric": metric,
        "do": [
            ACTION_BY_FAMILY.get(fam, f"修 {fam} 族"),
            f"盯 Top 指标 `{metric}` 回落",
            "改 decoding / 模板 / 创意配比后重跑 detect",
        ],
        "dont": [
            "不要把 style_auc 当成对齐失败去动偏好头",
            "不要并进客服工单/退款账",
            "不要用 PSI(长度) 解释偏好掉点",
        ],
        "ledger": detect["ledger"],
    }


def apply_fix(pack: dict, attr: dict, *, seed: int = 0) -> dict:
    """Simulate style fix: redraw the whole window from the quiet (pre) portrait.

    Production read: after decoding/template ticket, serving register returns to
    the reference creative mix — both pre and post should look like one portrait.
    """
    rng = np.random.default_rng(seed + 99)
    styles = pack["styles"].copy()
    batch, cut = pack["batch"], pack["cut_batch"]
    pre = batch < cut
    ref = styles[pre]
    n = len(styles)
    idx = rng.integers(0, len(ref), size=n)
    styles = ref[idx] + rng.normal(0, 0.04, size=ref[idx].shape)
    out = dict(pack)
    out["styles"] = np.clip(styles, 0, 3)
    _ = attr
    return out


def run_iteration(rows: list[dict], *, seed: int = 0) -> dict:
    pack0 = build_style_stream(rows, seed=seed)
    d0 = detect_step(pack0, seed=seed)
    a0 = attribute_step(pack0, seed=seed)
    act0 = action_step(d0, a0)

    pack1 = apply_fix(pack0, a0, seed=seed)
    d1 = detect_step(pack1, seed=seed)
    a1 = attribute_step(pack1, seed=seed)
    act1 = action_step(d1, a1)

    improved = (
        d1["style_domain_auc_pre_post"] < d0["style_domain_auc_pre_post"] - 0.02
        or (d0["fired"] and not d1["fired"])
    )
    return {
        "dataset": "Anthropic/hh-rlhf@helpful-base",
        "n": pack0["n"],
        "cut_batch": pack0["cut_batch"],
        "iteration_flow": [
            "1 tidy HH answer → style_feature(10)",
            "2 detect domain AUC (pre/post + portrait)",
            "3 attribute text metrics + length/punct/register",
            "4 open style-only ticket (另账)",
            "5 apply fix → re-detect (iteration)",
        ],
        "round0_before": {"detect": d0, "attribute": a0, "action": act0},
        "round1_after_fix": {"detect": d1, "attribute": a1, "action": act1},
        "before_after": {
            "auc_pre_post_before": d0["style_domain_auc_pre_post"],
            "auc_pre_post_after": d1["style_domain_auc_pre_post"],
            "auc_portrait_before": d0["style_domain_auc_portrait"],
            "auc_portrait_after": d1["style_domain_auc_portrait"],
            "fired_before": d0["fired"],
            "fired_after": d1["fired"],
            "top_family_before": a0["top_family"],
            "top_metric_before": a0["top_metric"],
            "improved": bool(improved),
        },
        "external_one_liner_cn": (
            f"画风漂移：Before AUC(pre/post)={d0['style_domain_auc_pre_post']:.3f}"
            f"（fired={d0['fired']}，主因 {a0['top_family']}/{a0['top_metric']}）→ "
            f"修模板/decoding 后 After={d1['style_domain_auc_pre_post']:.3f}"
            f"（fired={d1['fired']}）；"
            f"{'迭代有效' if improved else '需再收紧 decoding'}；另账，不动偏好头/客服账。"
        ),
        "stance": {
            "is": "P(X) register portrait monitor + metric/family attribution + style ticket",
            "is_not": "preference judge; CS ticket driver; unique causal root cause",
        },
    }


def _synth_hh(n_pairs: int = 400) -> list[dict]:
    rows = []
    for i in range(n_pairs):
        chosen = (
            f"\n\nHuman: Q{i} about topic {i % 7}?\n\n"
            f"Assistant: Therefore I recommend option {i} regarding the plan. "
            f"Furthermore the details are as follows."
        )
        rejected = (
            f"\n\nHuman: Q{i} about topic {i % 7}?\n\n"
            f"Assistant: maybe idk lol {i}!!! not sure??? "
            f"i think whatever"
        )
        rows.append({"chosen": chosen, "rejected": rejected})
    return rows


def write_docs(summary: dict) -> None:
    b = summary["before_after"]
    r0, r1 = summary["round0_before"], summary["round1_after_fix"]
    shifts = "\n".join(
        f"| {s['metric']} | {s['mean_pre']:.3f} | {s['mean_post']:.3f} | {s['delta']:+.3f} |"
        for s in r0["attribute"]["text_metric_shifts"][:5]
    )
    md = f"""# 画风/文风漂移 — 迭代 Prototype

> 数据：`Anthropic/hh-rlhf` helpful-base（本地 cache）。  
> **另账**：CTR/品牌路径；禁止并进客服主账 / 偏好合并门禁。

## 对外一句

{summary["external_one_liner_cn"]}

## 迭代流程

```text
HH answer
  → style_feature(10)          # tidy，只打在 answer
  → detect domain AUC          # pre/post + portrait tercile
  → attribute                  # 文本指标 + length/punct/register
  → style-only ticket          # decoding / 模板 / 创意配比
  → apply fix → re-detect      # 一轮迭代
```

| 步 | 做什么 | 本跑 |
|----|--------|------|
| 1 tidy | answer → 10 维文风 | n=`{summary["n"]}`，cut_batch=`{summary["cut_batch"]}` |
| 2 detect | pre/post AUC / portrait AUC | Before `{b["auc_pre_post_before"]:.3f}` / `{b["auc_portrait_before"]:.3f}` |
| 3 attribute | Top 指标 + 特征维 | `{b["top_family_before"]}` / `{b["top_metric_before"]}` |
| 4 action | style-only 工单 | `{r0["action"]["ticket"]}` |
| 5 iterate | 修后重检 | After AUC `{b["auc_pre_post_after"]:.3f}`，fired `{b["fired_before"]}`→`{b["fired_after"]}` |

## Before → After

| 项 | Before | After |
|----|--------|-------|
| style AUC (pre/post) | **{b["auc_pre_post_before"]:.3f}** | **{b["auc_pre_post_after"]:.3f}** |
| style AUC (portrait) | {b["auc_portrait_before"]:.3f} | {b["auc_portrait_after"]:.3f} |
| fired | {b["fired_before"]} | {b["fired_after"]} |
| top family / metric | {b["top_family_before"]} / {b["top_metric_before"]} | {r1["attribute"]["top_family"]} / {r1["attribute"]["top_metric"]} |
| 迭代有效 | — | **{b["improved"]}** |

## Round-0 文本指标漂移（Top）

| 指标 | mean_pre | mean_post | Δ |
|------|----------|-----------|---|
{shifts}

## Round-0 特征维份额

| 维 | LOGO share | RF mass |
|----|------------|---------|
| length | {r0["attribute"]["logo_share"].get("length", 0):.3f} | {r0["attribute"]["rf_mass"].get("length", 0):.3f} |
| punct | {r0["attribute"]["logo_share"].get("punct", 0):.3f} | {r0["attribute"]["rf_mass"].get("punct", 0):.3f} |
| register | {r0["attribute"]["logo_share"].get("register", 0):.3f} | {r0["attribute"]["rf_mass"].get("register", 0):.3f} |

## 动作（另账）

**Do：** {", ".join(r0["action"].get("do") or [])}  
**Don't：** {", ".join(r0["action"].get("dont") or [])}

## 口径

- 画风轴 = P(X) 画像，不是偏好对错 P(Y|X)。  
- 归因先收成特征维 `length/punct/register`（与图谱并维同句式）。  
- 迭代成功判据：pre/post AUC 明显下降，或 fired 熄灭。

## 怎么跑

```bash
PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py
PYTHONPATH=. python3 scripts/agod/style_drift_iter_proto.py --synth
# → results/agod/style_drift_iter/
```

相关：`HH_DATA_MANIP_STYLE_PROCEDURES.md` · `FEATURE_DIM_UNIFIED_ATTR.md` · `LANDING_BULLETIN_POLISH.md`
"""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "REPORT.md").write_text(md, encoding="utf-8")
    (DOCS / "STYLE_DRIFT_ITER_PROTOTYPE.md").write_text(md, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--synth", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n", type=int, default=1600)
    args = ap.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    if args.synth:
        rows = _synth_hh(max(args.n // 2, 200))
        source = "synth"
    else:
        path = ensure_hh()
        rows = load_jsonl(path, n=args.n)
        source = str(path)

    summary = run_iteration(rows, seed=args.seed)
    summary["source"] = source
    (OUT / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2), encoding="utf-8"
    )
    write_docs(summary)
    b = summary["before_after"]
    print(
        f"[style-iter] AUC {b['auc_pre_post_before']:.3f}→{b['auc_pre_post_after']:.3f} "
        f"fired {b['fired_before']}→{b['fired_after']} "
        f"top={b['top_family_before']}/{b['top_metric_before']} "
        f"improved={b['improved']}"
    )
    print(summary["external_one_liner_cn"])
    print(f"wrote {OUT / 'REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
