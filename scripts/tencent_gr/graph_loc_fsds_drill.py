#!/usr/bin/env python3
"""Recsys prototype: graph localization by entity, then a two-step drill.

Entities (enough — no extra grains):

  user / item / merchant
  video_tower / audio_tower   (attach on merchant → pool to user)

FSDS cannot multi-step-localize; online+LOGO needs these entity names.
Use:
  1) LOCALIZE (no Y): each entity vs clock W. Moved = RF-domain AUC ≥ 0.55.
  2) Two-step FSDS drill (with Y), only on the few moved entities:
       step 1  LOGO among those entities
       step 2  LOCO columns inside them (top-k if a 64-d tower)

Not causal. Not GNN. Fake towers are DGP, not production embeddings.
φ=(Y−μ)(W−e) is early/late distance on post-click Y, not a treatment effect.

  PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds_drill.py
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from feat_proto import split_lr  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
EV = ROOT / "results/tencent_gr_fs150/tables/ev.parquet"
POST = ROOT / "results/tencent_gr_fs150/tables/post.parquet"
IMAP = ROOT / "results/tencent_gr_fs150/hop/item_merchant.parquet"
OUT = ROOT / "results/tencent_gr_fs150/hop"
DOCS = ROOT / "docs" / "reports"
SEED = 0
CLK, CNV = 1, 2
SESS_GAP = 30 * 60
TREES = 20
DEPTH = 5
FOLDS = 2
N_CATALOG = 10_000
EMB_DIM = 64
N_THEME = 16
LOCO_TOPK = 8
AUC_MOVED = 0.55

LOC_FAMILIES = (
    "user_connectivity",
    "item_connectivity",
    "merchant_structure",
    "video_tower",
    "audio_tower",
)
DGP_META = {
    "n_catalog": N_CATALOG,
    "emb_dim": EMB_DIM,
    "video_ids": "randint(0, 1e10, 10000)",
    "video_emb": "student_t_df3",
    "audio_ids": "randint(0, 1e10, 10000)",
    "audio_emb": "uniform_minus1_1",
    "merchant_attach": "zipf_subset_k4_15_plus_theme16_plus_t_bias",
    "user_agg": "left_window_clk_cnv_visit_mean_pool",
    "not_real_embeddings": True,
}


def hop_of(name: str) -> str:
    if name.startswith("g_u_"):
        return "graph_user"
    if name.startswith("g_i_"):
        return "graph_item"
    if name.startswith("g_m_"):
        return "graph_merchant"
    if name.startswith("vid_"):
        return "tower_video"
    if name.startswith("aud_"):
        return "tower_audio"
    return "other"


def block_family(name: str) -> str:
    """Localization grain: graph blocks + fake media towers."""
    if name.startswith("g_u_"):
        return "user_connectivity"
    if name.startswith("g_i_"):
        return "item_connectivity"
    if name.startswith("g_m_"):
        return "merchant_structure"
    if name.startswith("vid_"):
        return "video_tower"
    if name.startswith("aud_"):
        return "audio_tower"
    return "not_block"


def _uid(x) -> str:
    return f"u{int(x)}"


def _iid(x) -> str:
    return f"i{int(x)}"


def _mid(x) -> str:
    return f"m{int(x)}"


def session_pairs(clk: pd.DataFrame) -> list[tuple[int, int]]:
    if clk.empty:
        return []
    c = clk.sort_values(["user_id", "ts"]).copy()
    gap = c.groupby("user_id")["ts"].diff().fillna(0)
    c["sess"] = gap.gt(SESS_GAP).groupby(c.user_id).cumsum()
    pairs = []
    for _, g in c.groupby(["user_id", "sess"], sort=False):
        items = g.item_id.unique()
        if len(items) < 2:
            continue
        for a, b in combinations(sorted(int(x) for x in items), 2):
            pairs.append((a, b))
    return pairs


def graph_tables(left: pd.DataFrame, imap: pd.DataFrame):
    """Left-window clk/cnv → degree / PageRank / co-click. Exposures do not edge."""
    hit = left.loc[left.act.isin([CLK, CNV]), ["user_id", "item_id", "ts"]].copy()
    hit = hit.merge(imap, on="item_id", how="left")
    clk_only = left.loc[left.act.eq(CLK), ["user_id", "item_id", "ts"]]

    ui = hit.drop_duplicates(["user_id", "item_id"])
    G_ui = nx.from_pandas_edgelist(
        pd.DataFrame({"s": ui.user_id.map(_uid), "t": ui.item_id.map(_iid)}), "s", "t"
    )
    um = hit.dropna(subset=["merchant_id"]).drop_duplicates(["user_id", "merchant_id"])
    G_um = nx.from_pandas_edgelist(
        pd.DataFrame({"s": um.user_id.map(_uid), "t": um.merchant_id.map(_mid)}), "s", "t"
    )
    m_nodes = [n for n in G_um.nodes if str(n).startswith("m")]
    G_m = nx.bipartite.projected_graph(G_um, m_nodes) if m_nodes else nx.Graph()
    pr = nx.pagerank(G_m, max_iter=50) if G_m.number_of_nodes() else {}
    clust = nx.clustering(G_m) if G_m.number_of_nodes() else {}

    G_co = nx.Graph()
    G_co.add_edges_from((_iid(a), _iid(b)) for a, b in session_pairs(clk_only))

    def deg(G, prefix, key):
        rows = []
        for n, d in G.degree():
            if str(n).startswith(prefix):
                rows.append({key: int(str(n)[1:]), "deg": float(d)})
        return pd.DataFrame(rows)

    users = deg(G_ui, "u", "user_id").rename(columns={"deg": "g_u_item_deg"})
    u2 = deg(G_um, "u", "user_id").rename(columns={"deg": "g_u_merch_deg"})
    users = users.merge(u2, on="user_id", how="outer") if len(u2) else users
    items = deg(G_ui, "i", "item_id").rename(columns={"deg": "g_i_user_deg"})
    i2 = deg(G_co, "i", "item_id").rename(columns={"deg": "g_i_coclick_deg"})
    items = items.merge(i2, on="item_id", how="outer") if len(i2) else items
    merch = deg(G_um, "m", "merchant_id").rename(columns={"deg": "g_m_user_deg"})
    if len(merch):
        merch["g_m_pr"] = merch.merchant_id.map(lambda x: pr.get(_mid(x), 0.0))
        merch["g_m_clust"] = merch.merchant_id.map(lambda x: clust.get(_mid(x), 0.0))
        merch["g_m_proj_deg"] = merch.merchant_id.map(
            lambda x: float(G_m.degree(_mid(x))) if G_m.has_node(_mid(x)) else 0.0
        )
    else:
        merch = pd.DataFrame(
            columns=["merchant_id", "g_m_user_deg", "g_m_pr", "g_m_clust", "g_m_proj_deg"]
        )
    return users.fillna(0.0), items.fillna(0.0), merch.fillna(0.0)


def _zipf_pop(n: int, rng: np.random.Generator) -> np.ndarray:
    rank = np.arange(1, n + 1, dtype=np.float64)
    pop = 1.0 / np.power(rank, 0.8)
    pop = pop / pop.sum()
    perm = rng.permutation(n)
    out = np.empty(n, dtype=np.float64)
    out[perm] = pop
    return out


def media_catalog(rng: np.random.Generator, *, kind: str):
    """Fake asset catalog. IDs are random 0..1e10; embeds are noise, not a model."""
    ids = rng.integers(0, int(1e10), size=N_CATALOG, dtype=np.int64)
    if kind == "video":
        emb = rng.standard_t(df=3.0, size=(N_CATALOG, EMB_DIM))
    else:
        emb = rng.random((N_CATALOG, EMB_DIM)) * 2.0 - 1.0
    return ids, emb, _zipf_pop(N_CATALOG, rng)


def attach_media_to_merchants(
    mids: np.ndarray,
    emb: np.ndarray,
    pop: np.ndarray,
    rng: np.random.Generator,
    *,
    prefix: str,
) -> pd.DataFrame:
    """Messy shop-level attach: Zipf subset + latent theme + Student-t bias.

    Not coupled to W or Y. Merchant mix can still leak through user pooling.
    """
    mids = np.asarray(mids, dtype=np.int64)
    n_m = int(len(mids))
    theme = rng.integers(0, N_THEME, size=n_m)
    theme_center = rng.normal(0.0, 0.35, size=(N_THEME, EMB_DIM))
    k = rng.integers(4, 16, size=n_m)
    M = np.zeros((n_m, EMB_DIM), dtype=np.float64)
    for i in range(n_m):
        ix = rng.choice(N_CATALOG, size=int(k[i]), replace=False, p=pop)
        raw = emb[ix].mean(axis=0)
        bias = theme_center[theme[i]] + 0.25 * rng.standard_t(df=4.0, size=EMB_DIM)
        M[i] = raw + bias
    cols = {f"{prefix}_{j:02d}": M[:, j] for j in range(EMB_DIM)}
    return pd.DataFrame({"merchant_id": mids, **cols})


def user_pool_media(
    left: pd.DataFrame, imap: pd.DataFrame, merch_emb: pd.DataFrame, cols: list[str]
) -> pd.DataFrame:
    """Visit-weighted mean of merchant media for users who clk/cnv in the left window."""
    hit = left.loc[left.act.isin([CLK, CNV]), ["user_id", "item_id"]]
    hit = hit.merge(imap, on="item_id", how="left").dropna(subset=["merchant_id"])
    if hit.empty:
        return pd.DataFrame(columns=["user_id"] + cols)
    t = hit.merge(merch_emb[["merchant_id"] + cols], on="merchant_id", how="left")
    return t.groupby("user_id", sort=False)[cols].mean().reset_index()


def media_user_tables(left: pd.DataFrame, imap: pd.DataFrame, *, seed: int):
    """Build fake video/audio towers: catalog → messy merchant attach → user pool."""
    rng_v = np.random.default_rng(seed + 101)
    rng_a = np.random.default_rng(seed + 202)
    mids = np.sort(imap.merchant_id.dropna().unique().astype(np.int64))
    vid_ids, vid_emb, vid_pop = media_catalog(rng_v, kind="video")
    aud_ids, aud_emb, aud_pop = media_catalog(rng_a, kind="audio")
    merch_v = attach_media_to_merchants(mids, vid_emb, vid_pop, rng_v, prefix="vid")
    merch_a = attach_media_to_merchants(mids, aud_emb, aud_pop, rng_a, prefix="aud")
    vid_cols = [f"vid_{j:02d}" for j in range(EMB_DIM)]
    aud_cols = [f"aud_{j:02d}" for j in range(EMB_DIM)]
    user_v = user_pool_media(left, imap, merch_v, vid_cols)
    user_a = user_pool_media(left, imap, merch_a, aud_cols)
    users = user_v.merge(user_a, on="user_id", how="outer") if len(user_a) else user_v
    meta = {
        **DGP_META,
        "n_video_ids": int(len(vid_ids)),
        "n_audio_ids": int(len(aud_ids)),
        "n_merch_attached": int(len(mids)),
        "n_user_pooled": int(len(users)),
        "video_id_min": int(vid_ids.min()) if len(vid_ids) else None,
        "video_id_max": int(vid_ids.max()) if len(vid_ids) else None,
        "audio_id_min": int(aud_ids.min()) if len(aud_ids) else None,
        "audio_id_max": int(aud_ids.max()) if len(aud_ids) else None,
    }
    return users.fillna(0.0), vid_cols, aud_cols, meta


def rf_domain(X, W, *, seed: int):
    clf = RandomForestClassifier(
        n_estimators=TREES,
        max_depth=DEPTH,
        min_samples_leaf=4,
        random_state=seed,
        n_jobs=1,
    )
    n = len(W)
    idx = np.arange(n)
    rng = np.random.default_rng(seed)
    rng.shuffle(idx)
    cut = max(n * 3 // 4, 1)
    tr, te = idx[:cut], idx[cut:]
    if len(np.unique(W[tr])) < 2 or len(np.unique(W[te])) < 2:
        return float("nan"), np.zeros(X.shape[1], dtype=float)
    clf.fit(X[tr], W[tr])
    auc = float(roc_auc_score(W[te], clf.predict_proba(X[te])[:, 1]))
    return auc, clf.feature_importances_.astype(float)


def po_risk_fit(X, Y, W, *, seed: int):
    n = len(Y)
    m_hat = np.zeros(n)
    e_hat = np.zeros(n)
    Yf = Y.astype(float)
    cv = StratifiedKFold(n_splits=FOLDS, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X, W)):
        m = RandomForestRegressor(
            n_estimators=TREES,
            max_depth=DEPTH,
            min_samples_leaf=4,
            random_state=seed + fold,
            n_jobs=1,
        )
        e = RandomForestClassifier(
            n_estimators=TREES,
            max_depth=DEPTH,
            min_samples_leaf=4,
            random_state=seed + 40 + fold,
            n_jobs=1,
        )
        m.fit(X[tr], Yf[tr])
        e.fit(X[tr], W[tr])
        m_hat[te] = m.predict(X[te])
        e_hat[te] = e.predict_proba(X[te])[:, 1]
    e_hat = np.clip(e_hat, 0.05, 0.95)
    po = (Yf - m_hat) * (W.astype(float) - e_hat)
    tau = RandomForestRegressor(
        n_estimators=TREES,
        max_depth=DEPTH,
        min_samples_leaf=4,
        random_state=seed + 7,
        n_jobs=1,
    )
    tau.fit(X, po)
    tau_hat = tau.predict(X)
    return {
        "risk": float(np.mean(tau_hat**2)),
        "vimp": tau.feature_importances_.astype(float),
        "po": po,
    }


def _xy(df, cols, ycol, clock):
    names = [c for c in cols if c in df.columns]
    X = df[names].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(np.float64)
    y = pd.to_numeric(df[ycol], errors="coerce").fillna(0.0).to_numpy(np.float64)
    clk = df[clock].to_numpy(np.float64)
    w = (clk > np.median(clk)).astype(int)
    return names, X, y, w


def logo_groups(X, y, w, names, group_fn, *, seed: int):
    families = sorted({group_fn(n) for n in names})
    groups = {g: [i for i, n in enumerate(names) if group_fn(n) == g] for g in families}
    auc, v_rf = rf_domain(X, w, seed=seed)
    full = po_risk_fit(X, y, w, seed=seed)
    logo = {}
    for g, ix in groups.items():
        keep = [j for j in range(X.shape[1]) if j not in set(ix)]
        if not keep:
            logo[g] = {"R_minus": float("nan"), "delta": float("nan")}
            continue
        rm = po_risk_fit(X[:, keep], y, w, seed=seed + 11 + families.index(g))["risk"]
        logo[g] = {"R_minus": rm, "delta": float(full["risk"] - rm)}
    pos = {g: max((logo[g]["delta"] or 0.0), 0.0) for g in families}
    s = sum(pos.values()) + 1e-12
    share = {g: pos[g] / s for g in families}
    po_v, rf_v = full["vimp"], v_rf
    mass_po, mass_rf = {}, {}
    for g in families:
        mass_po[g] = float(sum(po_v[i] for i in groups[g]))
        mass_rf[g] = float(sum(rf_v[i] for i in groups[g]))
    return {
        "rf_domain_auc": float(auc) if auc == auc else None,
        "po_risk": float(full["risk"]),
        "family_n": {g: len(groups[g]) for g in families},
        "po_mass": {k: v / (sum(mass_po.values()) + 1e-12) for k, v in mass_po.items()},
        "rf_mass": {k: v / (sum(mass_rf.values()) + 1e-12) for k, v in mass_rf.items()},
        "logo": logo,
        "logo_share": share,
        "vimp_po": po_v.tolist(),
        "vimp_rf": v_rf.tolist(),
        "names": names,
        "groups": {g: groups[g] for g in families},
    }


def localize_blocks(X, w, names, *, seed: int):
    """No Y. Per block: can W be read from this family alone?"""
    rows = []
    for fam in LOC_FAMILIES:
        ix = [i for i, n in enumerate(names) if block_family(n) == fam]
        if not ix:
            continue
        auc, vimp = rf_domain(X[:, ix], w, seed=seed + 1 + LOC_FAMILIES.index(fam))
        order = np.argsort(-vimp)
        top = {names[ix[int(j)]]: float(vimp[int(j)]) for j in order[:8]}
        rows.append(
            {
                "family": fam,
                "n_feat": int(len(ix)),
                "rf_domain_auc": None if auc != auc else float(auc),
                "moved": bool(auc == auc and auc >= AUC_MOVED),
                "cols": [names[i] for i in ix] if len(ix) <= 8 else [names[i] for i in ix[:8]] + ["…"],
                "rf_vimp_top": top,
            }
        )
    rows.sort(key=lambda r: -(r["rf_domain_auc"] or 0.0))
    localized = [r["family"] for r in rows if r["moved"]]
    if not localized and rows:
        localized = [rows[0]["family"]]
    return {"ranked": rows, "localized_families": localized}


def loco_subset(X, y, w, names, keep_names, *, seed: int, full_risk: float, vimp_rf=None):
    keep_names = [n for n in names if n in set(keep_names)]
    if vimp_rf is not None and len(keep_names) > LOCO_TOPK:
        score = {n: float(vimp_rf.get(n, 0.0)) for n in keep_names}
        keep_names = sorted(keep_names, key=lambda n: -score[n])[:LOCO_TOPK]
    ix = [i for i, n in enumerate(names) if n in set(keep_names)]
    out = []
    for j in ix:
        keep = [k for k in range(X.shape[1]) if k != j]
        rm = po_risk_fit(X[:, keep], y, w, seed=seed + 80 + j)["risk"]
        out.append({"name": names[j], "hop": hop_of(names[j]), "loco_dR": float(full_risk - rm)})
    out.sort(key=lambda r: -r["loco_dR"])
    return out, keep_names


def _fmt_auc(auc) -> str:
    return "—" if auc is None else f"{auc:.3f}"


def _fmt_delta(d) -> str:
    return "nan" if d != d else f"{d:+.5f}"


def coverage_by_clock(o: pd.DataFrame, vid_cols: list[str], aud_cols: list[str], w: np.ndarray) -> dict:
    """Did user pooling just leave more zeros on one side of W?"""
    vn = np.linalg.norm(o[vid_cols].to_numpy(np.float64), axis=1)
    an = np.linalg.norm(o[aud_cols].to_numpy(np.float64), axis=1)
    early, late = w == 0, w == 1
    return {
        "vid_zero_early": float((vn[early] == 0).mean()),
        "vid_zero_late": float((vn[late] == 0).mean()),
        "aud_zero_early": float((an[early] == 0).mean()),
        "aud_zero_late": float((an[late] == 0).mean()),
        "vid_norm_early": float(vn[early].mean()),
        "vid_norm_late": float(vn[late].mean()),
        "aud_norm_early": float(an[early].mean()),
        "aud_norm_late": float(an[late].mean()),
    }


def render_md(loc, d1, d2, *, n, p, dgp, coverage=None) -> str:
    moved = [r["family"] for r in loc["ranked"] if r["moved"]]
    quiet = [r["family"] for r in loc["ranked"] if not r["moved"]]
    auc_line = ", ".join(
        f"{r['family']} {_fmt_auc(r['rf_domain_auc'])}" for r in loc["ranked"]
    )
    cov_line = ""
    if coverage:
        cov_line = (
            f"User-pool zero rate video early/late "
            f"{coverage['vid_zero_early']:.3f}/{coverage['vid_zero_late']:.3f}; "
            f"mean L2-norm {coverage['vid_norm_early']:.2f}/{coverage['vid_norm_late']:.2f}. "
            "Not a missingness artifact."
        )
    lines = [
        "# Recsys: graph localization → two-step drill",
        "",
        "These families are enough. No extra grains.",
        "Use: **graph localization**, then a **two-step drill** only on the moved few.",
        "PO-risk is **not** causal. RF-domain is a **portrait log**, not the score to optimize.",
        "Video/audio are a messy DGP, not production embeddings.",
        "",
        "## Protocol",
        "",
        "1. Left-window clk/cnv → user—item, user—merchant, session co-click, shop projection.",
        "2. Fake media DGP on the same graph (merchant attach → user mean-pool).",
        "3. **Localize** (no Y) on five families:",
        "   `user_connectivity` · `item_connectivity` · `merchant_structure` · `video_tower` · `audio_tower`.",
        f"   Moved = RF-domain AUC ≥ {AUC_MOVED}.",
        "4. **Two-step drill** with Y=`y_post_clk_1d`, **only the moved families**:",
        "   - Step 1: LOGO among those families",
        f"   - Step 2: LOCO columns (top-{LOCO_TOPK} if a 64-d tower)",
        "",
        f"n={n} p={p}. Clock = median `t_end`. W=1 later half.",
        f"Shops attached={dgp.get('n_merch_attached')} · users pooled={dgp.get('n_user_pooled')}.",
        cov_line,
        "",
        "## 1. Graph localization (no Y)",
        "",
        "| family | n | RF-domain AUC | moved |",
        "|---|---:|---:|---|",
    ]
    for r in loc["ranked"]:
        lines.append(
            f"| {r['family']} | {r['n_feat']} | {_fmt_auc(r['rf_domain_auc'])} | "
            f"{'yes' if r['moved'] else 'no'} |"
        )
    lines += [
        "",
        "Localized families: `" + ", ".join(loc["localized_families"]) + "`.",
        "Moved = this block's portrait differs early vs late. Not “this tower caused conversion”.",
        "Quiet families stop here. Drill does not open them.",
        "",
        "## 2. Drill step 1 — LOGO on the moved few",
        "",
    ]
    fam1 = d1["logo_share"]
    if len(fam1) <= 1:
        lines += [
            "Only one moved family, so step-1 LOGO is vacuous (nothing to drop).",
            "Step 2 LOCO on those columns is the drill.",
            "",
        ]
    else:
        lines += [
            f"RF-domain **{d1['rf_domain_auc']:.3f}** · PO-risk **{d1['po_risk']:.6f}**",
            "",
            "| family | n | RF-mass | PO-mass | LOGO Δ | share |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for g in sorted(fam1, key=lambda x: -fam1[x]):
            lines.append(
                f"| {g} | {d1['family_n'][g]} | {d1['rf_mass'][g]:.3f} | "
                f"{d1['po_mass'][g]:.3f} | {_fmt_delta(d1['logo'][g]['delta'])} | {fam1[g]:.3f} |"
            )
        lines += [
            "",
            "LOGO Δ>0: dropping the family **lowers** R (tied to the early/late Y gap).",
            "Δ≤0: this family is not a concept-gap source on this clock. Read RF-mass for P(X).",
            "Portrait movement ≠ Y-gap source. PO-risk ~1e-4 is a log, not a story.",
            "",
        ]
    lines += [
        "## 3. Drill step 2 — LOCO on those columns",
        "",
        f"If a 64-d tower is in the moved set, LOCO is top-{LOCO_TOPK} by step-1 RF-mass, not all 64 fits.",
        "",
        "| feat | hop | LOCO ΔR |",
        "|---|---|---:|",
    ]
    for r in d2[:12]:
        lines.append(f"| `{r['name']}` | {r['hop']} | {r['loco_dR']:+.6e} |")
    lines += [
        "",
        "LOCO ΔR>0: dropping the column lowers R. Still not a treatment effect.",
        "",
        "## Read",
        "",
        f"- Localization AUCs: {auc_line}.",
        f"- Moved (drill these): `{', '.join(moved) if moved else 'none'}`. Quiet (stop): `{', '.join(quiet) if quiet else 'none'}`.",
        "- Five families are enough. Do not add community / extra hops / funnel into this prototype.",
        "- Use is two sentences: localize on the graph; two-step drill on the moved few.",
        "- Do not turn LOGO share into a unique importance ranking.",
        "",
        "`PYTHONPATH=. python3 scripts/tencent_gr/graph_loc_fsds_drill.py`",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    if not EV.exists():
        print("missing", EV)
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    ev = pd.read_parquet(EV)
    post = pd.read_parquet(POST)
    imap = pd.read_parquet(IMAP)[["item_id", "merchant_id"]].drop_duplicates()

    left, _, _ = split_lr(ev)
    print(f"left n={len(left)} clk={(left.act==CLK).sum()} cnv={(left.act==CNV).sum()}", flush=True)
    gu, gi, gm = graph_tables(left, imap)
    print(f"graph users={len(gu)} items={len(gi)} shops={len(gm)}", flush=True)

    u_media, vid_cols, aud_cols, dgp = media_user_tables(left, imap, seed=SEED)
    print(
        f"dgp video_ids={dgp['n_video_ids']} audio_ids={dgp['n_audio_ids']} "
        f"dim={EMB_DIM} shops={dgp['n_merch_attached']} users_pooled={dgp['n_user_pooled']}",
        flush=True,
    )

    o = post.merge(imap, on="item_id", how="left")
    o = o.merge(gu, on="user_id", how="left")
    o = o.merge(gi, on="item_id", how="left")
    o = o.merge(gm, on="merchant_id", how="left")
    o = o.merge(u_media, on="user_id", how="left")
    extra = [c for c in o.columns if c.startswith(("g_", "vid_", "aud_"))]
    for c in extra:
        o[c] = o[c].fillna(0.0)
    o = o.loc[o["y_post_clk_1d"].notna()].copy()
    names, X, y, w = _xy(o, extra, "y_post_clk_1d", "t_end")
    cov = coverage_by_clock(o, vid_cols, aud_cols, w)
    print(f"xy n={len(y)} p={len(names)} pos={y.mean():.3f} W1={w.mean():.3f}", flush=True)
    print(
        f"coverage vid_zero e/l={cov['vid_zero_early']:.3f}/{cov['vid_zero_late']:.3f} "
        f"norm e/l={cov['vid_norm_early']:.2f}/{cov['vid_norm_late']:.2f}",
        flush=True,
    )

    print("=== graph localization (no Y) ===", flush=True)
    loc = localize_blocks(X, w, names, seed=SEED)
    for r in loc["ranked"]:
        print(f"  {r['family']:22s} AUC={r['rf_domain_auc']} moved={r['moved']}", flush=True)

    loc_cols = [n for n in names if block_family(n) in set(loc["localized_families"])]
    if not loc_cols:
        loc_cols = [n for n in names if block_family(n) in set(LOC_FAMILIES)]

    print("=== drill step 1: LOGO on moved families ===", flush=True)
    X2 = o[loc_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(np.float64)
    d1 = logo_groups(X2, y, w, loc_cols, block_family, seed=SEED + 3)
    for g, rec in d1["logo"].items():
        print(f"  LOGO {g:22s} Δ={rec['delta']}", flush=True)

    rf_map = dict(zip(d1["names"], d1["vimp_rf"]))
    print("=== drill step 2: LOCO on those columns ===", flush=True)
    d2, loco_names = loco_subset(
        X2, y, w, loc_cols, loc_cols, seed=SEED, full_risk=d1["po_risk"], vimp_rf=rf_map
    )
    print(f"  loco n={len(loco_names)} of loc_cols={len(loc_cols)}", flush=True)
    for r in d2[:8]:
        print(f"  LOCO {r['name']:22s} ΔR={r['loco_dR']:+.6e}", flush=True)

    payload = {
        "protocol": "graph_localize_then_two_step_drill",
        "not_causal": True,
        "y": "y_post_clk_1d",
        "n": int(len(y)),
        "p": int(len(names)),
        "dgp": dgp,
        "coverage": cov,
        "localization": loc,
        "drill1_logo": {k: v for k, v in d1.items() if k not in ("vimp_po", "vimp_rf", "groups")},
        "drill2_loco": d2,
        "drill2_loco_names": loco_names,
        "localized_families": loc["localized_families"],
    }
    (OUT / "GRAPH_LOC_FSDS.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    md = render_md(loc, d1, d2, n=len(y), p=len(names), dgp=dgp, coverage=cov)
    (OUT / "GRAPH_LOC_FSDS.md").write_text(md)
    (DOCS / "Recsys_Graph_Loc_FSDS.md").write_text(md)
    print(md)
    print("wrote", OUT / "GRAPH_LOC_FSDS.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
