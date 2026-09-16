"""Feature catalog for the serving-gate tables already on disk.

Every coordinate is an `x_*` column. Raw prompts, wiki text, gold flags,
and HH chosen/rejected are not features.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# (column, meaning). Order matches the csv header.
JUDGE_X = [
    ("x_n_toks", "token count (scaled)"),
    ("x_n_chars", "character count (scaled)"),
    ("x_avg_word", "mean word length"),
    ("x_qmark", "question marks"),
    ("x_bang", "exclamation marks"),
    ("x_hedge", "hedge phrases (maybe / I think / …)"),
    ("x_formal", "formal connectives (therefore / however / …)"),
    ("x_i_count", "first-person I"),
    ("x_newlines", "newline count"),
    ("x_upper", "uppercase fraction"),
    ("x_refuse", "refusal phrases (I can't / I'm sorry / …)"),
    ("x_please", "please"),
    ("x_thank", "thank"),
]

GRAPH_X = [
    ("x_n_nodes", "title nodes in the pool"),
    ("x_n_edges", "shared-token co-mention edges"),
    ("x_mean_deg", "mean degree"),
    ("x_n_cc", "connected components"),
    ("x_n_q_seeds", "titles overlapping the query"),
    ("x_seed_frac", "seed fraction"),
    ("x_lcc_frac", "largest-component fraction"),
]

GRAPH_BATCH_X = GRAPH_X + [(f"{c}_std", f"window std of {c}") for c, _ in GRAPH_X]

HYBRID_X = [
    ("x_q_toks", "query token count"),
    ("x_q_chars", "query character count"),
    ("x_qmark", "query has '?'"),
    ("x_q_ents", "capitalized query tokens"),
    ("x_n_cand", "candidate pool size (10)"),
    ("x_js_overlap", "Jaccard of sparse vs dense top-k titles"),
    ("x_rank_corr", "Spearman of sparse vs dense ranks"),
    ("x_sparse_margin", "BM25 top1 − top2"),
    ("x_dense_margin", "dense top1 − top2"),
    ("x_rrf_top1_mass", "RRF mass on fused top-1"),
    ("x_fuse_uniq", "unique titles in fused top-k"),
    ("x_mean_q_overlap", "mean query–doc token overlap in fused pack"),
    ("x_ks", "sparse k"),
    ("x_kd", "dense k"),
]

PAIR_X = [
    ("x_bm25", "BM25 z-score of this title"),
    ("x_dense", "char-ngram cosine z-score of this title"),
    ("x_rrf", "RRF z-score of this title"),
    ("x_rank_sp", "inverse sparse rank"),
    ("x_rank_de", "inverse dense rank"),
    ("x_q_overlap", "query–paragraph token overlap"),
    ("x_title_seed", "title tokens overlap the query"),
    ("x_title_deg", "co-mention degree"),
    ("x_slot", "slot index in the 10-title pool"),
]

BLOCKS = [
    {
        "id": "judge",
        "facet": "审核 / judge",
        "y": "pass / fail",
        "table": "results/manuscript/llm_audit/xy_hh_helpful_consistent.csv",
        "n_rows": 1200,
        "features": JUDGE_X,
        "note": "Same 13 X on BeaverTails / WildGuard / ToxicChat. HH chosen is not Y.",
        "gate": True,
    },
    {
        "id": "graph_query",
        "facet": "Graph-RAG 子图（一问一行）",
        "y": "supporting titles sit in the served pack",
        "table": "results/manuscript/graph_rag_batches/xy_graph_query.csv",
        "n_rows": 1200,
        "features": GRAPH_X,
        "note": "The probe uses these 7. Community hop rewires edges and serves the LCC.",
        "gate": True,
    },
    {
        "id": "graph_batch",
        "facet": "Graph-RAG 子图（一窗一行）",
        "y": "pack-usable rate on the window",
        "table": "results/manuscript/graph_rag_batches/xy_graph_batch.csv",
        "n_rows": 15,
        "features": GRAPH_BATCH_X,
        "note": "Serving-table shape: mean and std of the 7 graph coordinates. Not the probe sample.",
        "gate": False,
    },
    {
        "id": "hybrid",
        "facet": "混合检索",
        "y": "all gold titles in fused top-k",
        "table": "results/manuscript/hybrid_retrieval/xy_hotpot_hybrid.csv",
        "n_rows": 1200,
        "features": HYBRID_X,
        "note": "Channel-agreement geometry. Gold hits are not X. Hop flips the dense channel.",
        "gate": True,
    },
    {
        "id": "pair",
        "facet": "Graph-RAG 快照（一问 × 10 标题）",
        "y": "this title is a supporting fact",
        "table": "results/manuscript/hybrid_retrieval/xy_hotpot_pairs.csv",
        "n_rows": 12000,
        "features": PAIR_X,
        "note": "Native Hotpot pool shape. batch = query index, not wall-clock time.",
        "gate": False,
    },
]


def cols(block_id: str) -> list[str]:
    for b in BLOCKS:
        if b["id"] == block_id:
            return [c for c, _ in b["features"]]
    raise KeyError(block_id)


def disk_x_cols(path: Path) -> list[str]:
    import csv

    with path.open() as f:
        header = next(csv.reader(f))
    return [c for c in header if c.startswith("x_")]
