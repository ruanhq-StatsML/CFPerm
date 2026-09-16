#!/usr/bin/env python3
"""Formulate hybrid-retrieval and Graph-RAG prediction tables from HotpotQA.

Native shape is one query × 10 wiki paras (Hotpot distractor pool), not one
row per question. Pair table: Y_j = 1 iff title j is a supporting fact.
Collapsed fused Y is an optional readout of that ranking.

Questions and paragraph text are not stored.

Usage::

    PYTHONPATH=. python3 scripts/prototype_hotpot_10para_shape.py
    python3 scripts/build_hybrid_retrieval_xy.py
"""

Hybrid serving (what production logs):
  sparse = BM25 on the distractor candidate pool
  dense  = char-ngram hashing cosine (a second channel, not a GPU embedder)
  fuse   = Reciprocal Rank Fusion
  Y      = 1 iff every gold supporting *title* is in the fused top-k
  X      = query geometry + channel agreement. Gold hits are not features.

Graph-RAG (same queries, extra X):
  nodes = candidate titles
  edges = titles that share a token (cheap co-mention graph)
  seeds = titles whose tokens overlap the query
  same Y as fused (subgraph pack is usable iff supporting titles made the pack)

Two-stream CFPerm: T=0 sparse-only usable, T=1 dense-only usable, same X.

Usage::

    python3 scripts/build_hybrid_retrieval_xy.py
"""
from __future__ import annotations

import csv
import json
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.metrics.pairwise import cosine_similarity

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "hf_cache" / "retrieval"
OUT = ROOT / "results" / "manuscript" / "hybrid_retrieval"

PARQUET_URL = (
    "https://huggingface.co/datasets/hotpotqa/hotpot_qa/resolve/"
    "refs%2Fconvert%2Fparquet/distractor/validation/0000.parquet"
)
UA = {"User-Agent": "cfperm-hybrid-xy"}
TOK = re.compile(r"[a-z0-9]+")
K_SPARSE = 5
K_DENSE = 5
K_FUSE = 5
RRF_K = 60
N_ROWS = 1200
N_PER = 80
N_CAND = 10  # Hotpot distractor pool: one query × 10 wiki paras
CUT_BATCH = 4
PAIR_DIMS = [
    "bm25",
    "dense",
    "rrf",
    "rank_sp",
    "rank_de",
    "q_overlap",
    "title_seed",
    "title_deg",
    "slot",
]
PAIR_COLS = [f"x_{d}" for d in PAIR_DIMS]

HYBRID_DIMS = [
    "q_toks",
    "q_chars",
    "qmark",
    "q_ents",
    "n_cand",
    "js_overlap",
    "rank_corr",
    "sparse_margin",
    "dense_margin",
    "rrf_top1_mass",
    "fuse_uniq",
    "mean_q_overlap",
    "ks",
    "kd",
]
GRAPH_DIMS = [
    "n_nodes",
    "n_edges",
    "mean_deg",
    "n_cc",
    "n_q_seeds",
    "seed_frac",
    "lcc_frac",
]
HYBRID_COLS = [f"x_{d}" for d in HYBRID_DIMS]
GRAPH_COLS = [f"x_{d}" for d in GRAPH_DIMS]


def tokenize(text: str) -> list[str]:
    return TOK.findall((text or "").lower())


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 10_000:
        return dest
    print(f"download {url}")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as f:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
    return dest


def bm25_scores(query: str, docs: list[str], k1: float = 1.2, b: float = 0.75) -> np.ndarray:
    q = tokenize(query)
    doc_toks = [tokenize(d) for d in docs]
    n = len(docs)
    if n == 0:
        return np.zeros(0)
    avgdl = float(np.mean([max(len(d), 1) for d in doc_toks]))
    df: Counter[str] = Counter()
    for d in doc_toks:
        df.update(set(d))
    qtf = Counter(q)
    scores = np.zeros(n, dtype=float)
    for i, d in enumerate(doc_toks):
        dl = max(len(d), 1)
        tf = Counter(d)
        s = 0.0
        for w, qf in qtf.items():
            n_w = df.get(w, 0)
            idf = math.log(1.0 + (n - n_w + 0.5) / (n_w + 0.5))
            freq = tf.get(w, 0)
            s += qf * idf * (freq * (k1 + 1.0)) / (freq + k1 * (1.0 - b + b * dl / avgdl))
        scores[i] = s
    return scores


_DENSE_VEC = HashingVectorizer(
    analyzer="char", ngram_range=(3, 5), n_features=2**18, alternate_sign=False, norm="l2"
)


def dense_scores(query: str, docs: list[str]) -> np.ndarray:
    if not docs:
        return np.zeros(0)
    mat = _DENSE_VEC.transform([query] + docs)
    return cosine_similarity(mat[0], mat[1:]).ravel()


def rrf_from_scores(scores: np.ndarray, k: int = RRF_K) -> np.ndarray:
    order = np.argsort(-scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    return 1.0 / (k + ranks)


def top_idx(scores: np.ndarray, k: int) -> np.ndarray:
    k = min(int(k), len(scores))
    if k <= 0:
        return np.zeros(0, dtype=int)
    return np.argsort(-scores, kind="mergesort")[:k]


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2:
        return 0.0
    ra = np.argsort(np.argsort(-a))
    rb = np.argsort(np.argsort(-b))
    if float(ra.std()) < 1e-12 or float(rb.std()) < 1e-12:
        return 0.0
    return float(np.corrcoef(ra.astype(float), rb.astype(float))[0, 1])


def margin(scores: np.ndarray) -> float:
    if len(scores) < 2:
        return 0.0
    s = np.sort(scores)[::-1]
    return float(s[0] - s[1])


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / max(len(a | b), 1)


def title_graph(titles: list[str]):
    tok_sets = [set(tokenize(t)) for t in titles]
    n = len(titles)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    edges = 0
    deg = np.zeros(n, dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if tok_sets[i] & tok_sets[j]:
                edges += 1
                deg[i] += 1
                deg[j] += 1
                pi, pj = find(i), find(j)
                if pi != pj:
                    parent[pj] = pi
    roots = [find(i) for i in range(n)]
    n_cc = len(set(roots)) if n else 0
    lcc = max(Counter(roots).values()) if roots else 0
    return {
        "n_nodes": float(n),
        "n_edges": float(edges),
        "mean_deg": float(deg.mean()) if n else 0.0,
        "n_cc": float(n_cc),
        "lcc_frac": float(lcc / n) if n else 0.0,
        "deg": deg,
    }


def pad_pool(titles, docs, n: int = N_CAND):
    """Rectangular pool: (n_cand,) titles/docs. Truncate or pad empties."""
    titles = [str(t) for t in list(titles)[:n]]
    docs = [str(d) for d in list(docs)[:n]]
    while len(titles) < n:
        titles.append("")
        docs.append("")
    return titles, docs


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    s = float(x.std())
    if s < 1e-12:
        return np.zeros_like(x)
    return (x - float(x.mean())) / s


def inv_rank(scores: np.ndarray) -> np.ndarray:
    order = np.argsort(-scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    return 1.0 / ranks


def query_pool_tensors(question, context, supporting, *, fracture_dense: bool = False):
    """The actual serving tensor: one query × 10 wiki paras.

    Returns arrays all length ``N_CAND``::

        titles (10,)   docs (10,)
        sparse (10,)   dense (10,)   rrf (10,)
        y_pair (10,)   1 iff this title is a supporting fact
    """
    titles, docs = pad_pool(
        context["title"],
        [" ".join(s) for s in context["sentences"]],
        N_CAND,
    )
    gold = set(supporting["title"])
    sp = bm25_scores(question, docs)
    de = dense_scores(question, docs)
    if fracture_dense:
        de = -de
    rrf = rrf_from_scores(sp) + rrf_from_scores(de)
    y = np.asarray([1 if t in gold and t else 0 for t in titles], dtype=int)
    return {
        "titles": titles,
        "docs": docs,
        "sparse": sp,
        "dense": de,
        "rrf": rrf,
        "y_pair": y,
        "n_cand": N_CAND,
        "n_gold": int(y.sum()),
    }


def pair_rows_for_query(question, context, supporting, *, batch: int, fracture_dense: bool = False):
    """10 rows: slot j is wiki para j. Y = supporting title, not fused Recall."""
    ten = query_pool_tensors(question, context, supporting, fracture_dense=fracture_dense)
    titles = ten["titles"]
    docs = ten["docs"]
    sp, de, rrf = ten["sparse"], ten["dense"], ten["rrf"]
    seeds = query_seeds(question, titles)
    deg = title_graph(titles)["deg"]
    qset = set(tokenize(question))
    bm25_z, dense_z, rrf_z = zscore(sp), zscore(de), zscore(rrf)
    rsp, rde = inv_rank(sp), inv_rank(de)
    rows = []
    for j in range(N_CAND):
        dt = set(tokenize(docs[j]))
        ov = len(qset & dt) / max(len(qset | dt), 1) if docs[j] else 0.0
        rows.append(
            {
                "y": int(ten["y_pair"][j]),
                "batch": int(batch),
                "x_bm25": float(bm25_z[j]),
                "x_dense": float(dense_z[j]),
                "x_rrf": float(rrf_z[j]),
                "x_rank_sp": float(rsp[j]),
                "x_rank_de": float(rde[j]),
                "x_q_overlap": float(ov),
                "x_title_seed": float(seeds[j]),
                "x_title_deg": float(deg[j]) / 6.0,
                "x_slot": (j + 1) / N_CAND,
            }
        )
    return rows, ten


def query_seeds(query: str, titles: list[str]) -> np.ndarray:
    q = set(tokenize(query))
    return np.asarray([1.0 if (q & set(tokenize(t))) else 0.0 for t in titles], dtype=float)


def hybrid_x(
    query: str,
    titles: list[str],
    docs: list[str],
    sp: np.ndarray,
    de: np.ndarray,
    fused: np.ndarray,
) -> dict:
    q_toks = tokenize(query)
    sp_top = {titles[i] for i in top_idx(sp, K_SPARSE)}
    de_top = {titles[i] for i in top_idx(de, K_DENSE)}
    fuse_i = top_idx(fused, K_FUSE)
    fuse_docs = [docs[i] for i in fuse_i]
    ov = []
    qset = set(q_toks)
    for d in fuse_docs:
        dt = set(tokenize(d))
        ov.append(len(qset & dt) / max(len(qset | dt), 1))
    rrf = fused
    mass = float(rrf[fuse_i[0]] / max(rrf.sum(), 1e-12)) if len(fuse_i) else 0.0
    ents = len(re.findall(r"\b[A-Z][a-zA-Z]+\b", query or ""))
    return {
        "x_q_toks": len(q_toks) / 40.0,
        "x_q_chars": len(query or "") / 160.0,
        "x_qmark": float((query or "").count("?") > 0),
        "x_q_ents": ents / 6.0,
        "x_n_cand": len(docs) / 10.0,
        "x_js_overlap": jaccard(sp_top, de_top),
        "x_rank_corr": spearman(sp, de),
        "x_sparse_margin": margin(sp),
        "x_dense_margin": margin(de),
        "x_rrf_top1_mass": mass,
        "x_fuse_uniq": len({titles[i] for i in fuse_i}) / max(K_FUSE, 1),
        "x_mean_q_overlap": float(np.mean(ov)) if ov else 0.0,
        "x_ks": float(K_SPARSE),
        "x_kd": float(K_DENSE),
    }


def graph_x(query: str, titles: list[str]) -> dict:
    g = title_graph(titles)
    seeds = query_seeds(query, titles)
    n = max(len(titles), 1)
    return {
        "x_n_nodes": g["n_nodes"] / 10.0,
        "x_n_edges": g["n_edges"] / 20.0,
        "x_mean_deg": g["mean_deg"] / 6.0,
        "x_n_cc": g["n_cc"] / 10.0,
        "x_n_q_seeds": float(seeds.sum()) / 6.0,
        "x_seed_frac": float(seeds.mean()) if len(seeds) else 0.0,
        "x_lcc_frac": g["lcc_frac"],
    }


def gold_in_topk(titles: list[str], scores: np.ndarray, gold: set[str], k: int) -> int:
    hit = {titles[i] for i in top_idx(scores, k)}
    return int(bool(gold) and gold <= hit)


def pack_rows(rows: list[dict], n: int = N_ROWS, n_per: int = N_PER):
    n_use = (min(len(rows), n) // n_per) * n_per
    rows = rows[:n_use]
    batch = np.repeat(np.arange(n_use // n_per), n_per)
    for i, r in enumerate(rows):
        r["batch"] = int(batch[i])
    return rows


def write_xy(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["y", "batch", *cols])
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "y": int(r["y"]),
                    "batch": int(r["batch"]),
                    **{c: f"{float(r[c]):.6g}" for c in cols},
                }
            )


def write_two_stream(path: Path, sparse_rows, dense_rows, cols: list[str]) -> None:
    out = []
    for t, block, name in ((0, sparse_rows, "sparse"), (1, dense_rows, "dense")):
        for r in block:
            out.append(
                {
                    "T": t,
                    "Y": int(r["y"]),
                    **{c: r[c] for c in cols},
                    "batch": int(r["batch"]),
                    "stream": name,
                }
            )
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["T", "Y", *cols, "batch", "stream"])
        w.writeheader()
        for r in out:
            w.writerow(
                {
                    "T": r["T"],
                    "Y": r["Y"],
                    "batch": r["batch"],
                    "stream": r["stream"],
                    **{c: f"{float(r[c]):.6g}" for c in cols},
                }
            )


def formulate_example(question, context, supporting, *, fracture_dense: bool):
    titles = list(context["title"])
    docs = [" ".join(s) for s in context["sentences"]]
    gold = set(supporting["title"])
    sp = bm25_scores(question, docs)
    de = dense_scores(question, docs)
    if fracture_dense:
        de = -de
    fused = rrf_from_scores(sp) + rrf_from_scores(de)
    hx = hybrid_x(question, titles, docs, sp, de, fused)
    gx = graph_x(question, titles)
    return {
        **hx,
        **gx,
        "y_fused": gold_in_topk(titles, fused, gold, K_FUSE),
        "y_sparse": gold_in_topk(titles, sp, gold, K_SPARSE),
        "y_dense": gold_in_topk(titles, de, gold, K_DENSE),
    }


def load_hotpot(n: int):
    path = download(PARQUET_URL, CACHE / "hotpot_distractor_val.parquet")
    tbl = pq.read_table(path, columns=["question", "supporting_facts", "context"])
    return (
        tbl.column("question").to_pylist()[:n],
        tbl.column("supporting_facts").to_pylist()[:n],
        tbl.column("context").to_pylist()[:n],
    )


def build_regime(questions, supporting, contexts, *, hop: bool):
    fused, graph, sparse, dense = [], [], [], []
    for i, (q, sf, ctx) in enumerate(zip(questions, supporting, contexts)):
        batch_i = i // N_PER
        rec = formulate_example(q, ctx, sf, fracture_dense=bool(hop and batch_i >= CUT_BATCH))
        base = {c: rec[c] for c in HYBRID_COLS}
        gcols = {c: rec[c] for c in GRAPH_COLS}
        fused.append({**base, "y": rec["y_fused"]})
        graph.append({**base, **gcols, "y": rec["y_fused"]})
        sparse.append({**base, "y": rec["y_sparse"]})
        dense.append({**base, "y": rec["y_dense"]})
    return {
        "fused": pack_rows(fused),
        "graph": pack_rows(graph),
        "sparse": pack_rows(sparse),
        "dense": pack_rows(dense),
    }


def summarize(name, rows, ykey="y"):
    y = np.asarray([int(r[ykey]) for r in rows], float)
    return {
        "name": name,
        "n": len(rows),
        "y_pass_rate": float(y.mean()) if len(y) else None,
        "n_batches": int(max(r["batch"] for r in rows) + 1) if rows else 0,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("load HotpotQA distractor validation")
    q, sf, ctx = load_hotpot(N_ROWS)
    native = build_regime(q, sf, ctx, hop=False)
    hop = build_regime(q, sf, ctx, hop=True)
    pairs, pairs_hop = [], []
    for i, (qi, sfi, ctxi) in enumerate(zip(q, sf, ctx)):
        pr, _ = pair_rows_for_query(qi, ctxi, sfi, batch=i, fracture_dense=False)
        ph, _ = pair_rows_for_query(
            qi, ctxi, sfi, batch=i, fracture_dense=bool(i // N_PER >= CUT_BATCH)
        )
        pairs.extend(pr)
        pairs_hop.extend(ph)

    write_xy(OUT / "xy_hotpot_hybrid.csv", native["fused"], HYBRID_COLS)
    write_xy(OUT / "xy_hotpot_hybrid_hop.csv", hop["fused"], HYBRID_COLS)
    write_xy(OUT / "xy_hotpot_graph.csv", native["graph"], HYBRID_COLS + GRAPH_COLS)
    write_xy(OUT / "xy_hotpot_pairs.csv", pairs, PAIR_COLS)
    write_xy(OUT / "xy_hotpot_pairs_hop.csv", pairs_hop, PAIR_COLS)
    write_two_stream(
        OUT / "xy_hotpot_two_stream.csv",
        native["sparse"],
        native["dense"],
        HYBRID_COLS,
    )

    manifest = {
        "note": (
            "Hybrid retrieval / Graph-RAG prediction tables. "
            "Y is fused (or channel) gold-title coverage. Questions are not stored. "
            "Gold hit indicators are not X."
        ),
        "source": "hotpotqa/hotpot_qa distractor validation",
        "k_sparse": K_SPARSE,
        "k_dense": K_DENSE,
        "k_fuse": K_FUSE,
        "dense_channel": "char-ngram HashingVectorizer cosine (not a GPU embedding)",
        "sparse_channel": "BM25 on the Hotpot distractor candidate pool",
        "y_fused": "1 iff all supporting titles are in RRF fused top-k",
        "hop": f"dense scores flipped after batch>={CUT_BATCH} (embedding-pack swap)",
        "hybrid_x": HYBRID_DIMS,
        "graph_x": GRAPH_DIMS,
        "n_cand": N_CAND,
        "pair_x": PAIR_DIMS,
        "y_pair": "1 iff this of the 10 wiki titles is a supporting fact",
        "files": {
            "xy_hotpot_hybrid.csv": summarize("fused hybrid native", native["fused"]),
            "xy_hotpot_hybrid_hop.csv": summarize("fused hybrid dense-fracture", hop["fused"]),
            "xy_hotpot_graph.csv": summarize("fused Y + title-graph X", native["graph"]),
            "xy_hotpot_pairs.csv": summarize("query×10 para pairs", pairs),
            "xy_hotpot_pairs_hop.csv": summarize("pairs, dense flipped after cut queries", pairs_hop),
            "xy_hotpot_two_stream.csv": {
                "n": len(native["sparse"]) + len(native["dense"]),
                "T": "0=sparse-only usable, 1=dense-only usable",
                "X": "same hybrid geometry",
            },
        },
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["files"], indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
