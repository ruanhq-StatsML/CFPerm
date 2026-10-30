"""Shift boards for Magnitude × Variance eval — one loader, many boards.

Each board is binary text (y∈{0,1}) with enough REF/LIVE mass to build
benign (uniform synonym) vs malignant (subset attack + Y-flip) shifts.

Boards
------
  hc3_all / hc3_finance / hc3_medicine / hc3_open_qa / hc3_wiki
      — AIGC integrity (HC3 verticals; open_qa / wiki from all.jsonl source)
  news20_sci / news20_talk / news20_comp / news20_rec
      — 20 Newsgroups topic pairs (sklearn)
  ag_news — AG News World vs Sci/Tech (Hugging Face fancyzhx/ag_news)
  imdb    — IMDB sentiment (Hugging Face stanfordnlp/imdb)
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.datasets import fetch_20newsgroups
from sklearn.model_selection import train_test_split

from as_dta_aigc import HC3_DEFAULT, HC3_INDUSTRIAL, ROOT, SEED, load_hc3_binary

BOARD_NAMES = (
    "hc3_all",
    "hc3_finance",
    "hc3_medicine",
    "hc3_open_qa",
    "hc3_wiki",
    "news20_sci",
    "news20_talk",
    "news20_comp",
    "news20_rec",
    "ag_news",
    "imdb",
)

# Concrete multi-dataset justify set (diverse families, same口径)
JUSTIFY_DATASETS = (
    "hc3_finance",
    "hc3_medicine",
    "hc3_open_qa",
    "hc3_wiki",
    "news20_sci",
    "news20_talk",
    "news20_comp",
    "ag_news",
    "imdb",
)

NEWS20_PAIRS = {
    "news20_sci": ("sci.med", "sci.space"),
    "news20_talk": ("talk.politics.guns", "talk.politics.mideast"),
    "news20_comp": ("comp.sys.ibm.pc.hardware", "comp.sys.mac.hardware"),
    "news20_rec": ("rec.autos", "rec.sport.hockey"),
}

HC3_SOURCE_BOARDS = {
    "hc3_open_qa": "open_qa",
    "hc3_wiki": "wiki_csai",
}


@dataclass
class ShiftBoard:
    name: str
    texts: np.ndarray
    labels: np.ndarray
    domain: str = "text"
    meta: dict = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(len(self.labels))


def _clip_texts(texts: Sequence[str], max_chars: int = 800) -> np.ndarray:
    out = []
    for t in texts:
        s = " ".join(str(t).split())
        if not s:
            continue
        out.append(s[:max_chars])
    return np.asarray(out, dtype=object)


def _balance_binary(
    texts: Sequence[str],
    labels: Sequence[int],
    max_per_class: int,
    seed: int,
    max_chars: int = 800,
) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(int(seed))
    by: Dict[int, List[str]] = {0: [], 1: []}
    for t, y in zip(texts, labels):
        y = int(y)
        if y not in by:
            continue
        s = " ".join(str(t).split())
        if not s:
            continue
        by[y].append(s[:max_chars])
    out_t, out_y = [], []
    for y in (0, 1):
        pool = by[y]
        if not pool:
            raise ValueError("class %d empty after clip" % y)
        take = min(int(max_per_class), len(pool))
        idx = rng.choice(len(pool), size=take, replace=False)
        for i in idx:
            out_t.append(pool[int(i)])
            out_y.append(y)
    order = rng.permutation(len(out_t))
    return (
        np.asarray([out_t[i] for i in order], dtype=object),
        np.asarray([out_y[i] for i in order], dtype=int),
    )


def load_hc3_board(name: str, max_pairs: int = 1200, seed: int = SEED) -> ShiftBoard:
    key = {
        "hc3_all": "all",
        "hc3_finance": "finance",
        "hc3_medicine": "medicine",
    }[name]
    path = HC3_INDUSTRIAL[key]
    data = load_hc3_binary(path, max_pairs=max_pairs, seed=seed)
    return ShiftBoard(
        name=name,
        texts=np.asarray(data["texts"], dtype=object),
        labels=np.asarray(data["labels"], dtype=int),
        domain="aigc",
        meta={"source": "HC3", "vertical": key, "path": str(path)},
    )


def load_hc3_source_board(
    name: str,
    max_pairs: int = 1200,
    seed: int = SEED,
    max_chars: int = 800,
) -> ShiftBoard:
    """HC3 verticals filtered from all.jsonl by ``source`` (open_qa / wiki_csai)."""
    source = HC3_SOURCE_BOARDS[name]
    path = Path(HC3_DEFAULT)
    if not path.exists():
        from as_dta_aigc import ensure_hc3

        ensure_hc3(path)
    rng = np.random.default_rng(int(seed))
    humans, machines = [], []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            if str(row.get("source", "")) != source:
                continue
            for ans in row.get("human_answers") or []:
                if isinstance(ans, str) and ans.strip():
                    humans.append(ans.strip()[:max_chars])
            for ans in row.get("chatgpt_answers") or []:
                if isinstance(ans, str) and ans.strip():
                    machines.append(ans.strip()[:max_chars])
    rng.shuffle(humans)
    rng.shuffle(machines)
    n = min(int(max_pairs), len(humans), len(machines))
    if n < 40:
        raise RuntimeError(
            "HC3 source %r too small (n=%d); need more rows in %s" % (source, n, path)
        )
    texts = humans[:n] + machines[:n]
    labels = [0] * n + [1] * n
    order = rng.permutation(len(texts))
    return ShiftBoard(
        name=name,
        texts=np.asarray([texts[i] for i in order], dtype=object),
        labels=np.asarray([labels[i] for i in order], dtype=int),
        domain="aigc",
        meta={"source": "HC3", "vertical": source, "path": str(path), "n_pairs": int(n)},
    )


def load_news20_board(
    name: str,
    max_per_class: int = 600,
    seed: int = SEED,
    max_chars: int = 800,
    max_pairs: Optional[int] = None,
) -> ShiftBoard:
    """Binary topic pair from 20 Newsgroups (local sklearn)."""
    if max_pairs is not None:
        max_per_class = max(40, int(max_pairs) // 2)
    c0, c1 = NEWS20_PAIRS[name]
    bun = fetch_20newsgroups(
        subset="all",
        categories=[c0, c1],
        remove=("headers", "footers", "quotes"),
        shuffle=True,
        random_state=int(seed),
    )
    texts, labels = _balance_binary(
        bun.data, bun.target, max_per_class=max_per_class, seed=seed, max_chars=max_chars
    )
    return ShiftBoard(
        name=name,
        texts=texts,
        labels=labels,
        domain="topic",
        meta={"source": "20newsgroups", "classes": [c0, c1]},
    )


def load_ag_news_board(
    name: str = "ag_news",
    max_per_class: int = 600,
    seed: int = SEED,
    max_chars: int = 800,
    max_pairs: Optional[int] = None,
) -> ShiftBoard:
    """AG News binary: World (0) vs Sci/Tech (3)."""
    if max_pairs is not None:
        max_per_class = max(40, int(max_pairs) // 2)
    from datasets import load_dataset

    # take a slice large enough to balance; full train is 120k
    ds = load_dataset("fancyzhx/ag_news", split="train[:20000]")
    texts, labels = [], []
    for row in ds:
        lab = int(row["label"])
        if lab not in (0, 3):
            continue
        texts.append(row["text"])
        labels.append(0 if lab == 0 else 1)
    texts, labels = _balance_binary(
        texts, labels, max_per_class=max_per_class, seed=seed, max_chars=max_chars
    )
    return ShiftBoard(
        name=name,
        texts=texts,
        labels=labels,
        domain="news",
        meta={
            "source": "fancyzhx/ag_news",
            "classes": ["World", "Sci/Tech"],
            "hf": "fancyzhx/ag_news",
        },
    )


def load_imdb_board(
    name: str = "imdb",
    max_per_class: int = 600,
    seed: int = SEED,
    max_chars: int = 800,
    max_pairs: Optional[int] = None,
) -> ShiftBoard:
    """IMDB binary sentiment."""
    if max_pairs is not None:
        max_per_class = max(40, int(max_pairs) // 2)
    from datasets import load_dataset

    # train is sorted by label (neg then pos); take both ends + shuffle via balance.
    neg = load_dataset("stanfordnlp/imdb", split="train[:4000]")
    pos = load_dataset("stanfordnlp/imdb", split="train[12500:16500]")
    texts = list(neg["text"]) + list(pos["text"])
    labels = list(neg["label"]) + list(pos["label"])
    texts, labels = _balance_binary(
        texts,
        labels,
        max_per_class=max_per_class,
        seed=seed,
        max_chars=max_chars,
    )
    return ShiftBoard(
        name=name,
        texts=texts,
        labels=labels,
        domain="sentiment",
        meta={"source": "stanfordnlp/imdb", "classes": ["neg", "pos"], "hf": "stanfordnlp/imdb"},
    )


_LOADERS = {
    "hc3_all": load_hc3_board,
    "hc3_finance": load_hc3_board,
    "hc3_medicine": load_hc3_board,
    "hc3_open_qa": load_hc3_source_board,
    "hc3_wiki": load_hc3_source_board,
    "news20_sci": load_news20_board,
    "news20_talk": load_news20_board,
    "news20_comp": load_news20_board,
    "news20_rec": load_news20_board,
    "ag_news": load_ag_news_board,
    "imdb": load_imdb_board,
}


def load_board(name: str, **kw) -> ShiftBoard:
    if name not in _LOADERS:
        raise KeyError("unknown board %r; choose from %s" % (name, BOARD_NAMES))
    # Drop kwargs the chosen loader does not accept (callers pass mixed bags).
    import inspect

    fn = _LOADERS[name]
    sig = inspect.signature(fn)
    accepted = {k: v for k, v in kw.items() if k in sig.parameters}
    return fn(name, **accepted)


def load_boards(names: Optional[Sequence[str]] = None, **kw) -> Dict[str, ShiftBoard]:
    names = list(names or BOARD_NAMES)
    return {n: load_board(n, **kw) for n in names}


def split_ref_live_board(
    board: ShiftBoard,
    live_frac: float = 0.35,
    seed: int = SEED,
):
    """Train REF vs LIVE holdout (stratified)."""
    X_ref, X_live, y_ref, y_live = train_test_split(
        board.texts,
        board.labels,
        test_size=float(live_frac),
        random_state=int(seed),
        stratify=board.labels,
    )
    return {
        "ref_texts": np.asarray(X_ref, dtype=object),
        "ref_labels": np.asarray(y_ref, dtype=int),
        "live_texts": np.asarray(X_live, dtype=object),
        "live_labels": np.asarray(y_live, dtype=int),
    }
