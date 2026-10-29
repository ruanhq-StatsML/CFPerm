"""Shift boards for Magnitude × Variance eval — one loader, many boards.

Each board is binary text (y∈{0,1}) with enough REF/LIVE mass to build
benign (uniform synonym) vs malignant (subset attack) shifts.

Boards
------
  hc3_all / hc3_finance / hc3_medicine — AIGC integrity (HC3)
  news20_sci / news20_talk             — 20 Newsgroups topic pairs (sklearn)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
from sklearn.datasets import fetch_20newsgroups
from sklearn.model_selection import train_test_split

from as_dta_aigc import HC3_INDUSTRIAL, SEED, load_hc3_binary

BOARD_NAMES = (
    "hc3_all",
    "hc3_finance",
    "hc3_medicine",
    "news20_sci",
    "news20_talk",
)


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


def load_hc3_board(name: str, max_pairs: int = 1200, seed: int = SEED) -> ShiftBoard:
    key = {"hc3_all": "all", "hc3_finance": "finance", "hc3_medicine": "medicine"}[name]
    path = HC3_INDUSTRIAL[key]
    data = load_hc3_binary(path, max_pairs=max_pairs, seed=seed)
    return ShiftBoard(
        name=name,
        texts=np.asarray(data["texts"], dtype=object),
        labels=np.asarray(data["labels"], dtype=int),
        domain="aigc",
        meta={"source": "HC3", "vertical": key, "path": str(path)},
    )


def load_news20_board(
    name: str,
    max_per_class: int = 600,
    seed: int = SEED,
    max_chars: int = 800,
) -> ShiftBoard:
    """Binary topic pair from 20 Newsgroups (local sklearn)."""
    pairs = {
        "news20_sci": ("sci.med", "sci.space"),
        "news20_talk": ("talk.politics.guns", "talk.politics.mideast"),
    }
    c0, c1 = pairs[name]
    bun = fetch_20newsgroups(
        subset="all",
        categories=[c0, c1],
        remove=("headers", "footers", "quotes"),
        shuffle=True,
        random_state=int(seed),
    )
    texts, labels = [], []
    rng = np.random.default_rng(int(seed))
    by = {0: [], 1: []}
    for t, y in zip(bun.data, bun.target):
        by[int(y)].append(t)
    for y in (0, 1):
        idx = rng.choice(len(by[y]), size=min(max_per_class, len(by[y])), replace=False)
        for i in idx:
            texts.append(by[y][int(i)])
            labels.append(y)
    texts = _clip_texts(texts, max_chars=max_chars)
    labels = np.asarray(labels[: len(texts)], dtype=int)
    # realign if clip dropped empties — re-sample already non-empty mostly
    n = min(len(texts), len(labels))
    return ShiftBoard(
        name=name,
        texts=texts[:n],
        labels=labels[:n],
        domain="topic",
        meta={"source": "20newsgroups", "classes": [c0, c1]},
    )


_LOADERS = {
    "hc3_all": load_hc3_board,
    "hc3_finance": load_hc3_board,
    "hc3_medicine": load_hc3_board,
    "news20_sci": load_news20_board,
    "news20_talk": load_news20_board,
}


def load_board(name: str, **kw) -> ShiftBoard:
    if name not in _LOADERS:
        raise KeyError("unknown board %r; choose from %s" % (name, BOARD_NAMES))
    return _LOADERS[name](name, **kw)


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
