#!/usr/bin/env python3
"""Build manuscript use-case 4 prediction tables: LLM moderation.

One binary Y (pass=1 / fail=0) and a pile of X (x_n_toks … x_thank).
HH chosen/rejected is never Y — those only mark which reviewer queue
the traffic came from (helpful vs harmless).

Sources
-------
- Anthropic HH-RLHF helpful/harmless: existing prototype tables (copied).
- BeaverTails ``is_safe``: real safety-auditor labels.
- WildGuard train ``response_harm_label``: real response-harm labels.
- ToxicChat ``toxicity`` (human-annotated): real moderation-queue labels.

Usage::

    python3 scripts/build_llm_audit_xy.py
"""
from __future__ import annotations

import csv
import json
import re
import urllib.request
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "hf_cache" / "audit"
OUT = ROOT / "results" / "manuscript" / "llm_audit"

X_DIMS = [
    "n_toks",
    "n_chars",
    "avg_word",
    "qmark",
    "bang",
    "hedge",
    "formal",
    "i_count",
    "newlines",
    "upper",
    "refuse",
    "please",
    "thank",
]
X_COLS = [f"x_{d}" for d in X_DIMS]

_HEDGE = ("maybe", "perhaps", "i think", "not sure", "probably", "might")
_FORMAL = ("therefore", "however", "furthermore", "regarding", "consequently")
_REFUSE = (
    "i can't",
    "i cannot",
    "i'm sorry",
    "i am sorry",
    "i won't",
    "i will not",
    "against my",
    "not able to",
)

PARQUET = {
    "beavertails": "https://huggingface.co/datasets/PKU-Alignment/BeaverTails/resolve/refs%2Fconvert%2Fparquet/default/30k_test/0000.parquet",
    "toxicchat": "https://huggingface.co/datasets/lmsys/toxic-chat/resolve/refs%2Fconvert%2Fparquet/toxicchat0124/train/0000.parquet",
    "wildguard": "https://huggingface.co/datasets/ToxicityPrompts/wildguard-train/resolve/refs%2Fconvert%2Fparquet/wildguardtrain/train/0000.parquet",
}

UA = {"User-Agent": "cfperm-llm-audit-xy"}


def style_vector(text: str) -> np.ndarray:
    t = (text or "").strip()
    low = t.lower()
    toks = re.findall(r"[a-zA-Z']+", low)
    chars = max(len(t), 1)
    avg_w = float(np.mean([len(w) for w in toks])) if toks else 0.0
    return np.asarray(
        [
            len(toks) / 200.0,
            chars / 800.0,
            avg_w / 10.0,
            low.count("?") / 5.0,
            low.count("!") / 5.0,
            sum(low.count(h) for h in _HEDGE) / 5.0,
            sum(low.count(h) for h in _FORMAL) / 5.0,
            low.count(" i ") / 10.0,
            t.count("\n") / 10.0,
            (sum(c.isupper() for c in t) / chars),
            sum(low.count(h) for h in _REFUSE) / 3.0,
            low.count("please") / 4.0,
            low.count("thank") / 4.0,
        ],
        dtype=float,
    )


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


def write_xy(path: Path, y: np.ndarray, X: np.ndarray, batch: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["y", "batch", *X_COLS])
        for i in range(len(y)):
            w.writerow([int(y[i]), int(batch[i]), *[f"{x:.6g}" for x in X[i]]])


def pack(y: np.ndarray, X: np.ndarray, n: int = 1200, n_per: int = 80):
    n_use = min(len(y), n)
    n_use = (n_use // n_per) * n_per
    y, X = y[:n_use], X[:n_use]
    batch = np.repeat(np.arange(n_use // n_per), n_per)
    return y, X, batch


def from_texts_and_y(texts, y_list, n=1200, n_per=80):
    y = np.asarray(y_list, dtype=int)
    X = np.vstack([style_vector(t) for t in texts])
    return pack(y, X, n=n, n_per=n_per)


def load_beavertails(n=1200):
    path = download(PARQUET["beavertails"], CACHE / "beavertails_30k_test.parquet")
    tbl = pq.read_table(path, columns=["response", "is_safe"])
    texts, ys = [], []
    for resp, safe in zip(tbl.column("response").to_pylist(), tbl.column("is_safe").to_pylist()):
        if resp is None:
            continue
        texts.append(str(resp))
        ys.append(1 if bool(safe) else 0)
        if len(ys) >= n:
            break
    return from_texts_and_y(texts, ys, n=n)


def load_wildguard(n=1200):
    path = download(PARQUET["wildguard"], CACHE / "wildguard_train.parquet")
    tbl = pq.read_table(
        path, columns=["response", "response_harm_label"]
    )
    texts, ys = [], []
    for resp, lab in zip(
        tbl.column("response").to_pylist(),
        tbl.column("response_harm_label").to_pylist(),
    ):
        if resp is None or lab is None:
            continue
        lab_s = str(lab).strip().lower()
        if lab_s not in {"harmful", "unharmful"}:
            continue
        texts.append(str(resp))
        ys.append(1 if lab_s == "unharmful" else 0)
        if len(ys) >= n:
            break
    return from_texts_and_y(texts, ys, n=n)


def load_toxicchat(n=1200):
    path = download(PARQUET["toxicchat"], CACHE / "toxicchat0124_train.parquet")
    tbl = pq.read_table(
        path, columns=["model_output", "toxicity", "human_annotation"]
    )
    texts, ys = [], []
    for out, tox, human in zip(
        tbl.column("model_output").to_pylist(),
        tbl.column("toxicity").to_pylist(),
        tbl.column("human_annotation").to_pylist(),
    ):
        if not bool(human):
            continue
        if out is None or tox is None:
            continue
        texts.append(str(out))
        ys.append(1 if int(tox) == 0 else 0)
        if len(ys) >= n:
            break
    return from_texts_and_y(texts, ys, n=n)


def summarize(name: str, y, batch) -> dict:
    return {
        "name": name,
        "n": int(len(y)),
        "n_batches": int(batch.max() + 1) if len(batch) else 0,
        "y_pass_rate": float(np.mean(y)) if len(y) else None,
        "y_fail_rate": float(1.0 - np.mean(y)) if len(y) else None,
        "schema": ["y", "batch", *X_COLS],
        "y_meaning": "1=pass (over / safe / unharmful / not toxic), 0=fail",
    }


def stack_two_xy(src0: Path, src1: Path, dest: Path, name0: str, name1: str) -> Path:
    """CFPerm-ready two-stream table. T=0/1 are reviewer queues, Y is the audit decision."""
    rows = []
    for t, src, stream in ((0, src0, name0), (1, src1, name1)):
        with src.open() as f:
            r = csv.DictReader(f)
            for row in r:
                rows.append(
                    {
                        "T": t,
                        "Y": int(row["y"]),
                        **{c: row[c] for c in X_COLS},
                        "batch": int(row["batch"]),
                        "stream": stream,
                    }
                )
    with dest.open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["T", "Y", *X_COLS, "batch", "stream"]
        )
        w.writeheader()
        w.writerows(rows)
    return dest


def stack_hh_streams() -> Path:
    """CFPerm-ready two-stream table. T=0 helpful, T=1 harmless. Y is auditor, not chosen."""
    return stack_two_xy(
        OUT / "xy_hh_helpful_consistent.csv",
        OUT / "xy_hh_harmless_consistent.csv",
        OUT / "xy_hh_two_stream.csv",
        "helpful",
        "harmless",
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    manifest = {
        "note": "Prediction tables for manuscript use-case 4 (LLM moderation). HH chosen is not Y.",
        "x_dims": X_DIMS,
        "files": {},
    }

    builders = [
        ("xy_beavertails.csv", "BeaverTails 30k_test is_safe", load_beavertails),
        ("xy_wildguard.csv", "WildGuard train response_harm_label", load_wildguard),
        ("xy_toxicchat.csv", "ToxicChat0124 human toxicity", load_toxicchat),
    ]
    for fname, title, fn in builders:
        print("build", title)
        y, X, batch = fn()
        path = OUT / fname
        write_xy(path, y, X, batch)
        info = summarize(title, y, batch)
        info["path"] = str(path.relative_to(ROOT))
        info["cite"] = {
            "xy_beavertails.csv": "ji2023beavertails",
            "xy_wildguard.csv": "han2024wildguard",
            "xy_toxicchat.csv": "lin2023toxicchat",
        }[fname]
        manifest["files"][fname] = info
        print(" ", info)

    for fname in (
        "xy_hh_helpful_consistent.csv",
        "xy_hh_helpful_hop.csv",
        "xy_hh_harmless_consistent.csv",
        "xy_hh_harmless_hop.csv",
    ):
        path = OUT / fname
        if not path.exists():
            continue
        with path.open() as f:
            rows = list(csv.DictReader(f))
        y = np.array([int(r["y"]) for r in rows])
        batch = np.array([int(r["batch"]) for r in rows])
        stream = "helpful" if "helpful" in fname else "harmless"
        regime = "hop" if "hop" in fname else "consistent"
        manifest["files"][fname] = {
            **summarize(f"HH-RLHF {stream} {regime} (prototype auditor Y)", y, batch),
            "path": str(path.relative_to(ROOT)),
            "cite": "bai2022hh-rlhf",
            "y_is_chosen": False,
            "regime": regime,
            "stream": stream,
        }

    dest = stack_hh_streams()
    with dest.open() as f:
        n = sum(1 for _ in f) - 1
    manifest["files"][dest.name] = {
        "name": "HH two-stream CFPerm table (T=helpful/harmless, Y=auditor)",
        "n": n,
        "path": str(dest.relative_to(ROOT)),
        "T": "0=helpful queue, 1=harmless queue",
        "y_is_chosen": False,
    }

    dest_real = stack_two_xy(
        OUT / "xy_beavertails.csv",
        OUT / "xy_toxicchat.csv",
        OUT / "xy_real_two_stream.csv",
        "beavertails",
        "toxicchat",
    )
    with dest_real.open() as f:
        n_real = sum(1 for _ in f) - 1
    manifest["files"][dest_real.name] = {
        "name": "Real-label two-stream CFPerm table (T=BeaverTails/ToxicChat)",
        "n": n_real,
        "path": str(dest_real.relative_to(ROOT)),
        "T": "0=BeaverTails is_safe, 1=ToxicChat human toxicity",
        "y_is_chosen": False,
    }

    man_path = OUT / "MANIFEST.json"
    man_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print("wrote", man_path)


if __name__ == "__main__":
    main()
