#!/usr/bin/env python3
"""Export 300-example infer-bench slices and write characterization JSON."""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPORT = ROOT / "data" / "hf_cache" / "infer_bench_export"
RESULTS = ROOT / "results" / "infer_bench"
N = 300
PAGE = 100
SEED_NOTE = "deterministic head/mid/tail windows of 100 over the source split"

DATASETS = [
    {
        "name": "HaluEval",
        "stem": "halueval",
        "hub": "pminervini/HaluEval",
        "config": "qa_samples",
        "split": "data",
        "cite": "li-etal-2023-halueval",
        "task": "hallucination recognition (QA)",
        "license": "see RUCAIBox/HaluEval",
    },
    {
        "name": "SQuAD",
        "stem": "squad",
        "hub": "rajpurkar/squad",
        "config": "plain_text",
        "split": "validation",
        "cite": "rajpurkar-etal-2016-squad",
        "task": "extractive reading comprehension",
        "license": "CC BY-SA 4.0",
    },
    {
        "name": "HotpotQA",
        "stem": "hotpotqa",
        "hub": "hotpotqa/hotpot_qa",
        "config": "distractor",
        "split": "validation",
        "cite": "yang-etal-2018-hotpotqa",
        "task": "multi-hop QA + supporting facts",
        "license": "CC BY-SA 4.0",
    },
    {
        "name": "TruthfulQA",
        "stem": "truthfulqa",
        "hub": "truthfulqa/truthful_qa",
        "config": "generation",
        "split": "validation",
        "cite": "lin-etal-2022-truthfulqa",
        "task": "truthfulness / imitative falsehoods",
        "license": "Apache-2.0",
    },
]


def api(path: str, **params):
    url = "https://datasets-server.huggingface.co/" + path + "?" + urllib.parse.urlencode(params)
    delay = 2.0
    last = None
    for _ in range(8):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            last = e
            if e.code not in (429, 500, 502, 503):
                raise
            time.sleep(delay)
            delay = min(delay * 2, 32)
    raise last


def fetch_indices(dataset, config, split, indices):
    wanted = set(indices)
    out = {}
    runs = []
    seq = sorted(indices)
    start = prev = seq[0]
    for i in seq[1:]:
        if i == prev + 1:
            prev = i
            continue
        runs.append((start, prev))
        start = prev = i
    runs.append((start, prev))
    for lo, hi in runs:
        offset = lo
        while offset <= hi:
            chunk = api(
                "rows",
                dataset=dataset,
                config=config,
                split=split,
                offset=offset,
                length=min(PAGE, hi - offset + 1),
            )
            for item in chunk["rows"]:
                idx = item["row_idx"]
                if idx in wanted:
                    out[idx] = item["row"]
            if not chunk["rows"]:
                break
            offset += PAGE
            time.sleep(0.5)
    missing = [i for i in indices if i not in out]
    if missing:
        raise RuntimeError(f"{dataset}: missing rows {missing[:8]}...")
    return [out[i] for i in indices]


def n_words(text: str) -> int:
    return len(str(text).split())


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def compact_row(stem: str, row: dict) -> dict:
    if stem == "halueval":
        return {
            "knowledge": row.get("knowledge"),
            "question": row.get("question"),
            "answer": row.get("answer"),
            "hallucination": row.get("hallucination"),
        }
    if stem == "squad":
        answers = row.get("answers") or {}
        texts = answers.get("text") or []
        return {
            "id": row.get("id"),
            "title": row.get("title"),
            "context": row.get("context"),
            "question": row.get("question"),
            "answer": texts[0] if texts else "",
        }
    if stem == "hotpotqa":
        ctx = row.get("context") or {}
        n_paras = len(ctx.get("title") or [])
        return {
            "id": row.get("id"),
            "question": row.get("question"),
            "answer": row.get("answer"),
            "type": row.get("type"),
            "level": row.get("level"),
            "n_context_titles": n_paras,
            "supporting_facts": row.get("supporting_facts"),
            "context": ctx,
        }
    if stem == "truthfulqa":
        return {
            "type": row.get("type"),
            "category": row.get("category"),
            "question": row.get("question"),
            "best_answer": row.get("best_answer"),
            "n_correct": len(row.get("correct_answers") or []),
            "n_incorrect": len(row.get("incorrect_answers") or []),
            "source": row.get("source"),
        }
    return row


def characterize(stem: str, rows: list[dict]) -> dict:
    q_lens = [n_words(r.get("question") or "") for r in rows]
    stats = {
        "n": len(rows),
        "question_words_mean": round(mean(q_lens), 2),
        "question_words_p50": sorted(q_lens)[len(q_lens) // 2],
    }
    if stem == "halueval":
        labels = [(r.get("hallucination") or "").lower() for r in rows]
        yes = sum(x == "yes" for x in labels)
        stats["hallucination_yes"] = yes
        stats["hallucination_yes_rate"] = round(yes / len(rows), 4)
        stats["answer_words_mean"] = round(mean(n_words(r.get("answer") or "") for r in rows), 2)
        stats["knowledge_words_mean"] = round(
            mean(n_words(r.get("knowledge") or "") for r in rows), 2
        )
    elif stem == "squad":
        ans = []
        ctx = []
        titles = set()
        for r in rows:
            titles.add(r.get("title"))
            ctx.append(n_words(r.get("context") or ""))
            answers = (r.get("answers") or {}).get("text") or [""]
            ans.append(n_words(answers[0] if answers else ""))
        stats["n_titles"] = len(titles)
        stats["context_words_mean"] = round(mean(ctx), 2)
        stats["answer_words_mean"] = round(mean(ans), 2)
    elif stem == "hotpotqa":
        types, levels, n_ctx = {}, {}, []
        for r in rows:
            types[r.get("type") or "unk"] = types.get(r.get("type") or "unk", 0) + 1
            levels[r.get("level") or "unk"] = levels.get(r.get("level") or "unk", 0) + 1
            ctx = r.get("context") or {}
            n_ctx.append(len(ctx.get("title") or []))
        stats["type_counts"] = types
        stats["level_counts"] = levels
        stats["n_context_titles_mean"] = round(mean(n_ctx), 2)
        stats["answer_words_mean"] = round(mean(n_words(r.get("answer") or "") for r in rows), 2)
    elif stem == "truthfulqa":
        cats, types = {}, {}
        for r in rows:
            cats[r.get("category") or "unk"] = cats.get(r.get("category") or "unk", 0) + 1
            types[r.get("type") or "unk"] = types.get(r.get("type") or "unk", 0) + 1
        stats["type_counts"] = types
        stats["n_categories"] = len(cats)
        stats["top_categories"] = dict(sorted(cats.items(), key=lambda kv: -kv[1])[:8])
        stats["best_answer_words_mean"] = round(
            mean(n_words(r.get("best_answer") or "") for r in rows), 2
        )
    return stats


def even_indices(n_total: int, n: int, windows: int = 3) -> list[int]:
    """Head / mid / tail windows of 100, to keep Hub calls small and the slice mixed."""
    if n_total <= n:
        return list(range(n_total))
    per = n // windows
    starts = [int(i * (n_total - per) / max(windows - 1, 1)) for i in range(windows)]
    idx = []
    for s in starts:
        idx.extend(range(s, min(n_total, s + per)))
    while len(idx) < n and idx[-1] + 1 < n_total:
        idx.append(idx[-1] + 1)
    return idx[:n]


def main():
    EXPORT.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    summary = {
        "n_per_dataset": N,
        "sampling": SEED_NOTE,
        "datasets": [],
        "online_rfperm": {
            "source": "results/agod/hf_landing on cursor/hf-llm-landing-protos-92f8",
            "note": "Existing OnlineRFPerm prototype numbers. HaluEval uses qa_samples, not the 300-export slice.",
            "halueval": {
                "dataset": "pminervini/HaluEval@qa_samples",
                "n": 2400,
                "n_batches": 24,
                "n_per": 100,
                "cut_batch": 4,
                "gate": 1.25,
                "halluc_rate_before": 0.080,
                "halluc_rate_after": 0.942,
                "fire_rate": 0.0909,
                "first_fire_t": 4,
                "cut_ratio": 11.714,
                "auroc_po_risk0_at_cut": 1.0,
                "precision_at_10_at_cut": 1.0,
                "instance_probe_auroc_post": 1.0,
                "domain_auc_pre_post": 0.811,
            },
            "hh_rlhf": {
                "dataset": "Anthropic/hh-rlhf@helpful-base",
                "n": 4000,
                "cite": "bai2022traininghelpfulharmlessassistant",
                "judge_err_pre": 0.147,
                "judge_err_post": 0.539,
                "judge_err_ratio": 3.671,
                "style_domain_auc": 0.999,
                "text_domain_auc": 0.494,
                "fire_rate": 0.104,
                "first_fire_t": 4,
            },
        },
    }
    for spec in DATASETS:
        info = api("size", dataset=spec["hub"], config=spec["config"])
        # size payload varies; fall back to rows
        n_total = None
        try:
            n_total = info["size"]["config"]["num_rows"]
        except Exception:
            pass
        preview = api(
            "rows",
            dataset=spec["hub"],
            config=spec["config"],
            split=spec["split"],
            offset=0,
            length=1,
        )
        n_split = preview.get("num_rows_total") or n_total
        idx = even_indices(int(n_split), N)
        raw = fetch_indices(spec["hub"], spec["config"], spec["split"], idx)
        compact = [compact_row(spec["stem"], r) for r in raw]
        outp = EXPORT / f"{spec['stem']}.jsonl"
        with outp.open("w", encoding="utf-8") as f:
            for row in compact:
                f.write(json.dumps(row, ensure_ascii=False) + "\n")
        stats = characterize(spec["stem"], raw)
        rec = {
            **spec,
            "n_source_split": int(n_split),
            "n_export": len(compact),
            "jsonl": str(outp.relative_to(ROOT)),
            "bytes": outp.stat().st_size,
            "stats": stats,
        }
        summary["datasets"].append(rec)
        print(spec["name"], "export", len(compact), "from", n_split, "->", outp)
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("wrote", RESULTS / "summary.json")


if __name__ == "__main__":
    main()
