# Infer-bench dataset citations and 300-example exports

Core four (this folder):

| Dataset | JSONL | Cite key | Venue |
|---|---|---|---|
| HaluEval | `halueval.jsonl` | `li-etal-2023-halueval` | EMNLP 2023 |
| SQuAD 1.1 | `squad.jsonl` | `rajpurkar-etal-2016-squad` | EMNLP 2016 |
| HotpotQA | `hotpotqa.jsonl` | `yang-etal-2018-hotpotqa` | EMNLP 2018 |
| TruthfulQA | `truthfulqa.jsonl` | `lin-etal-2022-truthfulqa` | ACL 2022 |

Each JSONL is a **deterministic $n=300$** head/mid/tail window over the Hugging Face split (see `results/infer_bench/summary.json`). Rebuild with:

```bash
python3 scripts/infer_bench/export_and_summarize.py
python3 scripts/infer_bench/write_tex.py
```

## BibTeX

- `recommended.bib` — the four venue keys (+ HH-RLHF if you cite the alignment stream)
- `infer_bench.bib` — venue + arXiv + Hub + GitHub + SQuAD 2.0 + HaluEval 2.0
- Per-dataset: `halueval.bib`, `squad.bib`, `hotpotqa.bib`, `truthfulqa.bib`
- Related only: `hh_rlhf.bib` (`bai2022hh-rlhf`)

HH-RLHF is **not** a fifth 300-export. It is the alignment counterpart already run in the HF landing prototype.

## LLM audit-logic hop

审核逻辑突变（judge / 政策 / 人审池一夜换制度）写在：

- `docs/reports/LLM_Judge_Audit_Hop.md`
- `docs/reports/LLM_Judge_Audit_Hop.tex`

Last-two 问的是审核链路还顺不顺；相对上线参考问的是这套审核还信不信得过。HH-RLHF 上 quiet \(t=2\) 不响，cut \(t=4\) 当窗 fire。

## LaTeX

Paste-ready tables: `docs/reports/Infer_Bench_tables_only.tex`  
Wrapped note: `docs/reports/Infer_Bench_Datasets.tex`

```latex
\bibliography{data/hf_cache/infer_bench_export/recommended}
\cite{li-etal-2023-halueval,rajpurkar-etal-2016-squad,yang-etal-2018-hotpotqa,lin-etal-2022-truthfulqa}
```

## Other public sets in the same pipeline

Only add these if you actually use them:

- **HH-RLHF** (Anthropic helpful/harmless prefs) — already has OnlineRFPerm judge/style numbers
- **HaluEval 2.0** — `li-etal-2024-dawn`
- **SQuAD 2.0** — `rajpurkar-etal-2018-know` (unanswerable questions)
- HaluEval's other tasks (dialogue / summarization) sit on OpenDialKG and CNN/DailyMail; not exported here
