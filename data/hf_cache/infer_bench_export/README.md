# Infer-bench dataset citations

BibTeX for the four 300-example exports under `data/hf_cache/infer_bench_export/`.

| Dataset | Local export | Recommended cite key | Venue |
|---|---|---|---|
| HaluEval | `halueval.jsonl` | `li-etal-2023-halueval` | EMNLP 2023 |
| SQuAD | `squad.jsonl` | `rajpurkar-etal-2016-squad` | EMNLP 2016 |
| HotpotQA | `hotpotqa.jsonl` | `yang-etal-2018-hotpotqa` | EMNLP 2018 |
| TruthfulQA | `truthfulqa.jsonl` | `lin-etal-2022-truthfulqa` | ACL 2022 |

## Files

- `infer_bench.bib` — combined, self-contained bibliography
- `halueval.bib` — HaluEval + HaluEval 2.0 + GitHub
- `squad.bib` — SQuAD 1.1 + SQuAD 2.0 + Hub / explorer
- `hotpotqa.bib` — HotpotQA venue, arXiv, Hub, homepage
- `truthfulqa.bib` — TruthfulQA venue, arXiv, Hub, GitHub

Venue entries include ACL Anthology metadata (DOI, pages, ISBN, editors), arXiv eprint, abstracts, and resource URLs. Dataset-card keys (`yang2018hotpotqa`, `lin2021truthfulqa`) are kept as aliases.

If the SQuAD export is v2 (unanswerable questions), also cite `rajpurkar-etal-2018-know`.

## LaTeX

```latex
\bibliographystyle{acl_natbib}
\bibliography{data/hf_cache/infer_bench_export/infer_bench}
```

```latex
\cite{li-etal-2023-halueval,rajpurkar-etal-2016-squad,yang-etal-2018-hotpotqa,lin-etal-2022-truthfulqa}
```
