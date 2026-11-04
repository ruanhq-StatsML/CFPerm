# Infer-bench results

Characterization of the four $n=300$ exports plus OnlineRFPerm numbers already shipped on `cursor/hf-llm-landing-protos-92f8`.

| File | Role |
|---|---|
| `summary.json` | Export stats + OnlineRFPerm payload |
| `Infer_Bench_tables_only.tex` | Paste-ready tables (copy of `docs/reports/Infer_Bench_tables_only.tex`) |

Rebuild:

```bash
python3 scripts/infer_bench/export_and_summarize.py
python3 scripts/infer_bench/write_tex.py
```

OnlineRFPerm:

- HaluEval QA stream ($n=2400$, not the 300-slice): halluc. rate $0.080\to 0.942$, fire at $t=4$, $\texttt{po\_risk0}$ AUROC $1.0$ at the cut
- HH-RLHF helpful-base ($n=4000$): judge err $0.147\to 0.539$, style-domain AUC $0.999$

SQuAD / HotpotQA / TruthfulQA 300-slices are posted as **dataset characterization** here; OnlineRFPerm has not been re-run on those three slices.
