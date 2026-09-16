# HF landing prototypes (today)

Two HuggingFace subsets + local `agod` OnlineRFPerm/PO-risk:

| Demo | Data | What fires |
|------|------|------------|
| RFPerm-as-Judge + 文风 | `Anthropic/hh-rlhf` helpful-base (2.5k) | preference-map hop + register portrait AUC |
| Hallucination regime + RAG/router | `pminervini/HaluEval` qa_samples (3k) | label-law hop at cut + `po_risk0` ranking |

```bash
# cache is already under data/hf_cache/; re-pull if missing
PYTHONPATH=. python3 scripts/agod/hf_landing_protos.py --gate 1.25
pytest -q tests/test_hf_landing_protos.py
```

Outputs: `results/agod/hf_landing/{REPORT.md,summary.json,hh_judge_style.json,halu_regime_rag.json}`

**Stance:** regime detection / ranking / gated reweight — not fact-checking, not unique causal attribution.
