# Gap-guided LLM inference (budgeted-block stand-in)

Frozen π / s_m(x) route a one-block expert — the same decision as enabling one tool or packing one modality into the prompt.
No vendor LLM is called. Wide text + concentrated valence shift is the regime where raw VIMP overweights text.

## Held-out domain AUC (one-block budget, mean over seeds)

| setting | oracle | blend | π | VIMP | random |
|---|---:|---:|---:|---:|---:|
| synthetic GT=valence | 0.959 | 0.957 | 0.959 | 0.959 | 0.642 |
| inject valence | 0.747 | 0.745 | 0.747 | 0.600 | 0.583 |
| shuffle valence (neg.) | 0.505 | 0.481 | 0.474 | 0.470 | 0.516 |

## P(selected block = GT)

| setting | blend | π | instance | VIMP | random |
|---|---:|---:|---:|---:|---:|
| synthetic GT=valence | 0.980 | 1.000 | 0.749 | 1.000 | 0.231 |
| inject valence | 0.824 | 1.000 | 0.465 | 0.000 | 0.245 |
| shuffle valence (neg.) | 0.128 | 0.333 | 0.191 | 0.000 | 0.256 |

Copy-paste prompt: `example_system_prompt.txt`. Packet JSON: `example_packet.json`.

```bash
python3 scripts/run_gap_guided_inference.py
```
