# Multimodal attribution — permute T, report blocks

Same CFPerm loop as group attribution. The unit is a **modality block**, not a single x_*.
X = concat(text, vision/dense, graph, fusion). T is the pack/group. Permute T.
No pixels, no raw prompts, no HH chosen.

```bash
PYTHONPATH=. python3 scripts/prototype_multimodal_attribution.py
```

| Stream | T | reject | block hits | top blocks | loc sig pairs |
|---|---|---|---|---|---:|
| hybrid: dense pack hop | T=0 before dense cut, T=1 after | no | — | fusion, query, sparse | 1 |
| Graph-RAG: local vs community | T=0 local pack, T=1 community pack | yes | graph | graph, query | 1 |

Post-hoc: subset indices, look at the **mean** (already computed), then pairwise **MMD** and **PO-risk**.

## Conditional mean (already computed — look at this)

| Stream | subset | n | mean Y | mean top x |
|---|---|---:|---:|---:|
| hybrid: dense pack hop | T0 | 320 | 0.572 | 0.11 |
| hybrid: dense pack hop | T1 | 880 | 0.450 | 0.107 |
| hybrid: dense pack hop | Q0 | 370 | 0.405 | 0.103 |
| hybrid: dense pack hop | Q1 | 277 | 0.401 | 0.104 |
| hybrid: dense pack hop | Q2 | 350 | 0.554 | 0.106 |
| hybrid: dense pack hop | Q3 | 203 | 0.611 | 0.125 |
| Graph-RAG: local vs community | T0 | 1200 | 0.534 | 0.43 |
| Graph-RAG: local vs community | T1 | 1200 | 0.425 | 0.43 |
| Graph-RAG: local vs community | Q0 | 716 | 0.235 | 0.0464 |
| Graph-RAG: local vs community | Q1 | 520 | 0.352 | 0.191 |
| Graph-RAG: local vs community | Q2 | 602 | 0.588 | 0.433 |
| Graph-RAG: local vs community | Q3 | 562 | 0.794 | 1.14 |

Read the mean first. Pairwise MMD and PO-risk are the significance next to it.

## Pairwise subset MMD / PO-risk

| Stream | pair | n | MMD | MMD p | PO-risk | PO p | mean Y |
|---|---|---|---:|---:|---:|---:|---|
| hybrid: dense pack hop | T0 vs T1 | 320/880 | 0.0505 | 0.0385 | 0.000206 | 0.962 | 0.572/0.450 |
| Graph-RAG: local vs community | T0 vs T1 | 1200/1200 | 0.0117 | 0.0385 | 0.00281 | 0.0385 | 0.534/0.425 |

## Synthetic check

| DGP | reject | top blocks |
|---|---|---|
| planted CATE on vision block | yes | vision, text, graph |
| null: Y depends on X, not T | no | text, graph, vision |

Planted: only the vision block interacts with T. Null: Y depends on text/graph, not on T.

Not this object: time-gate hop, single-channel Recall, HH chosen, raw image/audio.

