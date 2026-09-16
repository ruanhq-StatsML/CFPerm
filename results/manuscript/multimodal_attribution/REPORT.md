# Multimodal attribution — permute T, report blocks

Same CFPerm loop as group attribution. The unit is a **modality block**, not a single x_*.
X = concat(text, vision/dense, graph, fusion). T is the pack/group. Permute T.
No pixels, no raw prompts, no HH chosen.

```bash
PYTHONPATH=. python3 scripts/prototype_multimodal_attribution.py
```

| Stream | T | reject | block hits | top blocks |
|---|---|---|---|---|
| hybrid: dense pack hop | T=0 before dense cut, T=1 after | no | — | fusion, query, sparse |
| Graph-RAG: local vs community | T=0 local pack, T=1 community pack | yes | graph | graph, query |

## Synthetic check

| DGP | reject | top blocks |
|---|---|---|
| planted CATE on vision block | yes | vision, text, graph |
| null: Y depends on X, not T | no | text, graph, vision |

Planted: only the vision block interacts with T. Null: Y depends on text/graph, not on T.

Not this object: time-gate hop, single-channel Recall, HH chosen, raw image/audio.

