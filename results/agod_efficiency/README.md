# AGOD efficiency — Hugging Face datasets

## What this measures

**Online adaptation / distillation training cost**, not student **inference latency**.

| Phase | AGOD changes it? | Meaning |
|---|---|---|
| Serving / inference forward | **No** (default) | Deployed student still runs full modalities unless you add a separate sparse student. |
| Online distill step \(t\) | **Yes** | MSG → Softmax `α` → hard-gate: skip `L_m` (+ backward) when `α_m < θ`. |
| Stable-modality forgetting | **Indirectly** | Gated heads are not updated → less overwrite on non-drifting modalities. |

Metric: relative **adaptation** FLOPs = `|active update heads| / |M|`; efficiency = drift-coverage / rel-FLOPs.

## Datasets

| Dataset | Modalities | Shift axis |
|---|---|---|
| `PolyAI/minds14` | audio + text | locale: en-US → fr/de/es/it/nl/pt |
| `ChristophSchuhmann/MS_COCO_2017_URL_TEXT` | image + text | semantic cohort: person → vehicle/food |

See `AGOD_HF_Efficiency_Dashboard.png` and `AGOD_hf_efficiency_tables_only.tex`.
