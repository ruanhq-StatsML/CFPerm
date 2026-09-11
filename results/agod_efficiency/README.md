# AGOD efficiency — Hugging Face datasets

| Dataset | Modalities | Shift axis |
|---|---|---|
| `PolyAI/minds14` | audio + text | locale: en-US → fr/de/es/it/nl/pt |
| `ChristophSchuhmann/MS_COCO_2017_URL_TEXT` | image + text | semantic cohort: person → vehicle/food |

**Metric:** hard-gate modality heads with `α_m < θ`; relative FLOPs = |active|/|M|; efficiency = drift-coverage / rel-FLOPs.

See `AGOD_HF_Efficiency_Dashboard.png` and `AGOD_hf_efficiency_tables_only.tex`.
