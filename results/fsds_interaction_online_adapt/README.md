# FSDS × Interaction Value Online Adaptation

Localization names the sick modality; **clever covariates** treat at
inference without retrain:

```
score_adapted = score_raw / (1 + λ · W_shift)
```

PO-risk uses the FSDS crossfit with **RF outcome** m̂ and RF propensity ê:
`po=(Y−m̂)·(W−ê)`. W_shift mixes |po|, triangle inconsistency, and Φ drop.
`fsds_clever_cov` gates λ with typed hop `(ĉ, δ̂)` (locked TSS signs).

## Live numbers

| scenario | static | localize_only | clever_cov | **fsds_clever_cov** | lift vs localize |
| --- | ---: | ---: | ---: | ---: | ---: |
| healthy | 0.3036 | 0.3036 | 0.3039 | **0.3040** | -0.0003 |
| interaction_break | 0.3121 | 0.3121 | 0.3173 | **0.3178** | -0.0056 |
| alignment_attack | 0.3141 | 0.3141 | 0.3215 | **0.3255** | -0.0114 |
| concept_drift | 0.4458 | 0.4458 | 0.4040 | **0.3976** | 0.0483 |
| compound | 0.3122 | 0.3122 | 0.3137 | **0.3147** | -0.0025 |

## What each method does

- `static` — RF joint scores, no gate.
- `localize_only` — π / Φ / triangle / RF-PO (diagnosis) but leaves preds alone.
- `clever_cov` — KS on RF |po| & alignment → inject W into scores.
- `fsds_clever_cov` — same W, λ from `(ĉ,δ̂,Δπ)` + triangle/Φ collapse.
- `oracle_penalty` — upper bound using true |error| as W.

## Hop interpretations (`fsds_clever_cov`)

After each hop: RF-PO modality π, top Φ tokens, triangle, typed (c,δ), read.

### `healthy`

| hop | r_po | π_po dom | tri | top Φ_TI token | c | δ | read |
| ---: | ---: | --- | ---: | --- | ---: | ---: | --- |
| 1 | 0.0921 | text | 0.405 | k=9 (Φ=-0.0141) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 2 | 0.0746 | text | 0.428 | k=2 (Φ=-0.0046) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 3 | 0.0910 | text | 0.405 | k=2 (Φ=0.0060) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 4 | 0.0955 | text | 0.439 | k=2 (Φ=0.0001) | 1.000 | 0.003 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 5 | 0.0789 | text | 0.411 | k=1 (Φ=-0.0017) | 1.000 | 0.006 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 6 | 0.0805 | text | 0.417 | k=1 (Φ=0.0016) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 7 | 0.0848 | text | 0.416 | k=13 (Φ=0.0018) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |

### `interaction_break`

| hop | r_po | π_po dom | tri | top Φ_TI token | c | δ | read |
| ---: | ---: | --- | ---: | --- | ---: | ---: | --- |
| 1 | 0.0798 | text | 0.405 | k=9 (Φ=-0.0121) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 2 | 0.0733 | text | 0.428 | k=2 (Φ=-0.0035) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 3 | 0.0851 | text | 0.405 | k=2 (Φ=0.0072) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 4 | 0.0958 | text | 0.442 | k=2 (Φ=-0.0001) | 1.000 | 0.182 | concept/PO spike — RF PO-risk up; raise λ and review |
| 5 | 0.0994 | text | 0.383 | k=1 (Φ=0.0014) | 0.953 | 0.113 | concept/PO spike — RF PO-risk up; raise λ and review |
| 6 | 0.1072 | text | 0.385 | k=8 (Φ=0.0004) | 1.000 | 0.109 | concept/PO spike — RF PO-risk up; raise λ and review |
| 7 | 0.1029 | text | 0.392 | k=1 (Φ=0.0018) | 1.000 | 0.065 | concept/PO spike — RF PO-risk up; raise λ and review |

### `alignment_attack`

| hop | r_po | π_po dom | tri | top Φ_TI token | c | δ | read |
| ---: | ---: | --- | ---: | --- | ---: | ---: | --- |
| 1 | 0.0810 | text | 0.405 | k=2 (Φ=0.0219) | 0.982 | 0.001 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 2 | 0.0695 | text | 0.428 | k=2 (Φ=-0.0046) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 3 | 0.0756 | text | 0.405 | k=2 (Φ=0.0068) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 4 | 0.0596 | text | 0.456 | k=2 (Φ=0.0000) | 1.000 | 0.079 | alignment shift — triangle↑ / PO mass on text |
| 5 | 0.0819 | text | 0.467 | k=2 (Φ=0.0019) | 1.000 | 0.054 | alignment shift — triangle↑ / PO mass on text |
| 6 | 0.0883 | text | 0.447 | k=15 (Φ=-0.0035) | 1.000 | 0.023 | alignment shift — triangle↑ / PO mass on text |
| 7 | 0.0939 | text | 0.435 | k=11 (Φ=-0.0019) | 1.000 | 0.015 | concept/PO spike — RF PO-risk up; raise λ and review |

### `concept_drift`

| hop | r_po | π_po dom | tri | top Φ_TI token | c | δ | read |
| ---: | ---: | --- | ---: | --- | ---: | ---: | --- |
| 1 | 0.0935 | text | 0.405 | k=9 (Φ=-0.0138) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 2 | 0.0739 | text | 0.428 | k=2 (Φ=-0.0048) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 3 | 0.0925 | text | 0.405 | k=2 (Φ=0.0060) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 4 | 0.1238 | text | 0.439 | k=2 (Φ=0.0000) | 1.000 | 1.879 | concept/PO spike — RF PO-risk up; raise λ and review |
| 5 | 0.0913 | text | 0.405 | k=6 (Φ=-0.0006) | 0.981 | 0.606 | concept/PO spike — RF PO-risk up; raise λ and review |
| 6 | 0.0972 | text | 0.404 | k=6 (Φ=-0.0083) | 1.000 | 0.505 | concept/PO spike — RF PO-risk up; raise λ and review |
| 7 | 0.1042 | text | 0.437 | k=6 (Φ=0.0042) | 1.000 | 0.194 | concept/PO spike — RF PO-risk up; raise λ and review |

### `compound`

| hop | r_po | π_po dom | tri | top Φ_TI token | c | δ | read |
| ---: | ---: | --- | ---: | --- | ---: | ---: | --- |
| 1 | 0.0740 | text | 0.405 | k=2 (Φ=0.0202) | 0.994 | 0.001 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 2 | 0.0654 | text | 0.428 | k=2 (Φ=-0.0037) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 3 | 0.0741 | text | 0.405 | k=2 (Φ=0.0082) | 1.000 | 0.000 | covariate-leaning hop — shrink chase; bank/π not full relearn |
| 4 | 0.0554 | text | 0.447 | k=2 (Φ=0.0000) | 1.000 | 0.169 | concept/PO spike — RF PO-risk up; raise λ and review |
| 5 | 0.1080 | text | 0.456 | k=6 (Φ=-0.0003) | 1.000 | 0.139 | concept/PO spike — RF PO-risk up; raise λ and review |
| 6 | 0.1125 | text | 0.436 | k=13 (Φ=0.0046) | 1.000 | 0.092 | concept/PO spike — RF PO-risk up; raise λ and review |
| 7 | 0.1050 | text | 0.449 | k=1 (Φ=0.0024) | 1.000 | 0.064 | concept/PO spike — RF PO-risk up; raise λ and review |

## MSR-VTT hop probe

Source: `msrvtt` · n0=160 n1=160 · RF AUC=0.746 · r_po=0.0883

| π_text | π_video | π_audio | π_po | mean Φ_TI | triangle | mse_joint |
| ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 0.107 | 0.634 | 0.259 | text | 0.0119 | 0.457 | 2.8907 |

Read: mixed / mild hop — monitor Φ and PO (dom_po=text)

## Run

```
python3 scripts/run_fsds_interaction_online_adapt.py
python3 scripts/run_fsds_interaction_online_adapt.py --quick
```

Code: `Python/src/fsds_interaction_online_adapt.py`

Product line: RF Pseudo-Outcome Risk / Interaction Shift as **clever covariates**
closes the multimodal loop — diagnose with FSDS, treat in milliseconds.
