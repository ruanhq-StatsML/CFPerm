# Mag×Var multi-dataset justify (PermuCATE PO)

Concrete multi-dataset justify: same Mag×Var + Y-flip口径 on PermuCATE u across AIGC / topic / news / sentiment boards.

| board | family | sev(ben) | sev(mal) | gap | OK? |
| --- | --- | --- | --- | --- | --- |
| hc3_finance | AIGC integrity | 0.015 | 0.791 | 0.776 | yes |
| hc3_medicine | AIGC integrity | 0.498 | 0.800 | 0.302 | yes |
| hc3_open_qa | AIGC integrity | 0.255 | 1.044 | 0.789 | yes |
| hc3_wiki | AIGC integrity | 0.004 | 0.430 | 0.425 | yes |
| news20_sci | topic monitor | 0.252 | 0.312 | 0.060 | yes |
| news20_talk | topic monitor | -0.052 | 0.396 | 0.448 | yes |
| news20_comp | topic monitor | 0.015 | 0.816 | 0.801 | yes |
| ag_news | news topic | 0.130 | 0.511 | 0.381 | yes |
| imdb | sentiment | -0.032 | 0.418 | 0.450 | yes |

hit-rate: **1.00** (9/9)

## By family

| family | hit-rate |
| --- | --- |
| AIGC integrity | 1.00 (4/4) |
| topic monitor | 1.00 (3/3) |
| news topic | 1.00 (1/1) |
| sentiment | 1.00 (1/1) |

Across HC3 verticals, 20 Newsgroups pairs, AG News, and IMDB, malignant severity exceeds benign under the same Y-flip口径 — concrete evidence the product rule is not board-specific.

```bash
python3 scripts/run_as_dta_magvar_dataset_justify.py --quick
```
