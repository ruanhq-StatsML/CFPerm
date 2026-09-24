# TencentGR local subset (`data/tencent_subset`)

Working slice used by CFPerm TencentGR localization / continuous-time AD
prototypes. **Not** the full Hugging Face dump — only `seq` / `item_feat` /
`user_feat` (no `candidate`, no `mm_emb`).

## Download / browse (this repo)

After this commit lands on the branch:

```text
https://github.com/ruanhq-StatsML/CFPerm/tree/cursor/graph-localize-elab-abce/data/tencent_subset
```

Parquet blobs are stored with **Git LFS**. Clone with LFS:

```bash
git lfs install
git clone https://github.com/ruanhq-StatsML/CFPerm.git
cd CFPerm
git checkout cursor/graph-localize-elab-abce   # or main after merge
git lfs pull
```

Or fetch only this folder via sparse + LFS once the branch is available.

## Layout

```text
data/tencent_subset/
  seq/*.parquet          # user behavior sequences  (~261 MB, 3 parts)
  item_feat/*.parquet    # item side features        (~11 MB)
  user_feat/*.parquet    # user side features        (~1 MB)
  README.md              # this file
```

## Snapshot stats (this checkout)

| | |
|---|---|
| users in `seq` | **296 595** |
| events (seq length sum) | **≈ 26.7M** |
| timeline span | **≈ 231.3 days** |
| `item_feat` rows | **478 613** |
| `user_feat` rows | **100 315** |
| on-disk | **≈ 272 MB** |

Source full dataset: [TAAC2025/TencentGR-1M](https://huggingface.co/datasets/TAAC2025/TencentGR-1M)
(CC-BY-4.0). Labels: exposure(`action_type=0`) / click(`1`). Terminal success
for our localize scripts = **click**.

## Schema (same as upstream)

### `seq`
| field | type | notes |
|---|---|---|
| `user_id` | int64 | RID |
| `seq` | List\[Dict\] | each: `item_id`, `action_type`, `timestamp` |

Expanded edge: `(user_id, item_id, action_type, timestamp)`.

### `item_feat`
`item_id` + encrypted cols `100`…`122`. **Merchant proxy = column `122`.**

### `user_feat`
`user_id` + cols `103`…`110` (some List).

## Code interface

```python
from pathlib import Path
from scripts.tencent_gr.time_window_feats import (
    scan_time_range, feature_engineer, fsds_feature_columns,
)
from scripts.tencent_gr.run_three_step_subset_localize import (
    load_item_merchant_map, attach_merchant,
)

ROOT = Path("data/tencent_subset")
t_min, t_max, n = scan_time_range(ROOT, max_users=20_000)
panel = feature_engineer(ROOT, t_min, t_min + 45 * 86400, max_users=20_000)
grid = attach_merchant(panel["grid"], load_item_merchant_map(ROOT))
```

Drill recipe: **rank-average(MMD, PO, cmean)** → merchant → user → order → one FSDS.  
See `docs/summaries/Graph_CT_AD_Interface.md` and
`docs/summaries/Graph_Localize_Algorithm_Elaboration.md`.

## CLI smoke

```bash
PYTHONPATH=. python3 scripts/tencent_gr/run_three_step_subset_localize.py \
  --root data/tencent_subset --max-users 20000 --gap-days 30
```

## Provenance

Subset prepared for CFPerm agent runs (windowed FE + localize prototypes).
Does **not** claim to be an official TAAC release cut; for citation / license
follow the upstream HF card.
