#!/usr/bin/env bash
# 30-min TencentGR-10M subset. Skip mm_emb (tens of GB; not used by FS).
set -euo pipefail
ROOT="${1:-data/tencent_subset}"
export PATH="${HOME}/.local/bin:${PATH}"
hf download TAAC2025/TencentGR-10M \
  --repo-type dataset \
  --local-dir "${ROOT}" \
  --include "seq/part-00000*" \
  --include "item_feat/*" \
  --include "user_feat/*" \
  --include "indexer.pkl"
