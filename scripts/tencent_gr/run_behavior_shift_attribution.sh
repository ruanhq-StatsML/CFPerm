#!/usr/bin/env bash
# Behavior / purchase-intent shift attribution — one-command landing run.
# Usage:
#   bash scripts/tencent_gr/run_behavior_shift_attribution.sh
#   bash scripts/tencent_gr/run_behavior_shift_attribution.sh --gt-items path.csv --gt-orders path.csv
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
export PYTHONPATH="${PYTHONPATH:-}:."

ROOT_DATA="${TENCENT_ROOT:-data/tencent_subset}"
OUT="${OUT_DIR:-results/tencent_gr_w1w2_mmd_po_fsds}"
MAX_USERS="${MAX_USERS:-20000}"
GAP_DAYS="${GAP_DAYS:-30}"
LOCALIZE_K="${LOCALIZE_K:-200}"

extra=()
if [[ -n "${GT_ITEMS:-}" ]]; then extra+=(--gt-items "$GT_ITEMS"); fi
if [[ -n "${GT_ORDERS:-}" ]]; then extra+=(--gt-orders "$GT_ORDERS"); fi
# pass-through CLI args (e.g. --gt-items ...)
extra+=("$@")

echo "[behavior-shift] root=$ROOT_DATA out=$OUT max_users=$MAX_USERS gap=$GAP_DAYS k=$LOCALIZE_K"
python3 scripts/tencent_gr/run_w1w2_mmd_po_localize_fsds.py \
  --root "$ROOT_DATA" \
  --out-dir "$OUT" \
  --max-users "$MAX_USERS" \
  --gap-days "$GAP_DAYS" \
  --localize-k "$LOCALIZE_K" \
  --select-k 15 \
  "${extra[@]}"

echo "[behavior-shift] done → $OUT"
ls -la "$OUT" | head -20
