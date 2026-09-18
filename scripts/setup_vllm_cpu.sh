#!/usr/bin/env bash
# Install / verify vLLM CPU in .venv-vllm
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
VENV="${ROOT}/.venv-vllm"
python3 -m venv "${VENV}" 2>/dev/null || true
# shellcheck disable=SC1091
source "${VENV}/bin/activate"
python -m pip install -U pip wheel
python -m pip install 'torch==2.6.0+cpu' 'torchvision==0.21.0+cpu' 'torchaudio==2.6.0+cpu' \
  --index-url https://download.pytorch.org/whl/cpu
python -m pip install 'vllm-cpu==0.10.2' --extra-index-url https://download.pytorch.org/whl/cpu
python - <<'PY'
from vllm import LLM, SamplingParams
import vllm, torch
print("vllm", vllm.__version__, "torch", torch.__version__, "ok")
PY
