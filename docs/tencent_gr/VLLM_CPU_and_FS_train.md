# vLLM CPU + TencentGR 1000-feat → select → train

## vLLM CPU install (this machine)

```bash
# dedicated venv (already created as .venv-vllm)
source .venv-vllm/bin/activate
# or re-run:
bash scripts/setup_vllm_cpu.sh
```

Pinned stack that works here (no GPU):

- `torch==2.6.0+cpu`, `torchvision==0.21.0+cpu`, `torchaudio==2.6.0+cpu`
- `vllm-cpu==0.10.2`

Verify:

```bash
.venv-vllm/bin/python -c "from vllm import LLM, SamplingParams; import vllm, torch; print(vllm.__version__, torch.__version__)"
```

Notes:

- Do **not** install CUDA `torchaudio` into this venv (breaks with missing `libcudart`).
- First `LLM(...)` load will download weights; keep models small on CPU.

## Auto 1000 features → selection → train

```bash
python3 scripts/tencent_gr/auto_feats_select_train.py \
  --root data/tencent_subset \
  --max-users 6000 \
  --target-dim 1000 \
  --select-k 128 \
  --select-method f \
  --label future_cnv \
  --prefix-frac 0.75 \
  --out results/tencent_gr_fs
```

**Label `future_cnv` (recommended):** features from sequence **prefix**, label = whether **suffix** has a conversion. Avoids trivial leakage from `pay_user := has_any_cnv`.

Recipe:

1. Generate ~1000 combinatorial seq features (windows × funnel × decay × session × attribution × ARPU/deal × Markov × crosses × user side).
2. `VarianceThreshold` → `SelectKBest(f_classif|MI)` top-k.
3. Train `HistGradientBoosting` + `LogisticRegression`; report AUC / AP / Acc.
