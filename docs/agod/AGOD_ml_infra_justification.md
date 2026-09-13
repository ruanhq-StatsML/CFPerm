# AGOD in ML infra — why this tiny Amazon use-case belongs in the system

Amazon Reviews is a **smoke use-case**, not the product. The code we land is a
**control-plane primitive**: online sensors → per-modality LR actuators.

That is an ML-infra concern even when the demo is small.

## What we put in the system

```
agod/
  mmd.py            # sensor: unbiased RBF MMD²(P(X))
  shift.py          # sensor fusion: concept vs covariate
  lr_controller.py  # actuator: Softmax α → LR multipliers (+ EMA state)
  policies.py       # named policies B1–B5 (equal / RF / MMD / hybrid)
scripts/run_agod_amazon_mmd_lr.py   # thin Amazon harness + Acc telemetry
```

Amazon-specific IO (shards, ResNet tower, category stream) stays in
`scripts/run_agod_amazon_modality_lr.py`. The **policy math** lives in `agod/`.

## ML-infra justifications (concrete)

### 1. Adaptation controller ≠ training recipe

Online multimodal systems need a **scheduler that is not hardcoded**.
Hardcoding `lr_text=lr_image` (or a fixed decay) is a model-training choice.
Driving `LR_m` from live shift statistics is a **runtime control service**:

| plane | piece |
|-------|--------|
| sensor | MMD² / PO / RF AUC per modality per window |
| state | EMA α |
| actuator | AdamW param-group LR multipliers |

Same pattern as autoscaling / congestion control: measure → decide → act.
Amazon is just the first client of that API.

### 2. Sensor / actuator split (swap datasets without rewrite)

- **Sensors** (`mmd`, `shift`) only see feature blocks + labels.
- **Actuators** (`lr_controller`) only see modality names + scores.
- **Harness** owns data + model.

So Mind14/COCO/MSR-VTT can reuse the same controller; only the feature
extractor changes. That is the infra reuse story — not “another notebook”.

### 3. Acc-seeking routing is an **update-FLOPs** policy

AGOD’s efficiency claim is **online adaptation compute**, not inference latency.
Covariate↑ → LR↓ (and optionally gate) means fewer effective update steps on
channels that should not be chased. Concept↑ → LR↑ spends adapt budget where
`P(Y|X)` actually moved.

Infra owns **budgeted adaptation**: who gets step size under a fixed step
quota. That is why the controller lives next to training loops, not only in a paper appendix.

### 4. Telemetry contract for CI / dashboards

Each window emits a fixed schema:

```
t, category, cov_m, concept_m, alpha_m, lr_mult_m, acc_pre, acc_post, acc_lift
```

JSON + TeX + PNG under `results/agod_amazon/` are the regression artifacts.
PR smoke can assert `B5_minus_B1_lift > 0` without re-deriving the math.
That is standard ML-infra: **policy → metrics → gate**.

### 5. Why MMD specifically in infra (not just RF)

RF Domain AUC/VIMP is a heavy, saturating classifier sensor (often ~0.5/0.5
after the first jump). Unbiased RBF **MMD²** is:

- **O(n²d)** on a capped subsample (`max_n≈96`) → predictable CPU budget
- non-saturating graded `P(X)` signal → useful damping
- same estimator family as FSDS MMD-LOCO → one sensor stack across attribution
  and adaptation

Infra prefers **bounded-cost sensors** with stable numerical behavior over
“fit another RF every window” when the only job is to damp covariate channels.

### 6. Locked Acc policy (what to ship)

Empirical Amazon smoke (6 category windows):

| Policy | Rule | Mean Acc lift |
|--------|------|--------------:|
| B1 | equal LR | +0.058 |
| B2 | RF concept−cov | +0.072 |
| B3 | pure MMD concept−cov | +0.058 |
| B4 | MMD + intensity gain | +0.037 |
| **B5** | **MMD cov + PO concept** | **+0.087** |

**Ship B5 as default** (`agod.policies` / CLI `--policy` later): MMD for
covariate damping, PO for concept boost. Pure joint-MMD concept and intensity
gain are kept as ablations, not defaults.

## What this is *not*

- Not an inference-latency optimization.
- Not a full feature store / trainer rewrite.
- Not claiming Amazon SOTA — it is a **system-shaped smoke** that proves the
  controller API and Acc direction.

## Run

```bash
PYTHONPATH=. python3 scripts/run_agod_amazon_mmd_lr.py
# or: python3 -c "import agod; print(agod.__version__)"
```
