# Stable Diffusion evaluation — open the model, where it gets tricky

This note opens **SD 1.5** (`StableDiffusionPipeline`) and maps it onto our
DiffusionDB temporal-attribution combo (`X` prompt tokens → `Y` image
attribute → `T` time). Goal: say clearly what is *easy to narrate* vs what is
*evaluation-tricky*.

Inspected configs (no full weights):
`stable-diffusion-v1-5/stable-diffusion-v1-5` → `model_index.json` +
`text_encoder` / `unet` / `vae` / `scheduler` / `safety_checker`.

---

## 1. Open the pipeline (what actually runs)

```text
prompt
  → CLIPTokenizer          (max length 77)
  → CLIPTextModel          hidden (77, 768)     ← true token X in SD
  → UNet2DConditionModel   latent denoise
        ↑ cross-attn on text (cross_attention_dim=768)
        ↑ CFG: ε_θ(x|∅) vs ε_θ(x|c), scale = cfg
  → AutoencoderKL decode   → 512×512 RGB
  → SafetyChecker (CLIP)   → NSFW flag / filter   ← related to image_nsfw
```

| Block | Role | Shape / knobs |
|---|---|---|
| **Tokenizer** | BPE → ids | length ≤ **77** |
| **Text encoder** | CLIP text transformer | **(77, 768)**, 12 layers, vocab 49408 |
| **UNet** | conditional denoise in latent | latent 4×64×64; cross-attn 768 |
| **Scheduler** | noise schedule | train T=1000; inference `step` |
| **VAE** | latent ↔ RGB | out 512×512 |
| **Safety checker** | CLIP-based NSFW screen | projection_dim 768 |

So inside SD, **the real prompt embedding is already CLIP `(77, 768)`** —
exactly Definition 1’s `X`. Our TF-IDF `d≤512` is a *light surrogate* of that
token channel, not the UNet’s conditioning.

Hyperparameters that **also** enter the image (and thus any `Y`):

- `seed`, `step` (denoise steps), `cfg` (guidance scale), `sampler` / scheduler
- (implicit) which checkpoint: SD1.4 / 1.5 / 2.x / fine-tunes

DiffusionDB metadata stores many of these next to the prompt.

---

## 2. What “evaluation / attribution Y” can mean (three layers)

| Layer | `Y` example | Lives where | Reads as |
|---|---|---|---|
| **A. Safety / filter** | `image_nsfw` | SafetyChecker / LAION NSFW head on **pixels** | “how NSFW the *image* looks” |
| **B. Alignment** | CLIP score \(\langle e_{\mathrm{img}}, e_{\mathrm{txt}}\rangle\) | external CLIP on (image, prompt) | “does image match prompt?” |
| **C. Quality / style** | aesthetic 1–10, style logits | external aesthetic / attribute head | “how pretty / painterly / …” |

Our current smoke uses **Layer A** (`image_nsfw` in metadata) — convenient
(no PNG), but it is **not** CLIP alignment and **not** aesthetic quality.

---

## 3. Where it gets tricky (the important part)

### T1 — `Y` is not “the SD model score”
SD does not spit out a single scalar “quality”. `Y` is always an **external
or side-car head**. Mixing them changes the scientific question:

- NSFW ↑ over time ≠ prompts got “better”
- CLIP score ↑ ≠ aesthetics ↑
- Aesthetic ↑ with more `artstation` tokens may be **covariate** (prompt mix)
  not **concept** (same prompt → different look)

FSDS on tokens → `Y` answers: *which tokens associate with high `Y` under
shift* — only as honest as the chosen `Y`.

### T2 — Confounders sit next to the prompt
Even with fixed prompt text, image (and `Y`) moves with `cfg` / `step` /
`sampler` / `seed` / checkpoint. Temporal attribution that **only** puts
prompt tokens in `X` will **absorb hyperparam drift into “token effects”**
unless you:

- condition / stratify on `(cfg, step, sampler)`, or
- put those knobs into `X` as controls, or
- restrict to a fixed hyperparam slice

DiffusionDB rows are **user-actuated**: people change prompts *and* knobs
together over time. That is the main evaluation pitfall.

### T3 — SafetyChecker ≠ text encoder
`image_nsfw` is closer to a **vision CLIP safety head** on the decoded image.
`prompt_nsfw` is text-side. Our X is prompt tokens; predicting `image_nsfw`
is a **cross-modal** link (text tokens → image safety score). Intuitive, but:

- not the same as “which tokens the UNet attended to”
- not causal for generation; descriptive association under `(P_t, Y_t)`

### T4 — CLIP score is circular if misused
If `Y = CLIP(image, prompt)` and `X = CLIP text tokens`, early layers of the
story share the same embedding family. Still valid for *alignment drift*, but
you must say so — it is not an independent aesthetic outcome.

### T5 — Token position vs bag-of-tokens
True SD conditioning is **ordered** `(77, 768)` with BOS/EOS/pad. TF-IDF
bag-of-1–2grams loses position and prompt grammar (`a red car` vs
`car, red, …`). Good for “which words rose”, weak for “which slot in the 77
drove UNet cross-attn”.

### T6 — Time window vs model version
DiffusionDB 2M metadata span is short (~Aug 2022). “Late vs early” is mostly
**user-mix / prompt-fashion** within one gallery era, not SD1.5→SDXL.
SD Image Attribution (model identity as `Y`) is a **different** task:
checkpoint provenance, not temporal prompt drift.

---

## 4. Map onto our combo (what we did vs what SD actually conditions on)

| Piece | Our smoke | SD-native / stricter eval |
|---|---|---|
| `X` | TF-IDF tokens | CLIPText `(77,768)` from **same** tokenizer/encoder as the generator |
| `Y` | `image_nsfw` | pick one layer A/B/C and freeze it; report which |
| `T` | early/late count windows | calendar bins **or** checkpoint eras — don’t mix silently |
| Controls | none yet | `cfg`, `step`, `sampler` in `X` or strata |
| Drift views | Domain VIMP + cmean Δ + FSDS blend | same; optionally add PO-risk on `(X,T)→Y` |
| Direction | \(\mathrm{sign}(\Delta\bar Y)\), \(\mathrm{sign}(\delta_j)\) | keep — still required |

**Evaluation logic that stays clean:**

1. Fix the SD text pathway definition of `X` (light or CLIP).
2. Fix one `Y` head and name its layer (A/B/C).
3. Either control generation knobs or put them in `X`.
4. Run the same combo: covariate VIMP ∥ tip-cmean sign ∥ FSDS (and optional confirm timing later).
5. Never read FSDS tops as “UNet attention” without an attention / LOCO probe on the actual encoder.

---

## 5. One-sentence takeaway

Opening SD shows: **tokens enter only through CLIP text `(77,768)` + CFG**;
**every scalar `Y` is an add-on head**; temporal FSDS on prompts is intuitive
**only after** you fix which `Y`, and stop letting `cfg/step/sampler` pretend
to be “token drift”.

Our current DiffusionDB smoke is a valid **Layer-A + light-X** prototype;
the tricky upgrade path is CLIP-text `X` + controlled knobs + an explicit
Layer-B/C `Y`, not a different blend recipe.
