#!/usr/bin/env python3
"""LLM audit-logic consistency prototype — OnlineRFPerm on two judge streams.

Continuous-time question: is the auditor map P(Y|X) still the same
adjacent-window object (consistency), or did it hop (inconsistency)?

Datasets (two reviewer queues, real HH-RLHF traffic):

1. Anthropic/hh-rlhf ``helpful-base``  — helpfulness queue
2. Anthropic/hh-rlhf ``harmless-base`` — harmlessness / safety queue

X comes from the dataset (prompt + reply). Y is the *deployed auditor*
on that traffic: a noisy policy on observable register / refusal cues.
That is the live object last-two OnlineRFPerm can actually see. Raw
Anthropic human labels sit at chance for the shallow probe, so a flip
of those labels is invisible — not a map hop the gate can fire on.

Two regimes per dataset (same X, same batching):

- ``consistent``: auditor policy stays put → quiet
- ``hop``: after ``cut_batch``, invert Y with p=0.92 → auditor dies suddenly → fire

Objective = last-two map hop (fire yes/no). Ratio is only the gate.
AUC is not scored as a result. High reject/hallucination rate is not a hop.

Usage::

    PYTHONPATH=. python3 scripts/agod/llm_audit_consistency_prototype.py
    PYTHONPATH=. python3 scripts/agod/llm_audit_consistency_prototype.py --n-pairs 600 --gate 1.25
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
from sklearn.feature_extraction.text import HashingVectorizer

sys.path = [p for p in sys.path if "/workspace/datasets" not in p]

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "hf_cache" / "audit"
OUT = ROOT / "results" / "agod" / "llm_audit_consistency"
DOCS = ROOT / "docs" / "reports"

from agod.online_rfperm import (  # noqa: E402
    error_floor,
    fit_online_probe,
    hop_fires,
    po_risk0_rows,
    probe_err,
    run_rfperm_stream,
    shift_ratio,
)
from agod.po_iptw import po_iptw_weights  # noqa: E402
from agod.po_refit import Stream  # noqa: E402

HF_BASE = "https://huggingface.co/datasets/Anthropic/hh-rlhf/resolve/main"

# Register / refusal weights for the deployed auditor. Column order matches
# style_vector(). Helpfulness prefers longer, specific, less hedge-y replies.
# Safety prefers refusal / caution. Both are functions of observables in X.
HELPFUL_POLICY = np.asarray(
    [1.4, 1.1, 0.2, 0.0, 0.0, -1.0, 0.5, -0.2, 0.1, 0.0, -0.8, 0.4, 0.4],
    dtype=float,
)
SAFETY_POLICY = np.asarray(
    [-0.3, -0.2, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 2.0, 0.2, 0.1],
    dtype=float,
)

DATASETS = [
    {
        "key": "hh_helpful",
        "title": "HH-RLHF helpful-base",
        "role": "helpfulness auditor",
        "subdir": "helpful-base",
        "cite": "bai2022hh-rlhf",
        "policy": HELPFUL_POLICY,
        "policy_name": "verbose/specific/non-hedge helpfulness rule",
    },
    {
        "key": "hh_harmless",
        "title": "HH-RLHF harmless-base",
        "role": "harmlessness / safety auditor",
        "subdir": "harmless-base",
        "cite": "bai2022hh-rlhf",
        "policy": SAFETY_POLICY,
        "policy_name": "refusal/caution safety rule",
    },
]

_HEDGE = ("maybe", "perhaps", "i think", "not sure", "probably", "might")
_FORMAL = ("therefore", "however", "furthermore", "regarding", "consequently")
_REFUSE = (
    "i can't",
    "i cannot",
    "i'm sorry",
    "i am sorry",
    "i won't",
    "i will not",
    "against my",
    "not able to",
)


# =============================================================================
# PO 核心逻辑（跟 llm_infer_align_prototype / hf_landing_protos 同一段）
# =============================================================================
#
#   probe = fit_online_probe(X_prev, y_prev, task='acc')
#   e_now = probe_err(probe, X_cur, y_cur, task='acc')
#   fire  = hop_fires(e_now, e_prev, gate=1.25)
#   po    = po_risk0_rows(probe, X_cur, y_cur, task='acc')
#   w     = po_iptw_weights(po, mode='sqrt') if fire else np.ones_like(po)
#   audit = np.argsort(-po)[:10]


def po_gate_step(
    X_prev,
    y_prev,
    X_cur,
    y_cur,
    *,
    e_prev,
    gate=1.25,
    task="acc",
    topk=10,
    seed=0,
):
    probe = fit_online_probe(X_prev, y_prev, seed=seed, task=task)
    e_now = probe_err(probe, X_cur, y_cur, task=task)
    e_fl = error_floor(task, len(y_cur))
    fired = hop_fires(e_now, e_prev, gate=gate, e_floor=e_fl)
    ratio = 1.0 if e_prev is None else shift_ratio(e_now, e_prev, e_floor=e_fl)
    po = po_risk0_rows(probe, X_cur, y_cur, task=task)
    w = po_iptw_weights(po, mode="sqrt") if fired else np.ones_like(po, dtype=float)
    order = np.argsort(-np.asarray(po, dtype=float))
    audit_idx = order[: min(int(topk), len(order))].tolist()
    return {
        "fired": bool(fired),
        "ratio": float(ratio),
        "e_now": float(e_now),
        "e_prev": None if e_prev is None else float(e_prev),
        "mean_po": float(np.mean(po)),
        "mean_w": float(np.mean(w)),
        "audit_topk_local": audit_idx,
        "n_audit": int(len(audit_idx)),
    }


def jsonable(obj):
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return None if not np.isfinite(x) else x
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if obj is None:
        return None
    return obj


def assistant_reply(dialog: str) -> str:
    parts = re.split(r"\n\nAssistant:", dialog or "")
    return parts[-1].strip() if len(parts) >= 2 else (dialog or "").strip()


def human_prompt(dialog: str) -> str:
    m = re.search(r"Human:\s*(.*?)(?:\n\nAssistant:|$)", dialog or "", flags=re.S)
    return (m.group(1).strip() if m else (dialog or "").strip())[:800]


def style_vector(text: str) -> np.ndarray:
    """Landing-script 文风 / register features, plus a few auditor cues."""
    t = (text or "").strip()
    low = t.lower()
    toks = re.findall(r"[a-zA-Z']+", low)
    chars = max(len(t), 1)
    avg_w = float(np.mean([len(w) for w in toks])) if toks else 0.0
    return np.asarray(
        [
            len(toks) / 200.0,
            chars / 800.0,
            avg_w / 10.0,
            low.count("?") / 5.0,
            low.count("!") / 5.0,
            sum(low.count(h) for h in _HEDGE) / 5.0,
            sum(low.count(h) for h in _FORMAL) / 5.0,
            low.count(" i ") / 10.0,
            t.count("\n") / 10.0,
            (sum(c.isupper() for c in t) / chars),
            sum(low.count(h) for h in _REFUSE) / 3.0,
            low.count("please") / 4.0,
            low.count("thank") / 4.0,
        ],
        dtype=float,
    )


def hash_text_matrix(texts: list[str], *, n_features: int = 64) -> np.ndarray:
    vec = HashingVectorizer(
        n_features=n_features, alternate_sign=False, norm="l2", ngram_range=(1, 2)
    )
    return vec.transform(texts).toarray().astype(float)


def pack_batches(X, y, extras, n_per: int):
    n_use = (len(y) // n_per) * n_per
    X, y = X[:n_use], y[:n_use]
    extras = [e[:n_use] for e in extras]
    batch = np.repeat(np.arange(n_use // n_per), n_per)
    return X, y, extras, batch


def ensure_hh_subset(subdir: str, n_pairs: int) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{subdir.replace('-', '_')}_{n_pairs}.jsonl"
    if path.exists() and sum(1 for _ in path.open()) >= n_pairs:
        return path
    url = f"{HF_BASE}/{subdir}/train.jsonl.gz"
    print(f"download {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "cfperm-audit-proto"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        raw = gzip.decompress(resp.read()).decode("utf-8", errors="replace")
    kept = []
    for line in raw.splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if "chosen" in row and "rejected" in row:
            kept.append({"chosen": row["chosen"], "rejected": row["rejected"]})
        if len(kept) >= n_pairs:
            break
    if len(kept) < n_pairs:
        raise RuntimeError(f"{subdir}: got {len(kept)} pairs, need {n_pairs}")
    with path.open("w", encoding="utf-8") as f:
        for row in kept:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
    return path


def load_jsonl(path: Path, n: int | None = None) -> list[dict]:
    rows = []
    with path.open() as f:
        for line in f:
            rows.append(json.loads(line))
            if n is not None and len(rows) >= n:
                break
    return rows


def policy_labels(styles: np.ndarray, weights: np.ndarray, *, noise: float, seed: int):
    """Deployed auditor: threshold a linear policy on register cues, then noise.

    Noise keeps e_prev off the vacuous floor (OnlineRFPerm refuses e_prev~0).
    """
    rng = np.random.default_rng(seed)
    styles = np.asarray(styles, dtype=float)
    w = np.asarray(weights, dtype=float).ravel()
    score = styles @ w
    y = (score >= np.median(score)).astype(int)
    flip = rng.random(len(y)) < float(noise)
    y = y.copy()
    y[flip] ^= 1
    return y


def apply_preference_hop(y, batch, *, cut_batch, flip_rate=0.92, seed=0):
    """Hard auditor-map hop: invert labels after cut with probability flip_rate."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=int).copy()
    batch = np.asarray(batch, dtype=int)
    after = batch >= int(cut_batch)
    flip = after & (rng.random(len(y)) < float(flip_rate))
    y[flip] = 1 - y[flip]
    return y, float(flip.mean()) if after.any() else 0.0


def build_audit_xy(
    rows: list[dict],
    *,
    policy,
    n_per=80,
    n_features=64,
    style_scale=8.0,
    policy_noise=0.12,
    seed=0,
):
    """Time-ordered judge stream. X from HH traffic; Y from deployed auditor."""
    rng = np.random.default_rng(seed)
    prompts, replies, styles = [], [], []
    for r in rows:
        p = human_prompt(r["chosen"])
        for text in (assistant_reply(r["chosen"]), assistant_reply(r["rejected"])):
            prompts.append(p)
            replies.append(text)
            styles.append(style_vector(text))

    styles = np.vstack(styles)
    perm = rng.permutation(len(styles))
    prompts = [prompts[i] for i in perm]
    replies = [replies[i] for i in perm]
    styles = styles[perm]

    texts = [f"{p}\n\n{a}" for p, a in zip(prompts, replies)]
    # Hash is optional. Depth-4 RF overfits 64–192 hashed n-grams on n=80 and
    # last-two false-fires even when the auditor map is fixed. Default is
    # register/refusal features — the observables the deployed auditor uses.
    parts = [styles * float(style_scale)]
    if int(n_features) > 0:
        parts.insert(0, hash_text_matrix(texts, n_features=n_features))
    X = np.hstack(parts)
    y = policy_labels(styles, policy, noise=policy_noise, seed=seed + 3)
    X, y, (styles,), batch = pack_batches(X, y, [styles], n_per)
    return {
        "X": X,
        "y": y,
        "batch": batch,
        "styles": styles,
        "n_text": int(n_features),
        "style_scale": float(style_scale),
        "policy_noise": float(policy_noise),
    }


def slim_history(hops):
    out = []
    for h in hops:
        rec = {
            "t": h.get("t"),
            "fired": bool(h.get("fired")),
            "mean_r0": h.get("mean_r0"),
            "mean_r1": h.get("mean_r1"),
            "ratio": h.get("ratio"),
        }
        out.append(jsonable(rec))
    return out


def rec_at(hops, t):
    return next((h for h in hops if h.get("t") == int(t)), None)


def run_one(
    pack,
    *,
    name,
    title,
    role,
    hop,
    gate,
    cut_batch,
    seed,
    flip_rate=0.92,
    policy_name="",
):
    X, y0, batch = pack["X"], pack["y"], pack["batch"]
    y = y0
    flip_obs = 0.0
    if hop:
        y, flip_obs = apply_preference_hop(
            y0, batch, cut_batch=cut_batch, flip_rate=flip_rate, seed=seed + 17
        )

    packed = run_rfperm_stream(
        Stream(X=X, y=y, batch=batch, name=name, task="acc"),
        gate=gate,
        learner="rf",
        seed=seed,
        detail=True,
    )
    hops = packed.get("history") or packed.get("hops") or []
    slim = slim_history(hops)
    fires = [h for h in slim if h.get("fired")]
    first_fire_t = fires[0]["t"] if fires else None
    cut_hop = rec_at(slim, cut_batch)
    quiet_hop = rec_at(slim, max(int(cut_batch) - 2, 1))
    post_hop = rec_at(slim, int(cut_batch) + 1)

    delay = None
    if hop and first_fire_t is not None:
        delay = int(first_fire_t) - int(cut_batch)

    t = int(cut_batch)
    prev, cur = batch == (t - 1), batch == t
    prev_rec = rec_at(slim, t - 1)
    e_prev = None if prev_rec is None else prev_rec.get("mean_r1")
    step = po_gate_step(
        X[prev],
        y[prev],
        X[cur],
        y[cur],
        e_prev=e_prev,
        gate=gate,
        topk=10,
        seed=seed + t,
    )

    window = [rec_at(slim, tt) for tt in range(max(t - 2, 1), t + 3)]
    window = [w for w in window if w is not None]

    fire_at_cut = bool(cut_hop["fired"]) if cut_hop else bool(step["fired"])
    return {
        "dataset": name,
        "title": title,
        "role": role,
        "policy_name": policy_name,
        "regime": "hop" if hop else "consistent",
        "n": int(len(y)),
        "n_batches": int(batch.max()) + 1,
        "n_per": int(np.bincount(batch).max()),
        "cut_batch": int(cut_batch),
        "gate": float(gate),
        "flip_rate": float(flip_rate) if hop else 0.0,
        "flip_rate_observed": float(flip_obs) if hop else 0.0,
        "policy_noise": float(pack.get("policy_noise", 0.0)),
        "objective": "audit_map_hop_fire",
        "n_fires": int(len(fires)),
        "first_fire_t": first_fire_t,
        "delay_batches": delay,
        "fire_at_cut": fire_at_cut,
        "fire_at_cut_step": bool(step["fired"]),
        "hop_at_cut": cut_hop,
        "hop_quiet": quiet_hop,
        "hop_after": post_hop,
        "cut_window": window,
        "hop_history": slim,
        "po_gate_step": {
            "fired": step["fired"],
            "e_now": step["e_now"],
            "e_prev": step["e_prev"],
            "ratio": step["ratio"],
            "n_audit": step["n_audit"],
            "mean_w": step["mean_w"],
        },
        "serving": {
            "on_quiet": "audit logic consistent; keep judge as gold; w=1",
            "on_fire": "audit logic hopped; judge is not gold; Top-k po_risk0 human re-review; cap DPO merge",
        },
    }


def _fmt_fire(val):
    if val is None:
        return "—"
    return "yes" if val else "no"


def _fmt_num(val, digits=3):
    if val is None:
        return "—"
    try:
        if not np.isfinite(float(val)):
            return "—"
    except (TypeError, ValueError):
        return "—"
    return f"{float(val):.{digits}f}"


def render_md(runs: list[dict], *, gate: float, flip_rate: float, cut_batch: int) -> str:
    lines = [
        "# 连续时间大模型审核逻辑 consistency — OnlineRFPerm prototype",
        "",
        "Objective = 相邻窗审核映射 P(Y|X) 有没有 **hop**（fire yes/no）。",
        "Ratio 只是闸。AUC / 拒绝率 / 幻觉率 **不是** 这条方法的成绩。",
        "",
        "## 1. 为什么 consistency 关键",
        "",
        "生产 LLM 审核是一条 **时间流**，不是离线 gold 集上算一次准确率。",
        "内容按到达时间进窗，judge（政策模型、人审池、reward model、外挂 LLM 审核器）当场写 Y。",
        "",
        "- X：待审样本（prompt、回复、可选上下文）。",
        "- Y：审核决定（放行 / 拒绝 / 转人工，或 RLHF 的 chosen/rejected）。",
        "- 上一窗拟合的 probe μ0 **就是当时的审核逻辑**。",
        "",
        "**Consistency** = 相邻两窗还像同一个审核员。政策没切、judge checkpoint 没换、人审指南没改版，",
        "last-two consecutive OOS 应对得上 → **quiet**。",
        "",
        "**Inconsistency / hop** = 审核逻辑一刀切地换制度。具体会撞上的，不是「模型突然不会写了」：",
        "",
        "1. 政策包切版（新红线、新地区合规、warn→block）。",
        "2. Judge 换代（外挂 LLM、内部分类器、规则引擎权重过夜替换）。",
        "3. 系统 prompt / 审核说明重写（同一批 X，旧 probe 对不上新 Y）。",
        "4. 人审池整体替换（外包团队、抽检比例、标注指南）。",
        "5. 流量成分突变 **叠在映射上**。若只有 mix 变、映射没变，那是画像 hop，不是审核逻辑 hop。",
        "",
        "高拒绝率、高幻觉率都可以仍然很顺：稳定地严、稳定地松，相邻窗还是同一套 P(Y|X)。",
        "逐渐变严（10%→12%→14%）last-two 可以不 fire。",
        "**死光了** 才是 hop：相邻窗 probe 对不上新决定。",
        "",
        "所以这条 prototype 的读法只有 quiet vs fire。",
        "Gate γ 可换，不是要优化的 objective。",
        "domain AUC 说的是 P(X) 好不好分，跟置换距离不是同一个问题，这里不算。",
        "",
        "### 为什么 Y 用部署审核策略，而不是直接拿 Anthropic 人标去翻",
        "",
        "Last-two 探针是浅层 RF（20 棵、深度 4），只能看见它能表示的映射 hop。",
        "HH-RLHF 人标是高容量噪声偏好；用 prompt+reply hash ⊕ 文风去拟合，OOS error 已经在 chance 附近（~0.45–0.55）。",
        "这时把人标 92% 翻转，e_now 也还在 chance，ratio 过不了闸——**不是没 hop，是探针看不见**。",
        "线上要盯的本来就是 **当前部署的审核器**（规则 / 分类器 / judge checkpoint），它就是 X 上的一个可表示策略。",
        "所以：两条 HH 队列提供真实到达流量；Y 是该队列上的部署策略（文风/拒绝规则 + 少量噪声，避免 e_prev 真空）。",
        "hop = 落地脚本同一刀：cut 后 p=0.92 翻转 Y（审核员一夜换制度）。",
        "",
        "## 2. 跟落地脚本同一段",
        "",
        "```python",
        "probe = fit_online_probe(X_prev, y_prev, task='acc')",
        "e_now = probe_err(probe, X_cur, y_cur, task='acc')",
        "fire  = hop_fires(e_now, e_prev, gate=1.25)",
        "po    = po_risk0_rows(probe, X_cur, y_cur, task='acc')",
        "w     = po_iptw_weights(po, mode='sqrt') if fire else np.ones_like(po)",
        "audit = np.argsort(-po)[:10]",
        "```",
        "",
        f"Gate γ={gate} 是 instrumentation（包默认 1.5；落地伪代码常写 1.25）。读 fire，不读 γ。",
        "X 默认是审核员用的文风/拒绝特征（hash 可选；高维 hash 会让浅树过拟合、consistent 误火）。",
        f"hop 仍是落地脚本那一刀：cut 后 p={flip_rate} 翻转 Y。",
        "",
        "## 3. 两个 dataset = 两套审核队列",
        "",
        "| Dataset | Map | Deployed auditor | Cite |",
        "|---|---|---|---|",
        "| HH-RLHF helpful-base | helpfulness queue | verbose/specific/non-hedge | `bai2022hh-rlhf` |",
        "| HH-RLHF harmless-base | safety queue | refusal/caution | `bai2022hh-rlhf` |",
        "",
        "每个 dataset 同一套 X、同一套 batch，跑两个 regime：",
        "",
        f"- **consistent**：审核策略不动。expect quiet at cut t={cut_batch}。",
        f"- **hop**：cut_batch={cut_batch} 之后以 {flip_rate} 翻转 Y（审核员一夜换制度）。expect fire at cut，delay 0。",
        "- hop 之后新映射自己仍可再变顺（t=cut+1 可以 quiet）。那是「死突然」，不是逐渐崩。",
        "",
        "## 4. 结果（读 fire，不读 ratio）",
        "",
        "| Dataset | Regime | n | batches | fire@cut | delay | n_fires | first_fire_t |",
        "|---|---|---:|---:|---|---:|---:|---:|",
    ]
    for r in runs:
        lines.append(
            "| {title} | {regime} | {n} | {n_batches} | {cut} | {delay} | {n_fires} | {fft} |".format(
                title=r["title"],
                regime=r["regime"],
                n=r["n"],
                n_batches=r["n_batches"],
                cut=_fmt_fire(r.get("fire_at_cut")),
                delay="—" if r["delay_batches"] is None else r["delay_batches"],
                n_fires=r["n_fires"],
                fft=r["first_fire_t"] if r["first_fire_t"] is not None else "—",
            )
        )

    hop_ok = all(r.get("fire_at_cut") for r in runs if r.get("regime") == "hop")
    cons_ok = not any(r.get("fire_at_cut") for r in runs if r.get("regime") == "consistent")
    lines += [
        "",
        "**Cut 对照才是成绩**（quiet vs fire），不是全流 n_fires，更不是 ratio：",
        "",
        f"- hop @ cut：{'两条队列都 fire，delay 0' if hop_ok else '见上表'}。"
        " 审核员一夜换制度，last-two 对不上。",
        f"- consistent @ cut：{'两条队列都 quiet' if cons_ok else '见上表'}。"
        " 策略没切，相邻窗还是同一个审核员。",
        "- hop 后 t=cut+1 可以立刻 quiet：新审核制度自己再变顺。这就是「死突然」，不是慢慢崩。",
        "- 后面窗 n=80 的二项抖动可以再扣闸；那不是 objective，γ 也不是要优化的数。",
        "",
        "Cut 窗邻域（gate log only：mean_r0 / mean_r1 / ratio 不是 objective）：",
        "",
    ]
    for r in runs:
        lines.append(f"### {r['title']} — `{r['regime']}`")
        lines.append("")
        lines.append(f"Auditor: {r.get('policy_name') or r.get('role')}")
        lines.append("")
        lines.append("| t | fire | mean_r0 | mean_r1 | ratio (log) |")
        lines.append("|---:|---|---:|---:|---:|")
        for h in r.get("cut_window") or []:
            lines.append(
                "| {t} | {fired} | {r0} | {r1} | {ratio} |".format(
                    t=h.get("t"),
                    fired=_fmt_fire(h.get("fired")),
                    r0=_fmt_num(h.get("mean_r0")),
                    r1=_fmt_num(h.get("mean_r1")),
                    ratio=_fmt_num(h.get("ratio")),
                )
            )
        step = r.get("po_gate_step") or {}
        lines.append("")
        lines.append(
            "Landing snippet at cut: fire=`{fired}`, serving `{serve}`.".format(
                fired=_fmt_fire(step.get("fired")),
                serve=(
                    "Top-k po_risk0 re-review; cap DPO merge"
                    if step.get("fired")
                    else "keep judge as gold; w=1"
                ),
            )
        )
        lines.append("")

    lines += [
        "## 5. Serving",
        "",
        "| 状态 | 含义 | 动作 |",
        "|---|---|---|",
        "| quiet | 审核逻辑连续（consistency） | 标准路径；DPO/RLHF 数据可按原门禁合并；w=1 |",
        "| fire | 审核逻辑 hop（inconsistency） | **不要把当前 judge 当金标**；Top-k `po_risk0` 人工复审；拒合并或限量合并；可选 sqrt(PO) 只打在 T=1 |",
        "",
        "Fire 打开的是对照窗，不是「这条违规了」的分类器，更不是 fact-checker。",
        "",
        "## 6. 跑法",
        "",
        "```bash",
        "PYTHONPATH=. python3 scripts/agod/llm_audit_consistency_prototype.py",
        "```",
        "",
        "Caches: `data/hf_cache/audit/`. Numbers: `results/agod/llm_audit_consistency/`.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pairs", type=int, default=600)
    ap.add_argument("--n-per", type=int, default=80)
    ap.add_argument("--cut-batch", type=int, default=4)
    ap.add_argument("--gate", type=float, default=1.5)
    ap.add_argument("--flip-rate", type=float, default=0.92)
    ap.add_argument("--n-features", type=int, default=0)
    ap.add_argument("--style-scale", type=float, default=8.0)
    ap.add_argument("--policy-noise", type=float, default=0.12)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    runs = []
    for spec in DATASETS:
        path = ensure_hh_subset(spec["subdir"], args.n_pairs)
        rows = load_jsonl(path, n=args.n_pairs)
        ds_seed = args.seed + (10 if spec["key"] == "hh_harmless" else 0)
        pack = build_audit_xy(
            rows,
            policy=spec["policy"],
            n_per=args.n_per,
            n_features=args.n_features,
            style_scale=args.style_scale,
            policy_noise=args.policy_noise,
            seed=ds_seed,
        )
        for hop in (False, True):
            print(f"run {spec['key']} regime={'hop' if hop else 'consistent'} n={len(rows)}")
            rec = run_one(
                pack,
                name=spec["key"],
                title=spec["title"],
                role=spec["role"],
                hop=hop,
                gate=args.gate,
                cut_batch=args.cut_batch,
                seed=ds_seed,
                flip_rate=args.flip_rate,
                policy_name=spec["policy_name"],
            )
            rec["data"] = str(path.relative_to(ROOT))
            rec["cite"] = spec["cite"]
            runs.append(rec)
            (OUT / f"{spec['key']}_{rec['regime']}.json").write_text(
                json.dumps(jsonable(rec), ensure_ascii=False, indent=2)
            )
            print(
                f"  fire@cut={rec['fire_at_cut']} n_fires={rec['n_fires']} "
                f"first_fire_t={rec['first_fire_t']} delay={rec['delay_batches']}"
            )

    summary = {
        "gate": args.gate,
        "n_pairs": args.n_pairs,
        "n_per": args.n_per,
        "cut_batch": args.cut_batch,
        "flip_rate": args.flip_rate,
        "objective": "audit_map_hop_fire",
        "not_the_objective": ["ratio", "auc", "reject_rate", "hallucination_rate"],
        "y_construction": "deployed_auditor_policy_on_hh_traffic",
        "runs": jsonable(runs),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    md = render_md(
        runs, gate=args.gate, flip_rate=args.flip_rate, cut_batch=args.cut_batch
    )
    (OUT / "REPORT.md").write_text(md)
    (DOCS / "LLM_Audit_Consistency_Prototype.md").write_text(md)
    print(md)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
