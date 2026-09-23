#!/usr/bin/env python3
"""Iterate PO-risk metric versions → next-window training actuators.

Compares metric versions on Affec / Food-101 / Fashion-IQ / COCO (when present):
  equal | po_soft | po_minus_cov | po_gated | po_proto | po_delta | po_budget | po_next | po_fuse

Long⊗short concept emphasis: see docs/agod/AGOD_PO_Risk_PostTraining.tex §fuse
(freeze←long, step dump←short, LR←fused α).

Actuators applied next window:
  - modality LR multipliers from α
  - step budget allocation
  - optional freeze of low-α modalities
  - stack prior = α (KL on stacking weights)

Efficiency signal: relative online-update FLOPs proxy.
Default ``step_mode=per_mod`` (R2): ``step_alloc[m]`` is real optimizer
steps on modality ``m`` (block schedule, highest dump first).
Legacy ``step_mode=shared`` only scales a shared loop length.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn.functional as F

from agod.mmd import rbf_mmd2
from agod.online_stack import StackFusion, alpha_stack_kl
from agod.po_risk_train import (
    METRIC_VERSIONS,
    NextStepActuatorConfig,
    PORiskMetricConfig,
    continuous_gain_metrics,
    expand_step_schedule,
    freeze_flops_rel,
    metric_to_alpha,
    next_step_actuators,
    next_step_actuators_fused,
    opportunity_rank,
    pick_default_schedule_card,
    po_iptw_weights,
    realize_step_alloc,
    row_po_residual,
    step_flops_rel,
    stream_reject_proxy,
)
from agod.po_roi_export import export_roi_sqlite, rows_to_roi_records, write_roi_jsonl
from agod.proto_drift import ModalityPrototypeBank
from agod.reject_event import resolve_stream_reject
from agod.shift import residual_concept


# ---------------------------------------------------------------------------
# data helpers
# ---------------------------------------------------------------------------


def _pack_windows(
    feats: Dict[str, np.ndarray],
    y: np.ndarray,
    *,
    window: int,
    n_windows: int,
    seed: int,
) -> Tuple[List[Dict[str, np.ndarray]], List[np.ndarray]]:
    n = int(y.shape[0])
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    need = min(n_windows * window, (n // window) * window)
    order = order[:need]
    xs, ys = [], []
    for i in range(0, need, window):
        idx = order[i : i + window]
        xs.append({m: feats[m][idx] for m in feats})
        ys.append(y[idx])
    return xs, ys


def _load_affec(root: Path, max_n: int, seed: int):
    cache = root / "results/affec_fsds/affec_fsds_xyw_cache.npz"
    if not cache.is_file():
        return None
    z = np.load(cache, allow_pickle=True)
    X, Y = z["X"], z["Y"]
    mods = [str(m) for m in z["mods"].tolist()]
    raw_slices = z["block_slices"]
    # cache stores either (M,2) array or 0-d object dict {mod: (lo,hi)}
    if raw_slices.shape == ():
        sl_map = dict(raw_slices.item())
        slices = {m: (int(sl_map[m][0]), int(sl_map[m][1])) for m in mods}
    else:
        slices = {
            mods[i]: (int(raw_slices[i, 0]), int(raw_slices[i, 1]))
            for i in range(len(mods))
        }
    n = min(max_n, X.shape[0])
    rng = np.random.default_rng(seed)
    idx = rng.choice(X.shape[0], size=n, replace=False)
    X, Y = X[idx], Y[idx]
    feats = {
        m: X[:, slices[m][0] : slices[m][1]].astype(np.float32) for m in mods
    }
    y_bin = (Y > np.median(Y)).astype(np.int64)
    return feats, y_bin, mods, "affec"


def _load_img_txt(root: Path, name: str, max_n: int, seed: int):
    d = root / "data/img_txt" / name
    if not (d / "img_feats.npy").is_file():
        return None
    img = np.load(d / "img_feats.npy").astype(np.float32)
    txt = np.load(d / "txt_feats.npy").astype(np.float32)
    y_raw = np.load(d / "labels.npy").astype(np.int64)
    # multi-class packs → binary even/odd (mode-vs-rest too imbalanced)
    y = (y_raw % 2).astype(np.int64)
    n = min(max_n, img.shape[0])
    rng = np.random.default_rng(seed)
    idx = rng.choice(img.shape[0], size=n, replace=False)
    feats = {"img": img[idx], "txt": txt[idx]}
    return feats, y[idx], ["img", "txt"], name


def load_dataset(root: Path, name: str, max_n: int, seed: int):
    if name == "affec":
        return _load_affec(root, max_n, seed)
    # aliases for packs living under slightly different folder names
    aliases = {"coco": ["coco", "coco_outdoor_indoor", "coco_time_order", "coco_center_split"]}
    if name in aliases:
        for cand in aliases[name]:
            packed = _load_img_txt(root, cand, max_n, seed)
            if packed is not None:
                feats, y, mods, _ = packed
                return feats, y, mods, name
        return None
    return _load_img_txt(root, name, max_n, seed)


# ---------------------------------------------------------------------------
# sensors within a window
# ---------------------------------------------------------------------------


def window_sensors(
    xw: Mapping[str, np.ndarray],
    yw: np.ndarray,
    mods: Sequence[str],
    *,
    proto: ModalityPrototypeBank,
    prev_x: Optional[Mapping[str, np.ndarray]],
    prev_y: Optional[np.ndarray],
    seed: int,
) -> Dict[str, Any]:
    """PO-risk / MMD / proto / uni-acc sensors (numpy; external attributor)."""
    yw = np.asarray(yw).astype(np.int64).ravel()
    po: Dict[str, float] = {}
    mmd: Dict[str, float] = {}
    for i, m in enumerate(mods):
        Xm = np.asarray(xw[m], dtype=np.float64)
        if prev_x is not None and prev_y is not None:
            X0 = np.asarray(prev_x[m], dtype=np.float64)
            Y0 = np.asarray(prev_y, dtype=np.float64).ravel()
            Y1 = yw.astype(np.float64)
            po[m] = float(residual_concept(X0, Y0, Xm, Y1, seed=seed + 17 * i)) + 1e-6
            mmd[m] = float(rbf_mmd2(X0, Xm, max_n=min(64, len(X0), len(Xm)), seed=seed + 3 * i))
        else:
            # bootstrap: half-split residual within window
            mid = max(len(yw) // 2, 8)
            po[m] = float(
                residual_concept(
                    Xm[:mid],
                    yw[:mid].astype(np.float64),
                    Xm[mid:],
                    yw[mid:].astype(np.float64),
                    seed=seed + 17 * i,
                )
            ) + 1e-6
            pos = Xm[yw == 1]
            neg = Xm[yw == 0]
            if len(pos) > 2 and len(neg) > 2:
                mmd[m] = float(rbf_mmd2(pos, neg, max_n=min(48, len(pos), len(neg)), seed=seed + 5 * i))
            else:
                mmd[m] = 0.0

    proto.update(xw, yw)
    proto_d = proto.drift(xw, yw)

    # cheap uni-acc proxy via nearest class mean in feature space
    uni_acc: Dict[str, float] = {}
    for m in mods:
        Xm = np.asarray(xw[m], dtype=np.float64)
        acc = 0.5
        if (yw == 0).sum() >= 2 and (yw == 1).sum() >= 2:
            mu0 = Xm[yw == 0].mean(0)
            mu1 = Xm[yw == 1].mean(0)
            d0 = np.linalg.norm(Xm - mu0, axis=1)
            d1 = np.linalg.norm(Xm - mu1, axis=1)
            pred = (d1 < d0).astype(np.int64)
            acc = float((pred == yw).mean())
        uni_acc[m] = acc

    vimp = {m: float(uni_acc[m] * po[m]) for m in mods}
    return {
        "po": po,
        "mmd": mmd,
        "proto": proto_d,
        "uni_acc": uni_acc,
        "vimp": vimp,
    }


# ---------------------------------------------------------------------------
# one rollout
# ---------------------------------------------------------------------------


def _weighted_ce(logits: "torch.Tensor", yb: "torch.Tensor", w_row: "torch.Tensor | None"):
    """CE with optional per-row IPTW (R3)."""
    if w_row is None:
        return F.cross_entropy(logits, yb)
    per = F.cross_entropy(logits, yb, reduction="none")
    return (per * w_row).mean()


def run_version(
    feats: Dict[str, np.ndarray],
    y: np.ndarray,
    mods: Sequence[str],
    *,
    version: str,
    window: int,
    n_windows: int,
    seed: int,
    device: torch.device,
    metric_cfg: PORiskMetricConfig,
    act_cfg: NextStepActuatorConfig,
    step_mode: str = "per_mod",
    row_weight_mode: str = "sqrt",
    enable_row_iptw: bool = True,
) -> Dict[str, Any]:
    xs, ys = _pack_windows(feats, y, window=window, n_windows=n_windows, seed=seed)
    torch.manual_seed(seed + 101)
    np.random.seed(seed + 101)
    dims = {m: int(feats[m].shape[1]) for m in mods}
    model = StackFusion(dims, mods, fuse=64, n_class=2).to(device)
    model.train()
    proto = ModalityPrototypeBank(list(mods), ema=0.85)

    act = next_step_actuators({m: 1.0 / len(mods) for m in mods}, mods, cfg=act_cfg)
    prev_po = None
    po_ema = None
    alpha_hist: List[Dict[str, float]] = []
    prev_x = None
    prev_y = None
    # R3: reject at t → weights for t+1 (causal, same as modality actuators)
    next_w_np: Optional[np.ndarray] = None
    next_rejected = False
    prev_probe_err: Optional[float] = None

    accs, mses, flops, ents = [], [], [], []
    rows = []
    reject_mses = []
    calm_ws = []
    reject_sources: Dict[str, int] = {}

    for t, (xw, yw) in enumerate(zip(xs, ys)):
        xb = {m: torch.from_numpy(np.asarray(xw[m], dtype=np.float32)).to(device) for m in mods}
        yb = torch.from_numpy(np.asarray(yw, dtype=np.int64)).to(device)

        # row weights decided at end of previous window
        if enable_row_iptw and next_w_np is not None and next_rejected:
            w_row = torch.from_numpy(np.asarray(next_w_np, dtype=np.float32)).to(device)
            used_reject_w = True
        else:
            w_row = None
            used_reject_w = False
            if enable_row_iptw:
                calm_ws.append(1.0)

        # --- train THIS window with actuators decided at end of previous window ---
        base_lr = 1e-2
        for m in mods:
            frozen = bool(act["freeze_mask"].get(m, False))
            for p in list(model.projs[m].parameters()) + list(model.heads[m].parameters()):
                p.requires_grad_(not frozen)

        realized = realize_step_alloc(
            act["step_alloc"], act["freeze_mask"], mods, redistribute=False
        )
        steps_done = {m: 0 for m in mods}

        if step_mode == "shared":
            param_groups = []
            for m in mods:
                if act["freeze_mask"].get(m, False):
                    continue
                lr_m = base_lr * float(act["lr_mult"][m]) * float(act["lr_shared"])
                param_groups.append(
                    {
                        "params": list(model.projs[m].parameters())
                        + list(model.heads[m].parameters()),
                        "lr": max(lr_m, 1e-5),
                    }
                )
            param_groups.append(
                {"params": [model.stack_logits], "lr": base_lr * float(act["lr_shared"])}
            )
            if not any(not act["freeze_mask"].get(m, False) for m in mods):
                top = act["top_mod"]
                for p in list(model.projs[top].parameters()) + list(model.heads[top].parameters()):
                    p.requires_grad_(True)
                param_groups.insert(
                    0,
                    {
                        "params": list(model.projs[top].parameters())
                        + list(model.heads[top].parameters()),
                        "lr": base_lr * float(act["lr_shared"]),
                    },
                )
            opt = torch.optim.SGD(param_groups, momentum=0.0)
            n_steps = max(
                2, min(12, int(sum(act["step_alloc"].values()) // max(len(mods), 1)))
            )
            for _ in range(n_steps):
                logits, _h, _lm, w = model(xb, return_parts=True)
                ce = _weighted_ce(logits, yb, w_row)
                kl_pack = alpha_stack_kl(w, act["stack_prior"], mods, lambda_kl=0.10)
                loss = ce + kl_pack["loss"]
                opt.zero_grad(set_to_none=True)
                loss.backward()
                opt.step()
            flops_this = float(act["flops_rel"])
        else:
            schedule = expand_step_schedule(realized, mods, mode="block")
            if not schedule:
                top = act["top_mod"]
                for p in list(model.projs[top].parameters()) + list(model.heads[top].parameters()):
                    p.requires_grad_(True)
                schedule = [top]
                realized = {m: int(m == top) for m in mods}
            for m_upd in schedule:
                lr_m = base_lr * float(act["lr_mult"][m_upd]) * float(act["lr_shared"])
                groups = [
                    {
                        "params": list(model.projs[m_upd].parameters())
                        + list(model.heads[m_upd].parameters()),
                        "lr": max(lr_m, 1e-5),
                    },
                    {
                        "params": [model.stack_logits],
                        "lr": base_lr * float(act["lr_shared"]),
                    },
                ]
                opt_m = torch.optim.SGD(groups, momentum=0.0)
                logits, _h, _lm, w = model(xb, return_parts=True)
                ce = _weighted_ce(logits, yb, w_row)
                kl_pack = alpha_stack_kl(w, act["stack_prior"], mods, lambda_kl=0.10)
                loss = ce + kl_pack["loss"]
                opt_m.zero_grad(set_to_none=True)
                loss.backward()
                for m_o in mods:
                    if m_o == m_upd:
                        continue
                    for p in list(model.projs[m_o].parameters()) + list(
                        model.heads[m_o].parameters()
                    ):
                        if p.grad is not None:
                            p.grad = None
                opt_m.step()
                steps_done[m_upd] = steps_done.get(m_upd, 0) + 1
            flops_this = step_flops_rel(realized, total_steps=act_cfg.total_steps)

        with torch.no_grad():
            logits, _h, _lm, w = model(xb, return_parts=True)
            pred = logits.argmax(-1)
            acc = float((pred == yb).float().mean().item())
            mse = float(((pred.float() - yb.float()) ** 2).mean().item())
            proba = torch.softmax(logits, dim=-1)[:, 1].detach().cpu().numpy()

        if used_reject_w:
            reject_mses.append(mse)

        # --- sensors AFTER train → α / actuators + reject→weights for NEXT window ---
        sens = window_sensors(
            xw,
            yw,
            mods,
            proto=proto,
            prev_x=prev_x,
            prev_y=prev_y,
            seed=seed + 100 * t,
        )
        pack = metric_to_alpha(
            version,
            mods,
            po=sens["po"],
            mmd=sens["mmd"],
            proto=sens["proto"],
            uni_acc=sens["uni_acc"],
            vimp=sens["vimp"],
            po_prev=prev_po,
            po_ema=po_ema,
            alpha_hist=alpha_hist,
            cfg=metric_cfg,
        )
        alpha = pack["alpha"]
        if version == "po_fuse":
            next_act = next_step_actuators_fused(pack, mods, cfg=act_cfg)
        else:
            next_act = next_step_actuators(alpha, mods, cfg=act_cfg)
        alpha_hist.append(dict(alpha))
        prev_po_for_reject = prev_po
        prev_po = dict(sens["po"])
        if po_ema is None:
            po_ema = dict(sens["po"])
        else:
            for m in mods:
                po_ema[m] = 0.8 * float(po_ema[m]) + 0.2 * float(sens["po"][m])

        # R3 event-time row weights: RFPerm flag > hop OOS > proxy
        probe_err = float(1.0 - acc)  # shallow OOS proxy without sklearn RF
        rej = resolve_stream_reject(
            external_rejected=None,  # wire OnlineRFPerm flag here when available
            e_now=probe_err,
            e_prev=prev_probe_err,
            oos_gate=1.5,
            po_mods=sens["po"],
            mmd_mods=sens.get("mmd"),
            po_prev=prev_po_for_reject,
            mods=mods,
            use_proxy_fallback=True,
        )
        prev_probe_err = probe_err
        src = str(rej.get("source") or "none")
        reject_sources[src] = reject_sources.get(src, 0) + 1
        po_i = row_po_residual(np.asarray(yw), proba)
        if enable_row_iptw and rej["rejected"]:
            next_w_np = po_iptw_weights(
                po_i, mode=row_weight_mode, rejected=True
            )
            next_rejected = True
        else:
            next_w_np = np.ones(len(yw), dtype=float)
            next_rejected = False

        accs.append(acc)
        mses.append(mse)
        flops.append(flops_this)
        ents.append(float(act["alpha_entropy"]))
        row = {
            "t": t,
            "acc": acc,
            "mse": mse,
            "alpha": {m: float(alpha[m]) for m in mods},
            "lr_mult": {m: float(act["lr_mult"][m]) for m in mods},
            "step_alloc": {m: int(act["step_alloc"][m]) for m in mods},
            "step_realized": {m: int(realized.get(m, 0)) for m in mods},
            "steps_done": {m: int(steps_done.get(m, 0)) for m in mods},
            "step_mode": step_mode,
            "freeze": {m: bool(act["freeze_mask"][m]) for m in mods},
            "flops_rel": flops_this,
            "top_mod": act["top_mod"],
            "po": {m: float(sens["po"][m]) for m in mods},
            "rejected_for_next": bool(next_rejected),
            "reject_event": rej,
            "reject_proxy": rej,  # back-compat alias
            "row_weight_mode": row_weight_mode if next_rejected else "uniform",
            "used_reject_w": used_reject_w,
            "mean_row_w_next": float(np.mean(next_w_np)) if next_w_np is not None else 1.0,
            "probe_err": probe_err,
        }
        diag = dict(pack.get("diag") or {})
        if version == "po_fuse":
            row["fuse"] = {
                "omega_long": diag.get("omega_long"),
                "omega_short": diag.get("omega_short"),
                "spike_ratio": diag.get("spike_ratio"),
                "top_concept_mod": diag.get("top_concept_mod"),
                "top_spike_mod": diag.get("top_spike_mod"),
            }
            if "step_tilt" in next_act:
                row["step_tilt"] = {
                    m: float(next_act["step_tilt"].get(m, 0.0)) for m in mods
                }
        rows.append(row)
        act = next_act
        prev_x = {m: np.asarray(xw[m]) for m in mods}
        prev_y = np.asarray(yw)

    cg = continuous_gain_metrics(rows, mods=mods, acc_star=0.55)
    n_reject = sum(1 for r in rows if r.get("rejected_for_next"))
    return {
        "version": version,
        "n_windows": len(accs),
        "mean_acc_post": float(np.mean(accs)) if accs else 0.0,
        "mean_mse_post": float(np.mean(mses)) if mses else 0.0,
        "mean_flops_rel": float(np.mean(flops)) if flops else 1.0,
        "mean_alpha_entropy": float(np.mean(ents)) if ents else 0.0,
        "final_alpha": rows[-1]["alpha"] if rows else {},
        "rows": rows,
        "continuous": cg,
        "cum_flops": cg["cum_flops"],
        "t_to_acc_star": cg["t_to_acc_star"],
        "cum_flops_to_acc_star": cg["cum_flops_to_acc_star"],
        "mean_freeze_jaccard": cg["mean_freeze_jaccard"],
        "step_mode": step_mode,
        "row_iptw": {
            "enabled": enable_row_iptw,
            "mode": row_weight_mode,
            "n_reject_events": n_reject,
            "mean_mse_on_reject_fit": float(np.mean(reject_mses)) if reject_mses else None,
            "calm_w_ones": bool(calm_ws) and all(c == 1.0 for c in calm_ws),
            "reject_sources": reject_sources,
        },
    }


def summarize_cell(cell: Dict[str, Any], baseline: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "version": cell["version"],
        "mean_acc_post": cell["mean_acc_post"],
        "mean_mse_post": cell["mean_mse_post"],
        "mean_flops_rel": cell["mean_flops_rel"],
        "mean_alpha_entropy": cell["mean_alpha_entropy"],
        "mean_acc_lift": float(cell["mean_acc_post"] - baseline["mean_acc_post"]),
        "mean_mse_drop": float(baseline["mean_mse_post"] - cell["mean_mse_post"]),
        "final_alpha": cell["final_alpha"],
        "cum_flops": cell.get("cum_flops"),
        "t_to_acc_star": cell.get("t_to_acc_star"),
        "cum_flops_to_acc_star": cell.get("cum_flops_to_acc_star"),
        "mean_freeze_jaccard": cell.get("mean_freeze_jaccard"),
    }


# ---------------------------------------------------------------------------
# docs
# ---------------------------------------------------------------------------


def write_docs(out: Path, payload: Dict[str, Any]) -> None:
    lines = [
        "# AGOD PO-risk → next-step training metric compare",
        "",
        "Iterate metric versions that map **PO-risk sensors → α → next-window training actuators**",
        "(LR multipliers, step budget, freeze mask, stack prior).",
        "",
        "Causal loop: sensors/α at window `t` only set actuators for window `t+1`",
        "(no same-window leakage).",
        "",
        "## Metric versions",
        "",
        "| version | idea | next-step train role |",
        "|---|---|---|",
        "| `equal` | uniform α | dense baseline (all mods equal LR/steps) |",
        "| `po_soft` | Softmax(PO / τ) | steer LR/steps toward high residual-concept mods |",
        "| `po_minus_cov` | Softmax((PO − λ·MMD) / τ) | discount covariate mush before spending steps |",
        "| `po_gated` | drift-vs-noise gate → Softmax | freeze/skip noisy windows; FLOPs saver |",
        "| `po_proto` | Softmax(PO·(1+proto) − λ·MMD) | amplify true centroid moves |",
        "| `po_delta` | Softmax(EMA(PO) + γ·ΔPO) | anticipatory reallocation before Acc drops |",
        "| `po_budget` | floor + Softmax(PO) | keep all mods warm; soft reweight only |",
        "| `po_next` | α-hist ⊕ ΔPO forecast | next-α forecast for stack prior + LR |",
        "| `po_fuse` | long⊗short: freeze←L, steps←S, LR←α | concept-mod emphasis; spike-adaptive mix |",
        "",
        "Scenario gains summarize: [`AGOD_PO_Boost_PostTrain_Scenarios.tex`](AGOD_PO_Boost_PostTrain_Scenarios.tex).",
        "",
        "## Actuators (next window)",
        "",
        "- `lr_mult[m]` from simplex α (`alpha_to_lr`, β controls floor)",
        "- `step_alloc` proportional to α (sum ≈ `total_steps`)",
        "- `freeze_mask` when α_m < `freeze_theta` (BWD off; FWD still on)",
        "- `stack_prior = α` for KL(stack_w ‖ α)",
        "- efficiency proxy: `flops_rel = (#active mods) / M`",
        "",
    ]
    cross_votes: Dict[str, int] = {}
    for ds, block in payload["datasets"].items():
        lines += [
            f"## {ds}",
            "",
            f"mods = `{block.get('mods')}`",
            "",
            "| version | Acc↑ | MSE↓ | Acc lift vs equal | MSE drop vs equal | FLOPs_rel | H(α) | Jaccard | T(Acc*) | cumFLOPs@★ |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for v, c in block["summary"].items():
            t_star = c.get("t_to_acc_star")
            t_s = "—" if t_star is None else str(t_star)
            cum_star = c.get("cum_flops_to_acc_star")
            cum_s = "—" if cum_star is None else f"{float(cum_star):.2f}"
            jac = c.get("mean_freeze_jaccard")
            jac_s = "—" if jac is None else f"{float(jac):.3f}"
            lines.append(
                f"| `{v}` | {c['mean_acc_post']:.4f} | {c['mean_mse_post']:.4f} | "
                f"{c['mean_acc_lift']:+.4f} | {c['mean_mse_drop']:+.4f} | "
                f"{c['mean_flops_rel']:.3f} | {c['mean_alpha_entropy']:.3f} | "
                f"{jac_s} | {t_s} | {cum_s} |"
            )
        lines.append("")
        if block.get("ranking"):
            lines += ["### Opportunity ranking (vs `equal`)", ""]
            for i, r in enumerate(block["ranking"][:5], 1):
                lines.append(
                    f"{i}. `{r['version']}` opp_score={r['opp_score']:.4f} "
                    f"(acc_lift={r['mean_acc_lift']:+.4f}, mse_drop={r['mean_mse_drop']:+.4f}, "
                    f"flops={r['mean_flops_rel']:.3f})"
                )
                cross_votes[r["version"]] = cross_votes.get(r["version"], 0) + max(0, 6 - i)
            lines.append("")

    if cross_votes:
        ranked = sorted(cross_votes.items(), key=lambda kv: -kv[1])
        lines += [
            "## Cross-dataset opportunity votes",
            "",
            "| version | vote mass (top-5 ranks) |",
            "|---|---:|",
        ]
        for v, s in ranked:
            lines.append(f"| `{v}` | {s} |")
        lines.append("")

    lines += [
        "## Next-step training opportunities (actionable)",
        "",
        "### A. When PO-risk should change the train plan",
        "",
        "1. **Asymmetric multi-mod packs (Affec-like, M≥3)** — largest headroom.",
        "   Softmax(PO) peaking lets `freeze_mask` drop low-α towers → FLOPs↓ with Acc≈flat.",
        "   Best smoke winners: `po_gated` (FLOPs saver), `po_proto` (small Acc lift).",
        "2. **Balanced img/txt packs already easy (COCO ~0.90 Acc)** — LR reallocation is near-noop.",
        "   Keep `equal` or `po_budget` (warm floor); do not freeze; spend inventiveness on stack prior only.",
        "3. **Non-stationary streams** — prefer `po_delta` / `po_next` so next-window LR moves *before*",
        "   Acc collapses; pair with higher `gamma_delta` when ΔPO is reliable.",
        "",
        "### B. Actuator recipes for the next training loop",
        "",
        "| situation | metric | LR | steps | freeze | stack prior |",
        "|---|---|---|---|---|---|",
        "| concept spike on one mod | `po_soft` / `po_proto` | β≈0.2, boost top-α | dump steps to top-1/2 | freeze α<θ | KL→α |",
        "| high MMD, flat uni-acc | `po_minus_cov` / `po_gated` | damp shared LR | cut total_steps | freeze gated-off | weak KL |",
        "| need always-on towers | `po_budget` | soft only | equal-ish | never | KL→α |",
        "| forecast next risk | `po_next` | use forecast α | from forecast | optional | KL→forecast α |",
        "| already saturated Acc | `equal` | flat | flat | off | optional |",
        "",
        "### C. What not to do",
        "",
        "- Do **not** backprop into α (external attributor stays PO/FSDS/VIMP).",
        "- Do **not** treat KL(stack_w‖α) as MoE aux — prior is external, `w` is the only learned mixer.",
        "- Do **not** hard-gate FWD; freeze is BWD-only (efficiency = update FLOPs).",
        "- On 2-mod packs, freeze_theta must be high (~0.35+) to ever fire; soft LR is the real lever.",
        "",
        "### D. Suggested next training adjustments (priority)",
        "",
        "1. Affec online loop: switch default metric `equal` → `po_gated` (or `po_proto` if Acc-first).",
        "2. ~~Wire `step_alloc` into real optimizer step counts~~ **done (R2, `step_mode=per_mod`)**; ablate vs `--step-mode shared`.",
        "3. Add holdout Acc@FLOPs Pareto (freeze_theta sweep) before claiming efficiency wins.",
        "4. For Food-101 / Fashion-IQ: keep stacking KL prior = α, but keep freeze off; try `po_budget`.",
        "5. Log per-window (PO, α, lr_mult, step_realized, freeze, Acc) — use to tune τ / λ_cov / freeze_theta.",
        "",
        "### E. R2 step wiring",
        "",
        "- `realize_step_alloc`: freeze → 0 steps (FLOPs cut; no redistribute by default).",
        "- `expand_step_schedule(..., mode='block')`: highest budget modality dumped first.",
        "- `step_flops_rel = steps_used / total_steps`.",
        "- CLI: `--step-mode {per_mod,shared}` (default `per_mod`).",
        "",
        "### F. R3 reject-gated row IPTW (same stream)",
        "",
        "- `resolve_stream_reject`: external RFPerm flag → hop OOS (`e_now/e_prev`) → proxy.",
        "- Causal: reject at \(t\) weights Fit at \(t+1\); calm \(w=1\).",
        "- Orthogonal to modality freeze/steps; CLI `--row-weight-mode` / `--no-row-iptw`.",
        "",
    ]
    (out / "AGOD_po_risk_train_compare.md").write_text("\n".join(lines), encoding="utf-8")

    tex = [
        "% AGOD PO-risk → next-step train tables",
        "\\begin{table}[t]",
        "\\centering",
        "\\small",
        "\\caption{PO-risk metric versions $\\rightarrow$ next-window training actuators "
        "(causal: $\\alpha_t$ actuates window $t{+}1$).}",
        "\\label{tab:agod-po-risk-train}",
        "\\begin{tabular}{llrrrrr}",
        "\\toprule",
        "Dataset & Version & Acc & MSE & $\\Delta$Acc & FLOPs & $H(\\alpha)$ \\\\",
        "\\midrule",
    ]
    for ds, block in payload["datasets"].items():
        first = True
        for v, c in block["summary"].items():
            ds_cell = ds.replace("_", "\\_") if first else ""
            first = False
            vv = v.replace("_", "\\_")
            tex.append(
                f"{ds_cell} & \\texttt{{{vv}}} & {c['mean_acc_post']:.3f} & "
                f"{c['mean_mse_post']:.3f} & {c['mean_acc_lift']:+.3f} & "
                f"{c['mean_flops_rel']:.2f} & {c['mean_alpha_entropy']:.2f} \\\\"
            )
        tex.append("\\midrule")
    if tex[-1] == "\\midrule":
        tex[-1] = "\\bottomrule"
    tex += ["\\end{tabular}", "\\end{table}", ""]
    (out / "AGOD_po_risk_train_tables_only.tex").write_text("\n".join(tex), encoding="utf-8")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--datasets", nargs="+", default=["affec", "food101", "fashion_iq", "coco"])
    ap.add_argument("--max-n", type=int, default=640)
    ap.add_argument("--window", type=int, default=64)
    ap.add_argument("--n-windows", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--versions", nargs="+", default=list(METRIC_VERSIONS))
    ap.add_argument("--out", type=Path, default=Path("results/agod_po_risk_train"))
    ap.add_argument(
        "--step-mode",
        choices=("per_mod", "shared"),
        default="per_mod",
        help="R2: per_mod = real step_alloc budgets; shared = legacy loop-length scale",
    )
    ap.add_argument(
        "--row-weight-mode",
        choices=("sqrt", "cbrt", "prop", "uniform"),
        default="sqrt",
        help="R3: IPTW on reject windows only (calm stays w=1)",
    )
    ap.add_argument(
        "--no-row-iptw",
        action="store_true",
        help="Disable R3 reject-gated sample weights",
    )
    ap.add_argument(
        "--no-schedule-card",
        action="store_true",
        help="Disable R4 M>=3/M=2 schedule card defaults",
    )
    ap.add_argument(
        "--no-export-roi",
        action="store_true",
        help="Disable R5 ROI sqlite/jsonl export",
    )
    args = ap.parse_args()
    apply_card = not args.no_schedule_card
    do_export = not args.no_export_roi

    device = torch.device("cpu")
    metric_cfg = PORiskMetricConfig()
    act_cfg = NextStepActuatorConfig(beta_lr=0.20, freeze_theta=0.14, total_steps=40, min_steps=2)

    payload: Dict[str, Any] = {
        "datasets": {},
        "versions": list(args.versions),
        "config": {
            "max_n": args.max_n,
            "window": args.window,
            "n_windows": args.n_windows,
            "seed": args.seed,
            "step_mode": args.step_mode,
            "row_weight_mode": args.row_weight_mode,
            "enable_row_iptw": not args.no_row_iptw,
            "apply_schedule_card": apply_card,
            "export_roi": do_export,
        },
    }
    payload["config"]["metric_cfg"] = {
        "tau": metric_cfg.tau,
        "lam_cov": metric_cfg.lam_cov,
        "gamma_delta": metric_cfg.gamma_delta,
        "ema_po": metric_cfg.ema_po,
        "budget_floor": metric_cfg.budget_floor,
    }
    payload["config"]["act_cfg"] = {
        "beta_lr": act_cfg.beta_lr,
        "freeze_theta": act_cfg.freeze_theta,
        "total_steps": act_cfg.total_steps,
        "min_steps": act_cfg.min_steps,
    }

    schema_sql = args.root / "docs/agod/po_posttrain_roi_map.sql"
    roi_db = args.out / "po_posttrain_roi.sqlite"
    all_roi: List[Dict[str, Any]] = []

    for ds in args.datasets:
        packed = load_dataset(args.root, ds, args.max_n, args.seed)
        if packed is None:
            print(f"[skip] {ds}")
            continue
        feats, y, mods, name = packed
        print(f"== {name} mods={mods} n={y.shape[0]} ==")
        mcfg = PORiskMetricConfig(
            tau=metric_cfg.tau,
            lam_cov=metric_cfg.lam_cov,
            gamma_delta=metric_cfg.gamma_delta,
            ema_po=metric_cfg.ema_po,
            budget_floor=metric_cfg.budget_floor,
            omega_long=metric_cfg.omega_long,
            omega_short0=metric_cfg.omega_short0,
            spike_gain=metric_cfg.spike_gain,
        )
        acfg = NextStepActuatorConfig(
            beta_lr=act_cfg.beta_lr,
            freeze_theta=act_cfg.freeze_theta,
            total_steps=act_cfg.total_steps,
            min_steps=act_cfg.min_steps,
        )
        card = pick_default_schedule_card(len(mods)) if apply_card else None
        if card:
            mcfg.ema_po = float(card["ema_po"])
            mcfg.omega_long = float(card["omega_long"])
            mcfg.omega_short0 = float(card["omega_short0"])
            mcfg.spike_gain = float(card["spike_gain"])
            mcfg.tau = float(card["tau"])
            mcfg.budget_floor = float(card["budget_floor"])
            acfg.freeze_theta = float(card["freeze_theta"])
            print(f"  schedule_card={card['card_id']} freeze_theta={acfg.freeze_theta}")

        cells: Dict[str, Any] = {}
        for ver in args.versions:
            print(f"  [{ver}] ...", flush=True)
            cells[ver] = run_version(
                feats,
                y,
                mods,
                version=ver,
                window=args.window,
                n_windows=args.n_windows,
                seed=args.seed + 17,
                device=device,
                metric_cfg=mcfg,
                act_cfg=acfg,
                step_mode=args.step_mode,
                row_weight_mode=args.row_weight_mode,
                enable_row_iptw=not args.no_row_iptw,
            )
            print(
                f"    acc={cells[ver]['mean_acc_post']:.4f} "
                f"mse={cells[ver]['mean_mse_post']:.4f} "
                f"flops={cells[ver]['mean_flops_rel']:.3f}"
            )
            if do_export:
                recs = rows_to_roi_records(
                    cells[ver].get("rows") or [],
                    run_id=f"{name}:{ver}",
                    pack=name,
                    mods=mods,
                    acc_equal=None,
                    schedule_card_id=(card or {}).get("card_id"),
                )
                for rec in recs:
                    rec["version"] = ver
                    if recs and rec["window_t"] == recs[-1]["window_t"]:
                        rec["t_to_acc_star"] = cells[ver].get("t_to_acc_star")
                all_roi.extend(recs)

        if do_export and "equal" in cells:
            acc_eq = float(cells["equal"]["mean_acc_post"])
            for rec in all_roi:
                if rec.get("pack") == name:
                    rec["acc_equal"] = acc_eq

        base = cells["equal"]
        summary = {v: summarize_cell(cells[v], base) for v in cells}
        ranking = opportunity_rank(list(summary.values()), baseline="equal")
        light = {v: {k: cells[v][k] for k in cells[v] if k != "rows"} for v in cells}
        payload["datasets"][name] = {
            "mods": mods,
            "summary": summary,
            "ranking": ranking,
            "cells": light,
            "schedule_card": card,
        }
        top = ranking[0]["version"] if ranking else "n/a"
        print(f"  top opportunity: {top}")

    args.out.mkdir(parents=True, exist_ok=True)
    if do_export and all_roi and schema_sql.is_file():
        info = export_roi_sqlite(roi_db, all_roi, schema_sql=schema_sql)
        write_roi_jsonl(args.out / "po_roi_window_log.jsonl", all_roi)
        payload["roi_export"] = info
        print(f"ROI export: {info}")

    (args.out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    docs = args.root / "docs/agod"
    docs.mkdir(parents=True, exist_ok=True)
    write_docs(docs, payload)
    write_docs(args.out, payload)
    print(f"wrote {args.out}/summary.json")


if __name__ == "__main__":
    main()
