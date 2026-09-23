"""Synthetic PO-boost stream smoke (no Affec / no torch).

Builds a multi-mod PO trajectory with a mid-stream spike, runs
metric→α→actuators + resolve_stream_reject, reports continuous gains
and reject sources. Fills R1 stub table until Affec cache exists.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path
from typing import Any

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]


def _boot():
    if "agod" not in sys.modules or not getattr(sys.modules["agod"], "__path__", None):
        pkg = types.ModuleType("agod")
        pkg.__path__ = [str(_ROOT / "agod")]
        sys.modules["agod"] = pkg

    def load(name: str, rel: str):
        if name in sys.modules and hasattr(sys.modules[name], "__file__"):
            return sys.modules[name]
        spec = importlib.util.spec_from_file_location(name, _ROOT / rel)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        assert spec.loader is not None
        spec.loader.exec_module(mod)
        return mod

    load("agod.lr_controller", "agod/lr_controller.py")
    load("agod.next_step", "agod/next_step.py")
    load("agod.smooth_router", "agod/smooth_router.py")
    pr = load("agod.po_risk_train", "agod/po_risk_train.py")
    rej = load("agod.reject_event", "agod/reject_event.py")
    return pr, rej


_pr, _rej = _boot()


def synth_po_stream(
    mods: list[str],
    n_windows: int = 12,
    *,
    spike_t: int = 4,
    spike_mod: str | None = None,
    seed: int = 0,
) -> list[dict[str, float]]:
    """Chronic imbalance + one spike window on ``spike_mod``."""
    rng = np.random.default_rng(seed)
    spike_mod = spike_mod or mods[0]
    # chronic: mods[0] elevated, mods[-1] near-floor (freeze candidate)
    stream = []
    for t in range(n_windows):
        po = {}
        for i, m in enumerate(mods):
            base = 0.08 - 0.012 * i  # decreasing chronic
            po[m] = float(max(0.01, base + 0.005 * rng.random()))
        if t == spike_t:
            po[spike_mod] = float(0.35 + 0.05 * rng.random())
        elif t > spike_t:
            po[spike_mod] = float(max(po[spike_mod], 0.12 + 0.02 * rng.random()))
        stream.append(po)
    return stream


def run_synthetic_compare(
    *,
    versions: list[str] | None = None,
    n_windows: int = 12,
    seed: int = 0,
) -> dict[str, Any]:
    mods = ["m0", "m1", "m2", "m3", "m4"]
    card = _pr.pick_default_schedule_card(len(mods))
    metric_cfg = _pr.PORiskMetricConfig(
        ema_po=float(card["ema_po"]),
        omega_long=float(card["omega_long"]),
        omega_short0=float(card["omega_short0"]),
        spike_gain=float(card["spike_gain"]),
        tau=float(card["tau"]),
        budget_floor=float(card["budget_floor"]),
    )
    act_cfg = _pr.NextStepActuatorConfig(
        freeze_theta=0.12,  # smoke: allow freeze under chronic imbalance
        total_steps=40,
        min_steps=2,
        beta_lr=0.20,
    )
    versions = versions or ["equal", "po_soft", "po_gated", "po_fuse"]
    po_stream = synth_po_stream(mods, n_windows, seed=seed)
    out: dict[str, Any] = {"mods": mods, "schedule_card": card, "versions": {}}

    for ver in versions:
        rows = []
        prev_po = None
        po_ema = None
        alpha_hist: list[dict[str, float]] = []
        act = _pr.next_step_actuators({m: 1.0 / len(mods) for m in mods}, mods, cfg=act_cfg)
        prev_probe = None
        reject_sources: dict[str, int] = {}
        acc = 0.40
        for t, po in enumerate(po_stream):
            realized = _pr.realize_step_alloc(
                act["step_alloc"], act["freeze_mask"], mods, redistribute=False
            )
            active_n = sum(1 for m in mods if not act["freeze_mask"].get(m, False))
            # FLOPs = BWD-active fraction (honest freeze accounting)
            flops = float(active_n / max(len(mods), 1))
            top = max(mods, key=lambda m: po[m])
            focus = float(realized.get(top, 0)) / max(act_cfg.total_steps, 1)
            # Acc climbs with focus; equal is slow uniform climb
            if ver == "equal":
                acc = float(np.clip(0.42 + 0.018 * t, 0, 0.95))
            else:
                acc = float(
                    np.clip(
                        acc + 0.055 * max(focus, 0.12) + 0.03 * (1.0 - flops),
                        0,
                        0.95,
                    )
                )

            # inject OOS hop at spike window so hop_oos can fire once
            probe_err = float(1.0 - acc)
            if t == 4:
                probe_err = float((prev_probe or 0.3) * 2.2)
            rej = _rej.resolve_stream_reject(
                e_now=probe_err,
                e_prev=prev_probe,
                oos_gate=1.5,
                po_mods=po,
                po_prev=prev_po,
                mods=mods,
                use_proxy_fallback=(t != 4),  # force hop path visibility at spike
            )
            prev_probe = probe_err
            src = str(rej.get("source") or "none")
            reject_sources[src] = reject_sources.get(src, 0) + 1
            po_i = np.array([po[m] for m in mods], float)
            _w = _pr.po_iptw_weights(po_i, mode="sqrt", rejected=bool(rej["rejected"]))

            pack = _pr.metric_to_alpha(
                ver,
                mods,
                po=po,
                mmd={m: 0.02 for m in mods},
                proto={m: 0.1 if m == top else 0.0 for m in mods},
                uni_acc={m: 0.55 for m in mods},
                vimp={m: 0.2 for m in mods},
                po_prev=prev_po,
                po_ema=po_ema,
                alpha_hist=alpha_hist,
                cfg=metric_cfg,
            )
            alpha = pack["alpha"]
            if ver == "po_fuse":
                next_act = _pr.next_step_actuators_fused(pack, mods, cfg=act_cfg)
            else:
                next_act = _pr.next_step_actuators(alpha, mods, cfg=act_cfg)

            rows.append(
                {
                    "t": t,
                    "acc": acc,
                    "flops_rel": flops,
                    "freeze": dict(act["freeze_mask"]),
                    "alpha": dict(alpha),
                    "rejected_for_next": bool(rej["rejected"]),
                    "reject_event": rej,
                    "mean_row_w": float(np.mean(_w)),
                }
            )
            alpha_hist.append(dict(alpha))
            prev_po = dict(po)
            if po_ema is None:
                po_ema = dict(po)
            else:
                po_ema = {
                    m: 0.8 * float(po_ema[m]) + 0.2 * float(po[m]) for m in mods
                }
            act = next_act

        cg = _pr.continuous_gain_metrics(rows, mods=mods, acc_star=0.55)
        out["versions"][ver] = {
            "continuous": cg,
            "mean_flops_rel": cg["mean_flops_rel"],
            "t_to_acc_star": cg["t_to_acc_star"],
            "cum_flops_to_acc_star": cg["cum_flops_to_acc_star"],
            "mean_freeze_jaccard": cg["mean_freeze_jaccard"],
            "final_acc": rows[-1]["acc"] if rows else 0.0,
            "n_reject": sum(1 for r in rows if r.get("rejected_for_next")),
            "reject_sources": reject_sources,
            "rows": rows,
        }
    # lifts vs equal
    base = out["versions"].get("equal", {})
    for ver, cell in out["versions"].items():
        cell["delta_acc"] = float(cell["final_acc"] - base.get("final_acc", cell["final_acc"]))
        cell["ship_pass"] = bool(
            cell["mean_flops_rel"] < 1.0 and cell["delta_acc"] >= -0.005
        )
    return out


def write_synthetic_docs(payload: dict[str, Any], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    lines = [
        "# PO-boost synthetic continuous gains (no Affec)",
        "",
        f"schedule_card=`{payload['schedule_card'].get('card_id')}` mods=`{payload['mods']}`",
        "",
        "| version | FLOPs_rel | T(Acc*) | cumFLOPs@★ | Jaccard | final Acc | ΔAcc | rejects | ship |",
        "|---|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for ver, c in payload["versions"].items():
        t = c.get("t_to_acc_star")
        cum = c.get("cum_flops_to_acc_star")
        lines.append(
            f"| `{ver}` | {c['mean_flops_rel']:.3f} | "
            f"{'—' if t is None else t} | "
            f"{'—' if cum is None else f'{cum:.2f}'} | "
            f"{c['mean_freeze_jaccard']:.3f} | {c['final_acc']:.3f} | "
            f"{c['delta_acc']:+.3f} | {c['n_reject']} | "
            f"{'Y' if c['ship_pass'] else 'N'} |"
        )
    lines += ["", "## Reject sources", ""]
    for ver, c in payload["versions"].items():
        lines.append(f"- `{ver}`: {c.get('reject_sources')}")
    (out / "continuous_gains_synthetic.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    light = {
        "mods": payload["mods"],
        "schedule_card": payload["schedule_card"],
        "versions": {
            v: {k: c[k] for k in c if k != "rows"} for v, c in payload["versions"].items()
        },
    }
    (out / "summary_synthetic.json").write_text(json.dumps(light, indent=2), encoding="utf-8")


def main() -> None:
    payload = run_synthetic_compare()
    out = Path("results/agod_po_boost_synthetic")
    write_synthetic_docs(payload, out)
    print(f"wrote {out}/continuous_gains_synthetic.md")
    for ver, c in payload["versions"].items():
        print(
            f"  {ver}: flops={c['mean_flops_rel']:.3f} "
            f"T*={c['t_to_acc_star']} ship={c['ship_pass']}"
        )


if __name__ == "__main__":
    main()
