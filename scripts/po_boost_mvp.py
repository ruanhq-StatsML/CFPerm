#!/usr/bin/env python3
"""PO-Boost MVP entry — one command for advisor demos.

Runs the no-Affec / no-torch synthetic continuous-gains smoke and prints
the MVP one-liner. Details: docs/agod/PO_Boost_MVP.md
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def _load_smoke():
    spec = importlib.util.spec_from_file_location(
        "smoke_po_boost_synthetic", _ROOT / "scripts/smoke_po_boost_synthetic.py"
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    smoke = _load_smoke()
    payload = smoke.run_synthetic_compare(seed=0)
    out = _ROOT / "results/agod_po_boost_synthetic"
    smoke.write_synthetic_docs(payload, out)

    print("=" * 60)
    print("PO-Boost MVP (synthetic M=5; no Affec / no torch)")
    print("  modality weights every window | row √PO only on reject")
    print("  spike → step dump | chronic L → freeze | ship: FLOPs<1 & ΔAcc≥-0.5%")
    print("=" * 60)
    fuse = payload["versions"].get("po_fuse", {})
    equal = payload["versions"].get("equal", {})
    print(
        f"equal:  FLOPs={equal.get('mean_flops_rel', 1):.3f}  "
        f"T*={equal.get('t_to_acc_star')}  ship={equal.get('ship_pass')}"
    )
    print(
        f"po_fuse: FLOPs={fuse.get('mean_flops_rel', 1):.3f}  "
        f"T*={fuse.get('t_to_acc_star')}  "
        f"ΔAcc={fuse.get('delta_acc', 0):+.3f}  ship={fuse.get('ship_pass')}"
    )
    print(f"reject_sources(fuse): {fuse.get('reject_sources')}")
    print(f"wrote {out}/continuous_gains_synthetic.md")
    print("MVP doc: docs/agod/PO_Boost_MVP.md")
    ok = bool(fuse.get("ship_pass"))
    print("MVP check:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
