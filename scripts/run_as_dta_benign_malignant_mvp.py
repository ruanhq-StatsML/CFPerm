#!/usr/bin/env python3
"""MVP: Benign vs Malignant OOD via PO magnitude + variance (+ anti-patterns).

  python3 scripts/run_as_dta_benign_malignant_mvp.py
  python3 scripts/run_as_dta_benign_malignant_mvp.py --quick
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_benign_malignant_mvp import run_mvp  # noqa: E402

OUT = ROOT / "results" / "as_dta_benign_malignant_mvp"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--max-pairs", type=int, default=900)
    ap.add_argument("--live-n", type=int, default=220)
    ap.add_argument("--K", type=int, default=36)
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    max_pairs = 500 if args.quick else int(args.max_pairs)
    live_n = 120 if args.quick else int(args.live_n)
    K = 20 if args.quick else int(args.K)

    payload = run_mvp(
        max_pairs=max_pairs, live_n=live_n, K=K, seed=int(args.seed), out_dir=OUT
    )
    for k in ("ref", "benign_ood", "malignant_ood"):
        b = payload["batches"][k]
        d = b.get("delta_vs_ref") or {}
        print(
            "[bm-mvp] %-14s mean=%.4f subStd=%.4f peak=%.4f cvar=%.4f  "
            "verdict=%s route=%s  d_mean=%s d_subStd=%s"
            % (
                k,
                b["diag"]["mean"],
                b["po_subenv"]["sub_std"],
                b["diag"]["peak_gap"],
                b["diag"]["cvar"],
                b["verdict"],
                b["route"],
                ("%.4f" % d["d_mean"]) if d else "-",
                ("%.4f" % d["d_sub_std"]) if d else "-",
            ),
            flush=True,
        )
    print("[bm-mvp] justify:", payload["justify"]["rule"], flush=True)
    print("[bm-mvp] anti:", *payload["justify"]["anti_patterns"], sep="\n  - ", flush=True)

    DOCS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "README.md", DOCS / "AS_DTA_benign_vs_malignant_ood.md")
    note = DOCS / "AS_DTA_benign_vs_malignant_ood.md"
    extra = """

## Scope

Resolve **benign vs malignant shift** (higher resolution than a single OOD score).
Headline axes: **Magnitude × Variance** — not Peak-first.

See `docs/method/AS_DTA_magnitude_variance_logic.md`.

## Hard filter & pseudo-label (why not)

| Anti-pattern | What it does | Failure mode |
| --- | --- | --- |
| Hard filter | Drop top-PO / high-OOD rows | Selection bias; deletes malignant *variance* signal |
| Pseudo-label | Flip Y on high PO | Label noise; benign shift looks "fixed"; fake uplift |

## Separability (merge)

1. Δμ quiet → `hold`
2. Δμ↑ & Δs flat → **benign shift** → `conformal_calibrate`
3. Δμ↑ & Δs↑ → **malignant shift** → `freeze_lexical_adv_aug`
4. Optional ops scalar: $J=\Delta\mu\cdot(1+\kappa\Delta s_{+})$

See also `results/as_dta_benign_malignant_mvp/`.
"""
    note.write_text(note.read_text() + extra)

    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "benign_vs_malignant_ood.png", ARTIFACTS / "as_dta_benign_vs_malignant_ood.png")
    print("[bm-mvp] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
