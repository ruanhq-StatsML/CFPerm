#!/usr/bin/env python3
"""PermuCATE RF Mag×Var + multimodal attribution MVP.

  python3 scripts/run_as_dta_permucate_magvar.py
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python" / "src"))

from as_dta_permucate_magvar import run_mvp  # noqa: E402

OUT = ROOT / "results" / "as_dta_permucate_magvar"
DOCS = ROOT / "docs" / "method"
ARTIFACTS = Path("/opt/cursor/artifacts")


def main():
    payload = run_mvp(out_dir=OUT)
    mag = payload["magvar_eval"]
    c = mag["contingency_2x2"]
    print("[pc-mvp] estimator:", mag["estimator"]["name"], mag["estimator"]["model_m"], flush=True)
    print(
        "[pc-mvp] used_permucate_rf=",
        mag["estimator"]["used_permucate_rf"],
        " unit=",
        mag["estimator"]["unit_po"],
        flush=True,
    )
    print(
        "[pc-mvp] 2x2 [[%d,%d],[%d,%d]] acc=%.2f"
        % (c["tn"], c["fp"], c["fn"], c["tp"], mag["accuracy"]),
        flush=True,
    )
    for k, row in mag["classified"].items():
        print(
            "[pc-mvp] %-10s pred=%s  d_mu=%.4f d_s=%.4f J=%.3f"
            % (k, row["pred"], row["d_mu"], row["d_s"], row["J"]),
            flush=True,
        )
    mm = payload["multimodal"]
    print("[pc-mvp] multimodal shares permute:", mm["modality_shares_permute"], flush=True)
    print("[pc-mvp] multimodal shares vimp:", mm["modality_shares_vimp"], flush=True)

    DOCS.mkdir(parents=True, exist_ok=True)
    # keep professor doc stable; only refresh results README pointer file
    (DOCS / "AS_DTA_permucate_magvar_results.md").write_text(
        (OUT / "README.md").read_text()
    )
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    shutil.copy(OUT / "permucate_magvar_mm.png", ARTIFACTS / "as_dta_permucate_magvar_mm.png")
    print("[pc-mvp] wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
