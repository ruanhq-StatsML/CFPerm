"""Ensure CFPerm Python/src is importable (PermuCATE helpers, shift boards)."""
from __future__ import annotations

import sys
from pathlib import Path

_PKG_ROOT = Path(__file__).resolve().parents[1]  # Python/po_risk_var
_REPO = _PKG_ROOT.parents[1]  # repo root
_SRC = _REPO / "Python" / "src"


def ensure_cfperm_src() -> Path:
    if _SRC.is_dir() and str(_SRC) not in sys.path:
        sys.path.insert(0, str(_SRC))
    return _SRC


REPO_ROOT = _REPO
PKG_ROOT = _PKG_ROOT
