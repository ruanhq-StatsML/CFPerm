#!/usr/bin/env python3
"""Run the disentangled style-transfer + FSDS/LOGO prototype."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "Python" / "src"
sys.path.insert(0, str(SRC))

from style_transfer_fsds import Config, train  # noqa: E402


if __name__ == "__main__":
    cfg = Config()
    train(cfg)
