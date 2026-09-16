#!/usr/bin/env python3
"""Print the serving-gate feature catalog and lock it against on-disk csvs.

Usage::

    PYTHONPATH=. python3 scripts/list_serving_features.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.serving_features import BLOCKS, disk_x_cols  # noqa: E402

OUT = ROOT / "results" / "manuscript" / "serving_gates"


def render_md() -> str:
    lines = [
        "# Serving-gate features (on disk)",
        "",
        "These `x_*` columns are the prototype. No prompts, no wiki text, no gold flags, no HH chosen.",
        "Gates: `PYTHONPATH=. python3 scripts/run_serving_gates.py`.",
        "",
    ]
    for b in BLOCKS:
        n = len(b["features"])
        gate = "live gate" if b["gate"] else "table shape only"
        lines += [
            f"## {b['facet']}",
            "",
            f"`{b['table']}` · n={b['n_rows']} · {n} X · {gate}",
            "",
            f"Y = {b['y']}. {b['note']}",
            "",
            "| column | meaning |",
            "|---|---|",
        ]
        for col, meaning in b["features"]:
            lines.append(f"| `{col}` | {meaning} |")
        lines.append("")
    lines += [
        "Not features: raw prompt/question, wiki paragraph text, gold/supporting flags, HH chosen/rejected, single-channel Recall, hallucination rate.",
        "",
    ]
    return "\n".join(lines)


def check_disk() -> list[str]:
    errors = []
    for b in BLOCKS:
        path = ROOT / b["table"]
        if not path.exists():
            errors.append(f"missing {b['table']}")
            continue
        got = disk_x_cols(path)
        want = [c for c, _ in b["features"]]
        if got != want:
            errors.append(f"{b['table']}: disk {got} != catalog {want}")
    return errors


def main() -> int:
    errors = check_disk()
    md = render_md()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "FEATURES.md").write_text(md, encoding="utf-8")
    payload = [
        {
            "id": b["id"],
            "facet": b["facet"],
            "table": b["table"],
            "n_x": len(b["features"]),
            "x": [c for c, _ in b["features"]],
            "gate": b["gate"],
        }
        for b in BLOCKS
    ]
    (OUT / "features.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(md)
    if errors:
        print("schema mismatches:", file=sys.stderr)
        for e in errors:
            print(" ", e, file=sys.stderr)
        return 1
    print("catalog matches disk.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
