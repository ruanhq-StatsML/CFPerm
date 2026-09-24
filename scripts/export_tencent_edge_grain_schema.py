#!/usr/bin/env python3
"""Export TencentGR feature_grid edge-grain column schema (for board docs).

  PYTHONPATH=. python3 scripts/export_tencent_edge_grain_schema.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scripts.run_sample_chunk_adjacent_board import TENCENT_CONVERT_LEAK, TENCENT_VOLUME

ROOT = Path(__file__).resolve().parents[1]
GRID = ROOT / "results/tencent_gr_tabular_ui/feature_grid.parquet"
OUT = ROOT / "results/sample_chunk_adjacent_board/edge_grain_schema.json"


def main() -> None:
    df = pd.read_parquet(GRID)
    keys = {"user_id", "item_id", "y_convert", "e_last_ts"}
    num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    ops = [c for c in num if c not in keys | TENCENT_CONVERT_LEAK]
    content = [c for c in ops if c not in TENCENT_VOLUME]
    sorted_df = df.sort_values("e_last_ts")
    spans = []
    for n in (1000, 2000):
        for i in range(min(5, len(sorted_df) // n)):
            sl = sorted_df.iloc[i * n : (i + 1) * n]
            spans.append(
                {
                    "N": n,
                    "chunk": i,
                    "n_edges": int(len(sl)),
                    "span_hours": float(
                        (sl["e_last_ts"].max() - sl["e_last_ts"].min()) / 3600.0
                    ),
                    "y_convert_rate": float(sl["y_convert"].mean()),
                }
            )
    payload = {
        "grain": "one_row_equals_one_user_item_edge",
        "table": str(GRID.relative_to(ROOT)),
        "n_rows": int(len(df)),
        "n_users": int(df["user_id"].nunique()),
        "n_items": int(df["item_id"].nunique()),
        "dup_user_item": int(df.duplicated(["user_id", "item_id"]).sum()),
        "sort_key": "e_last_ts",
        "label": "y_convert",
        "windowing": "sort_by_e_last_ts then chunk_id = row_index // N",
        "columns": {
            "keys_time_label": sorted(keys),
            "edge_e": sorted(c for c in df.columns if c.startswith("e_")),
            "user_u": sorted(c for c in df.columns if c.startswith("u_")),
            "item_i": sorted(
                c for c in df.columns if c.startswith("i_") or c.startswith("ui_")
            ),
            "ops_X": ops,
            "content_X": content,
            "convert_leak_dropped": sorted(TENCENT_CONVERT_LEAK),
            "volume_dropped_in_content": sorted(TENCENT_VOLUME),
        },
        "equal_count_calendar_span_smoke": spans,
        "note": (
            "Equal edge-count windows have unequal calendar width; "
            "that is intended for ad-log burstiness."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    md = [
        "# Edge grain schema (auto)",
        "",
        f"- rows={payload['n_rows']} · users={payload['n_users']} · items={payload['n_items']}",
        f"- grain: `{payload['grain']}`",
        f"- sort: `{payload['sort_key']}` · Y: `{payload['label']}`",
        f"- ops |X|={len(ops)} · content |X|={len(content)}",
        "",
        "## Equal-count calendar spans (smoke)",
        "",
        "| N | chunk | span_h | y_rate |",
        "|---:|---:|---:|---:|",
    ]
    for s in spans:
        md.append(
            f"| {s['N']} | {s['chunk']} | {s['span_hours']:.1f} | {s['y_convert_rate']:.4f} |"
        )
    md += [
        "",
        "## ops_X",
        "",
        ", ".join(f"`{c}`" for c in ops),
        "",
        "## content_X",
        "",
        ", ".join(f"`{c}`" for c in content),
        "",
    ]
    (OUT.parent / "EDGE_GRAIN_SCHEMA.md").write_text("\n".join(md) + "\n")
    print(json.dumps({"out": str(OUT), "n_ops": len(ops), "n_content": len(content)}, indent=2))


if __name__ == "__main__":
    main()
