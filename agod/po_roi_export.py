"""Export PO-boost window decisions into ROI SQL schema (R5).

Maps compare ``rows`` → ``po_roi_window_log`` records; business reads
``docs/agod/po_posttrain_roi_map.sql`` views only.
"""
from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Mapping, Sequence


def rows_to_roi_records(
    rows: Sequence[Mapping[str, Any]],
    *,
    run_id: str,
    pack: str,
    mods: Sequence[str],
    acc_equal: float | None = None,
    schedule_card_id: str | None = None,
) -> list[dict]:
    """Convert compare trajectory rows into po_roi_window_log dicts."""
    mods = list(mods)
    n_mods = len(mods)
    out: list[dict] = []
    for r in rows:
        freeze = r.get("freeze") or {}
        active = [m for m in mods if not freeze.get(m, False)]
        if not active and mods:
            active = [r.get("top_mod") or mods[0]]
        fuse = r.get("fuse") or {}
        rej = r.get("reject_proxy") or {}
        rejected = int(
            bool(r.get("rejected_for_next") or rej.get("rejected") or r.get("used_reject_w"))
        )
        top_c = fuse.get("top_concept_mod") or r.get("top_mod")
        top_s = fuse.get("top_spike_mod") or r.get("top_mod")
        steps = r.get("steps_done") or r.get("step_realized") or {}
        steps_used = int(sum(int(v) for v in steps.values())) if steps else None
        rec = {
            "run_id": run_id,
            "window_t": int(r.get("t", 0)),
            "pack": pack,
            "top_concept_mod": top_c,
            "top_spike_mod": top_s,
            "active_mods_csv": ",".join(active),
            "n_active_mods": len(active),
            "n_mods": n_mods,
            "rejected": rejected,
            "n_rows": None,
            "n_hard_top20": None,
            "flops_rel": float(r.get("flops_rel", 1.0)),
            "wall_clock_s": None,
            "steps_used": steps_used,
            "t_to_acc_star": None,
            "acc": float(r.get("acc", 0.0)),
            "acc_equal": float(acc_equal) if acc_equal is not None else None,
            "next_mse": float(r["mse"]) if rejected and "mse" in r else None,
            "next_mse_uniform": None,
            "hard_p_at_20": None,
            "row_weight_mode": r.get("row_weight_mode"),
            "mean_row_w_next": r.get("mean_row_w_next"),
            "schedule_card_id": schedule_card_id,
            "version": r.get("version"),
        }
        out.append(rec)
    return out


def apply_roi_schema(conn: sqlite3.Connection, schema_sql: Path | str) -> None:
    path = Path(schema_sql)
    conn.executescript(path.read_text(encoding="utf-8"))
    # R5 additive columns (IF NOT EXISTS via try)
    for stmt in (
        "ALTER TABLE po_roi_window_log ADD COLUMN row_weight_mode TEXT",
        "ALTER TABLE po_roi_window_log ADD COLUMN mean_row_w_next REAL",
        "ALTER TABLE po_roi_window_log ADD COLUMN schedule_card_id TEXT",
        "ALTER TABLE po_roi_window_log ADD COLUMN version TEXT",
    ):
        try:
            conn.execute(stmt)
        except sqlite3.OperationalError:
            pass
    conn.commit()


def insert_roi_records(conn: sqlite3.Connection, records: Sequence[Mapping]) -> int:
    cols = [
        "run_id",
        "window_t",
        "pack",
        "top_concept_mod",
        "top_spike_mod",
        "active_mods_csv",
        "n_active_mods",
        "n_mods",
        "rejected",
        "n_rows",
        "n_hard_top20",
        "flops_rel",
        "wall_clock_s",
        "steps_used",
        "t_to_acc_star",
        "acc",
        "acc_equal",
        "next_mse",
        "next_mse_uniform",
        "hard_p_at_20",
        "row_weight_mode",
        "mean_row_w_next",
        "schedule_card_id",
        "version",
    ]
    sql = (
        f"INSERT OR REPLACE INTO po_roi_window_log ({', '.join(cols)}) "
        f"VALUES ({', '.join('?' for _ in cols)})"
    )
    n = 0
    for r in records:
        conn.execute(sql, [r.get(c) for c in cols])
        n += 1
    conn.commit()
    return n


def export_roi_sqlite(
    db_path: Path | str,
    records: Sequence[Mapping],
    *,
    schema_sql: Path | str,
) -> dict:
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    try:
        apply_roi_schema(conn, schema_sql)
        n = insert_roi_records(conn, records)
        # ship-gate snapshot
        cur = conn.execute(
            "SELECT COUNT(*) FROM v_roi_ship_gate"
            if _view_exists(conn, "v_roi_ship_gate")
            else "SELECT 0"
        )
        ship_rows = int(cur.fetchone()[0])
    finally:
        conn.close()
    return {"db": str(db_path), "n_inserted": n, "ship_gate_rows": ship_rows}


def _view_exists(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='view' AND name=?", (name,)
    ).fetchone()
    return row is not None


def write_roi_jsonl(path: Path | str, records: Sequence[Mapping]) -> int:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    return len(records)
