"""R5 ROI export tests (no torch — bootstrap agod.po_roi_export only)."""
from __future__ import annotations

import importlib.util
import sqlite3
import sys
import types
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def _load_export():
    if "agod" not in sys.modules or not getattr(sys.modules["agod"], "__path__", None):
        pkg = types.ModuleType("agod")
        pkg.__path__ = [str(_ROOT / "agod")]
        sys.modules["agod"] = pkg
    name = "agod.po_roi_export"
    if name in sys.modules and hasattr(sys.modules[name], "rows_to_roi_records"):
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _ROOT / "agod/po_roi_export.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


_exp = _load_export()


def test_rows_to_roi_and_sqlite(tmp_path: Path):
    mods = ["img", "txt", "aud"]
    rows = [
        {
            "t": 0,
            "acc": 0.5,
            "mse": 0.4,
            "flops_rel": 0.7,
            "freeze": {"img": False, "txt": True, "aud": True},
            "top_mod": "img",
            "fuse": {"top_concept_mod": "img", "top_spike_mod": "txt"},
            "rejected_for_next": True,
            "row_weight_mode": "sqrt",
            "mean_row_w_next": 1.0,
            "steps_done": {"img": 20, "txt": 0, "aud": 0},
        },
        {
            "t": 1,
            "acc": 0.6,
            "mse": 0.3,
            "flops_rel": 0.7,
            "freeze": {"img": False, "txt": True, "aud": True},
            "top_mod": "img",
            "rejected_for_next": False,
            "row_weight_mode": "uniform",
            "steps_done": {"img": 18, "txt": 0, "aud": 0},
        },
    ]
    recs = _exp.rows_to_roi_records(
        rows,
        run_id="affec:po_fuse",
        pack="affec",
        mods=mods,
        acc_equal=0.5,
        schedule_card_id="M_ge_3",
    )
    assert len(recs) == 2
    assert recs[0]["n_active_mods"] == 1
    assert recs[0]["rejected"] == 1
    assert recs[0]["top_concept_mod"] == "img"

    schema = _ROOT / "docs/agod/po_posttrain_roi_map.sql"
    db = tmp_path / "roi.sqlite"
    info = _exp.export_roi_sqlite(db, recs, schema_sql=schema)
    assert info["n_inserted"] == 2
    conn = sqlite3.connect(str(db))
    n = conn.execute("SELECT COUNT(*) FROM po_roi_window_log").fetchone()[0]
    assert n == 2
    gate = conn.execute(
        "SELECT pass_modality_emphasis FROM v_roi_ship_gate WHERE pack='affec'"
    ).fetchone()
    assert gate is not None
    assert int(gate[0]) == 1
    conn.close()

    jl = tmp_path / "log.jsonl"
    assert _exp.write_roi_jsonl(jl, recs) == 2
    assert jl.read_text(encoding="utf-8").count("\n") == 2
