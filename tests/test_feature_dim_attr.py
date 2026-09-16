"""Tests for unified feature-dimension attribution (style + graph)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from scripts.agod.feature_dim_attr import (
    GRAPH_DIMS,
    LLM_DIM_ALIAS,
    STYLE_SUBFAMILIES,
    _synth_halu,
    _synth_hh,
    attribute_pack,
    build_audit_feature_dims,
    build_graph_feature_dims,
    build_style_as_feature_dims,
    dim_mass,
    main,
    pick_top_dims,
)

OUT = Path(__file__).resolve().parents[1] / "results" / "agod" / "feature_dim_attr"


def test_catalog_aliases_cover_graph_dims():
    assert set(LLM_DIM_ALIAS) == set(GRAPH_DIMS)
    assert set(STYLE_SUBFAMILIES) == {"length", "punct", "register"}


def test_pick_top_dims_falls_back_to_rf_mass():
    logo = {"a": 0.0, "b": 0.0, "c": 0.0}
    rf = {"a": 0.1, "b": 0.7, "c": 0.2}
    assert pick_top_dims(logo, rf, k=2) == ["b", "c"]


def test_dim_mass_sums_to_one():
    v = np.array([0.5, 0.3, 0.2])
    m = dim_mass(v, {"x": [0, 1], "y": [2]})
    assert abs(sum(m.values()) - 1.0) < 1e-9


def test_style_and_graph_packs_share_attribute_api():
    hh, hu = _synth_hh(240), _synth_halu(400)
    style = attribute_pack(build_style_as_feature_dims(hh), seed=0)
    graph = attribute_pack(build_graph_feature_dims(hu, seed=0, alias="llm"), seed=0, top_k=2)
    audit = attribute_pack(build_audit_feature_dims(hh, seed=0), seed=0)

    assert set(style["top_dims"]).issubset({"length", "punct", "register"})
    assert len(graph["top_dims"]) == 2
    for d in graph["top_dims"]:
        assert d in set(LLM_DIM_ALIAS.values())
    assert "within_dim" in graph and graph["top_dims"][0] in graph["within_dim"]
    assert set(audit["rf_mass"]) == {"preference_text", "style_register"}
    # style folded into graph author alias — not a separate product key
    assert "style_register" in graph["rf_mass"] or "author" in graph["rf_mass"]


def test_graph_alias_business_names():
    hu = _synth_halu(400)
    pack = build_graph_feature_dims(hu, seed=1, alias="graph")
    assert set(pack.dims) == set(GRAPH_DIMS)
    out = attribute_pack(pack, seed=1, top_k=2)
    assert set(out["top_dims"]).issubset(set(GRAPH_DIMS))


def test_main_synth_writes(tmp_path):
    assert main(["--synth", "--seed", "0"]) == 0
    summary = json.loads((OUT / "summary.json").read_text())
    assert summary["stance"].startswith("style + graph")
    assert "catalog" in summary
    assert (OUT / "REPORT.md").exists()
