"""Unit tests for xenium_tcr_ecology.annotation.marker_registry (`06_cell_type_annotation/00_compile_marker_and_reference_registry.py`)."""

from __future__ import annotations

import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.annotation.marker_registry import (
    CELL_TYPE_MARKER_REGISTRY,
    build_marker_registry_table,
    validate_registry,
)


class TestValidateRegistry:
    def test_restricts_markers_to_available_genes(self):
        available = {"CD3D", "CD3E"}
        registry = [
            {
                "cell_identity": "T_cell",
                "hierarchy_level": "major_lineage",
                "parent_identity": None,
                "markers": ["CD3D", "CD3E", "NOT_PRESENT"],
                "confidence_tier": "high",
                "rationale": "x",
            }
        ]
        result = validate_registry(available, registry)
        assert result[0]["markers"] == ["CD3D", "CD3E"]
        assert result[0]["n_markers_in_panel"] == 2

    def test_raises_if_identity_loses_all_markers(self):
        available = {"UNRELATED_GENE"}
        registry = [
            {
                "cell_identity": "T_cell",
                "hierarchy_level": "major_lineage",
                "parent_identity": None,
                "markers": ["CD3D"],
                "confidence_tier": "high",
                "rationale": "x",
            }
        ]
        with pytest.raises(PipelineError, match="0 markers"):
            validate_registry(available, registry)

    def test_every_curated_identity_has_at_least_one_marker_defined(self):
        for entry in CELL_TYPE_MARKER_REGISTRY:
            assert len(entry["markers"]) >= 1, entry["cell_identity"]

    def test_every_substate_has_a_parent_that_exists(self):
        identities = {e["cell_identity"] for e in CELL_TYPE_MARKER_REGISTRY}
        for entry in CELL_TYPE_MARKER_REGISTRY:
            if entry["hierarchy_level"] == "substate":
                assert entry["parent_identity"] in identities, entry["cell_identity"]


class TestBuildMarkerRegistryTable:
    def test_produces_expected_columns(self):
        available = {g for entry in CELL_TYPE_MARKER_REGISTRY for g in entry["markers"]}
        table = build_marker_registry_table(available)
        for col in [
            "cell_identity",
            "hierarchy_level",
            "parent_identity",
            "markers",
            "n_markers_in_panel",
            "confidence_tier",
            "rationale",
            "registry_version",
        ]:
            assert col in table.columns
        assert len(table) == len(CELL_TYPE_MARKER_REGISTRY)

    def test_no_duplicate_cell_identities(self):
        available = {g for entry in CELL_TYPE_MARKER_REGISTRY for g in entry["markers"]}
        table = build_marker_registry_table(available)
        assert table["cell_identity"].is_unique
