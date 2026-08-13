"""Unit tests for xenium_tcr_ecology.clone_ecology.taxonomy_freeze (`11_clone_spatial_descriptors/07_freeze_provisional_taxonomy_version.py`)."""

from __future__ import annotations

import pandas as pd

from xenium_tcr_ecology.clone_ecology.taxonomy_freeze import (
    TAXONOMY_VERSION,
    build_taxonomy_version_entry,
    merge_taxonomy_log,
)


class TestBuildTaxonomyVersionEntry:
    def test_entry_has_expected_fields(self):
        entry = build_taxonomy_version_entry(
            structure_type="continuous",
            n_units=261,
            dominant_feature="cycling_fraction",
            dominant_loading=-0.998,
            frozen_date="2026-07-11",
        )
        assert entry["taxonomy_version"] == TAXONOMY_VERSION
        assert entry["structure_type"] == "continuous"
        assert entry["n_clone_sections"] == 261
        assert entry["dominant_feature"] == "cycling_fraction"
        assert entry["status"] == "pending_phase12_external_validation"


class TestMergeTaxonomyLog:
    def test_no_existing_log_creates_single_row(self):
        entry = build_taxonomy_version_entry(
            "continuous", 261, "cycling_fraction", -0.998, "2026-07-11"
        )
        result = merge_taxonomy_log(None, entry)
        assert len(result) == 1
        assert result.iloc[0]["taxonomy_version"] == TAXONOMY_VERSION

    def test_existing_different_version_is_preserved(self):
        existing = pd.DataFrame(
            [{"taxonomy_version": "v0_pilot", "structure_type": "discrete", "n_clone_sections": 10}]
        )
        entry = build_taxonomy_version_entry(
            "continuous", 261, "cycling_fraction", -0.998, "2026-07-11"
        )
        result = merge_taxonomy_log(existing, entry)
        assert len(result) == 2
        assert "v0_pilot" in result["taxonomy_version"].tolist()
        assert TAXONOMY_VERSION in result["taxonomy_version"].tolist()

    def test_existing_same_version_is_replaced_not_duplicated(self):
        old_entry = build_taxonomy_version_entry(
            "continuous", 200, "old_feature", 0.5, "2026-01-01"
        )
        existing = pd.DataFrame([old_entry])
        new_entry = build_taxonomy_version_entry(
            "continuous", 261, "cycling_fraction", -0.998, "2026-07-11"
        )
        result = merge_taxonomy_log(existing, new_entry)
        assert len(result) == 1
        assert result.iloc[0]["n_clone_sections"] == 261
        assert result.iloc[0]["dominant_feature"] == "cycling_fraction"
