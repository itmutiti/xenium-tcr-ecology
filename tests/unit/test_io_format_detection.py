"""Unit tests for xenium_tcr_ecology.io.format_detection (`03_spatialdata_import/00_detect_xenium_format_version.py`).

The OME-XML extraction and full report-building functions need a real (or
realistically-shaped) gzip-compressed OME-TIFF and are exercised against
the real project data in test_real_project_format_report below; the
standard-reader compatibility check is pure logic and tested directly.
"""

from __future__ import annotations

from xenium_tcr_ecology.io.format_detection import check_standard_reader_compatibility


class TestCheckStandardReaderCompatibility:
    def test_incompatible_when_required_files_absent(self, tmp_path):
        section_dir = tmp_path / "P01_run1"
        section_dir.mkdir()
        for f in ["cell_boundaries.parquet.gz", "cell_feature_matrix.h5"]:
            (section_dir / f).write_bytes(b"")

        result = check_standard_reader_compatibility(section_dir)
        assert result["compatible"] is False
        assert "experiment.xenium" in result["missing"]
        assert "cells.zarr.zip" in result["missing"]

    def test_compatible_when_required_files_present(self, tmp_path):
        section_dir = tmp_path / "P01_run1"
        section_dir.mkdir()
        for f in ["experiment.xenium", "cells.zarr.zip"]:
            (section_dir / f).write_bytes(b"")

        result = check_standard_reader_compatibility(section_dir)
        assert result["compatible"] is True
        assert result["missing"] == []


class TestRealProjectFormatReport:
    def test_real_standardised_sections_confirm_custom_reader_needed(self):
        """Exercises the real data.standardised/ sections produced during
        this project's Raw Data Ingestion run. Confirms, against the actual GEO
        release, that the standard spatialdata_io reader is unusable here,
        and that pixel size is consistent across sections."""
        import pytest

        from xenium_tcr_ecology.infra.paths import find_project_root
        from xenium_tcr_ecology.io.format_detection import build_format_report

        project_root = find_project_root()
        standardised_root = project_root / "data" / "standardised"
        if not standardised_root.is_dir() or not any(standardised_root.iterdir()):
            pytest.skip(
                "data/standardised/ not populated yet (Raw Data Ingestion not run in this environment)"
            )

        summary = build_format_report(
            standardised_root, project_root / "metadata" / "_test_format_versions.tsv", project_root
        )
        (project_root / "metadata" / "_test_format_versions.tsv").unlink()

        assert summary["sections_checked"] == 18
        assert summary["standard_reader_incompatible"] == 18  # every section
        assert (
            summary["distinct_pixel_size_um_combos"] == 1
        )  # one consistent scale across the cohort
