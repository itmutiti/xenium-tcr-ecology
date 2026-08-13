"""Unit tests for xenium_tcr_ecology.ingest.xenium_inventory (`02_raw_data_ingestion/04_inventory_xenium_files.py`)."""

from __future__ import annotations

import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.ingest.xenium_inventory import build_xenium_file_inventory

MANDATORY_SUFFIXES = [
    "_cell_boundaries.parquet.gz",
    "_cell_feature_matrix.h5",
    "_cells.parquet.gz",
    "_morphology.ome.tif.gz",
    "_nucleus_boundaries.parquet.gz",
    "_transcripts.parquet.gz",
]


def _make_complete_sample(staged_root, gsm: str, extra_files: list[str] | None = None):
    d = staged_root / gsm
    d.mkdir(parents=True)
    for suffix in MANDATORY_SUFFIXES:
        (d / f"{gsm}_x{suffix}").write_bytes(b"data")
    for extra in extra_files or []:
        (d / extra).write_bytes(b"data")


class TestBuildXeniumFileInventory:
    def test_passes_when_all_samples_complete(self, tmp_path):
        staged_root = tmp_path / "staged"
        _make_complete_sample(staged_root, "GSM1")
        _make_complete_sample(staged_root, "GSM2", extra_files=["GSM2_raw.tif.gz"])

        summary = build_xenium_file_inventory(
            staged_root, tmp_path / "inv.tsv", project_root=tmp_path
        )
        assert summary["samples_scanned"] == 2
        assert summary["samples_complete"] == 2

        report = (tmp_path / "inv.tsv").read_text()
        assert "raw_microscopy_optional" in report

    def test_raises_on_missing_mandatory_file(self, tmp_path):
        staged_root = tmp_path / "staged"
        d = staged_root / "GSM1"
        d.mkdir(parents=True)
        # Omit the transcripts file -- one of six mandatory roles.
        for suffix in MANDATORY_SUFFIXES:
            if "transcripts" in suffix:
                continue
            (d / f"GSM1_x{suffix}").write_bytes(b"data")

        with pytest.raises(PipelineError, match="missing mandatory files"):
            build_xenium_file_inventory(staged_root, tmp_path / "inv.tsv", project_root=tmp_path)

    def test_raises_on_missing_staged_dir(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            build_xenium_file_inventory(
                tmp_path / "nope", tmp_path / "inv.tsv", project_root=tmp_path
            )

    def test_raises_on_empty_staged_dir(self, tmp_path):
        (tmp_path / "staged").mkdir()
        with pytest.raises(PipelineError, match="No sample directories"):
            build_xenium_file_inventory(
                tmp_path / "staged", tmp_path / "inv.tsv", project_root=tmp_path
            )
