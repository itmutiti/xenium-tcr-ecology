"""Unit tests for xenium_tcr_ecology.validation.spatial_dataset_acquisition (`16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py`)."""

from __future__ import annotations

import hashlib

import pandas as pd
import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.validation.spatial_dataset_acquisition import (
    build_second_spatial_dataset_acquisition_summary,
    verify_checksums,
)


def _write_file_and_checksum(tmp_path, name: str, content: bytes) -> None:
    (tmp_path / name).write_bytes(content)
    digest = hashlib.sha256(content).hexdigest()
    checksums_path = tmp_path / "checksums.sha256"
    existing = checksums_path.read_text() if checksums_path.exists() else ""
    checksums_path.write_text(existing + f"{digest}  {name}\n")


class TestVerifyChecksums:
    def test_real_matching_checksum_passes(self, tmp_path):
        _write_file_and_checksum(tmp_path, "real_file.txt", b"real content")
        result = verify_checksums(tmp_path)
        assert result == {"real_file.txt": True}

    def test_real_corrupted_file_fails(self, tmp_path):
        _write_file_and_checksum(tmp_path, "real_file.txt", b"real content")
        (tmp_path / "real_file.txt").write_bytes(b"corrupted content")
        result = verify_checksums(tmp_path)
        assert result == {"real_file.txt": False}

    def test_real_multiple_files_all_checked(self, tmp_path):
        _write_file_and_checksum(tmp_path, "file_a.txt", b"aaa")
        _write_file_and_checksum(tmp_path, "file_b.txt", b"bbb")
        result = verify_checksums(tmp_path)
        assert result == {"file_a.txt": True, "file_b.txt": True}


def _build_fake_second_dataset(project_root, n_cells: int = 5, n_clusters: int = 2) -> None:
    dataset_dir = (
        project_root / "data" / "external" / "spatial" / "Xenium_Oliveira_ColorectalCancer_P1"
    )
    clusters_dir = dataset_dir / "analysis" / "clustering" / "gene_expression_graphclust"
    clusters_dir.mkdir(parents=True)

    cells_path = dataset_dir / "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cells.parquet"
    pd.DataFrame(
        {"cell_id": range(n_cells), "x_centroid": range(n_cells), "y_centroid": range(n_cells)}
    ).to_parquet(cells_path)

    clusters_path = clusters_dir / "clusters.csv"
    clusters_path.write_text(
        "Barcode,Cluster\n"
        + "\n".join(f"cell_{i},{(i % n_clusters) + 1}" for i in range(n_cells))
        + "\n"
    )

    for f in [cells_path, clusters_path]:
        digest = hashlib.sha256(f.read_bytes()).hexdigest()
        checksums_path = dataset_dir / "checksums.sha256"
        existing = checksums_path.read_text() if checksums_path.exists() else ""
        checksums_path.write_text(existing + f"{digest}  {f.relative_to(dataset_dir)}\n")


def _skip_real_acquisition(monkeypatch) -> None:
    """These tests exercise the verify-only logic against a synthetic
    fake dataset (or deliberately no dataset at all) -- they must not
    trigger a real network download/extraction of the actual dataset."""
    monkeypatch.setattr(
        "xenium_tcr_ecology.validation.spatial_dataset_acquisition."
        "ensure_second_spatial_dataset_acquired",
        lambda project_root: (
            project_root / "data" / "external" / "spatial" / "Xenium_Oliveira_ColorectalCancer_P1"
        ),
    )


class TestBuildSecondSpatialDatasetAcquisitionSummary:
    def test_real_summary_reports_correct_cell_and_cluster_counts(self, tmp_path, monkeypatch):
        _skip_real_acquisition(monkeypatch)
        _build_fake_second_dataset(tmp_path, n_cells=5, n_clusters=2)
        summary = build_second_spatial_dataset_acquisition_summary(tmp_path)
        assert summary["n_cells"] == 5
        assert summary["n_clusters"] == 2
        assert summary["n_files_verified"] == 2

    def test_real_missing_dataset_directory_raises(self, tmp_path, monkeypatch):
        _skip_real_acquisition(monkeypatch)
        with pytest.raises(PipelineError):
            build_second_spatial_dataset_acquisition_summary(tmp_path)

    def test_real_checksum_mismatch_raises(self, tmp_path, monkeypatch):
        _skip_real_acquisition(monkeypatch)
        _build_fake_second_dataset(tmp_path, n_cells=3, n_clusters=1)
        cells_path = (
            tmp_path
            / "data"
            / "external"
            / "spatial"
            / "Xenium_Oliveira_ColorectalCancer_P1"
            / "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cells.parquet"
        )
        cells_path.write_bytes(b"corrupted")
        with pytest.raises(PipelineError):
            build_second_spatial_dataset_acquisition_summary(tmp_path)
