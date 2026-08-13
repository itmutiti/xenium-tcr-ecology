"""Unit tests for xenium_tcr_ecology.ingest.standardize (`02_raw_data_ingestion/05_standardise_sample_directory_layout.py`)."""

from __future__ import annotations

import csv

import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.ingest.standardize import CANONICAL_FILENAMES, standardize_layout

MANDATORY_SUFFIXES = {
    "cell_boundaries": "_cell_boundaries.parquet.gz",
    "cell_feature_matrix": "_cell_feature_matrix.h5",
    "cells": "_cells.parquet.gz",
    "morphology": "_morphology.ome.tif.gz",
    "nucleus_boundaries": "_nucleus_boundaries.parquet.gz",
    "transcripts": "_transcripts.parquet.gz",
}


def _write_manifest(tmp_path, rows: list[dict]):
    path = tmp_path / "sample_manifest.tsv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["gsm_accession", "section_id"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def _make_staged_sample(staged_root, gsm: str):
    d = staged_root / gsm
    d.mkdir(parents=True)
    for role, suffix in MANDATORY_SUFFIXES.items():
        (d / f"{gsm}_weirdname{suffix}").write_bytes(f"{role}-data".encode())
    return d


class TestStandardizeLayout:
    def test_creates_symlinked_canonical_layout(self, tmp_path):
        staged_root = tmp_path / "staged"
        _make_staged_sample(staged_root, "GSM1")
        manifest_path = _write_manifest(
            tmp_path, [{"gsm_accession": "GSM1", "section_id": "P01_run1"}]
        )
        standardised_root = tmp_path / "standardised"

        summary = standardize_layout(
            staged_root, standardised_root, manifest_path, project_root=tmp_path
        )
        assert summary["sections_standardised"] == 1

        section_dir = standardised_root / "P01_run1"
        assert section_dir.is_dir()
        for canonical_name in CANONICAL_FILENAMES.values():
            link = section_dir / canonical_name
            assert link.is_symlink(), f"{canonical_name} should be a symlink, not a copy"
            assert link.read_bytes()  # resolves and reads through the symlink

    def test_raises_on_unmapped_gsm(self, tmp_path):
        staged_root = tmp_path / "staged"
        _make_staged_sample(staged_root, "GSM_UNKNOWN")
        manifest_path = _write_manifest(
            tmp_path, [{"gsm_accession": "GSM_OTHER", "section_id": "P01_run1"}]
        )

        with pytest.raises(PipelineError, match="no section_id"):
            standardize_layout(
                staged_root, tmp_path / "standardised", manifest_path, project_root=tmp_path
            )

    def test_raises_on_missing_manifest(self, tmp_path):
        staged_root = tmp_path / "staged"
        _make_staged_sample(staged_root, "GSM1")
        with pytest.raises(PipelineError, match="Missing"):
            standardize_layout(
                staged_root, tmp_path / "standardised", tmp_path / "nope.tsv", project_root=tmp_path
            )

    def test_raises_on_missing_staged_dir(self, tmp_path):
        manifest_path = _write_manifest(tmp_path, [])
        with pytest.raises(PipelineError, match="not found"):
            standardize_layout(
                tmp_path / "nope", tmp_path / "standardised", manifest_path, project_root=tmp_path
            )
