"""Downloads (if not already present) and verifies the two independent
Xenium spatial datasets used for external validation
(`16_external_validation_and_generalisation/01_acquire_independent_spatial_dataset.py`,
`16_external_validation_and_generalisation/08_acquire_second_independent_spatial_dataset.py`):
Janesick et al. 2023 (breast cancer) and de Oliveira et al. 2025
(colorectal cancer), both public 10x Genomics datasets, CC BY 4.0, no
access restriction. See `data/external/spatial/
Xenium_Janesick_BreastCancer_Rep1/README.md` and `data/external/spatial/
Xenium_Oliveira_ColorectalCancer_P1/README.md` for the full acquisition
provenance, citations and license (the same READMEs this module's exact
URLs, filenames and file sizes are taken from) -- including the
tar-vs-tar.gz quirk handled below: the breast-cancer dataset's
`analysis` bundle is served at a `.tar.gz`-suffixed URL but its content
is plain (unzipped) POSIX tar (verified with `file` during
acquisition); the colorectal dataset's own `.tar.gz` is a
real gzipped tar. Both are extracted with Python's standard-library
`tarfile`, which auto-detects gzip compression regardless of filename.
"""

from __future__ import annotations

import tarfile
from pathlib import Path
from typing import Mapping

import pandas as pd

from xenium_tcr_ecology.infra.download import (
    download_file,
    verify_checksums,
)  # noqa: F401 -- re-exported
from xenium_tcr_ecology.infra.exceptions import PipelineError

DATASET_DIR_NAME = "Xenium_Janesick_BreastCancer_Rep1"
SECOND_DATASET_DIR_NAME = "Xenium_Oliveira_ColorectalCancer_P1"

_JANESICK_BASE_URL = (
    "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1"
)
_JANESICK_FILES = {
    # (dest filename, remote filename if different, expected size in bytes)
    "Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.h5": (None, 12_148_885),
    "Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.parquet": (None, 3_453_894),
    # Served at a .tar.gz-suffixed URL but the content is plain tar
    # (verified with `file`) -- saved locally without the
    # .gz suffix to match its real content, per checksums.sha256.
    "Xenium_FFPE_Human_Breast_Cancer_Rep1_analysis.tar": (
        "Xenium_FFPE_Human_Breast_Cancer_Rep1_analysis.tar.gz",
        64_317_440,
    ),
}

_OLIVEIRA_BASE_URL = (
    "https://cf.10xgenomics.com/samples/xenium/2.0.0/"
    "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE"
)
_OLIVEIRA_FILES = {
    "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cell_feature_matrix.h5": (None, 14_168_777),
    "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cells.parquet": (None, 5_503_158),
    "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_analysis.tar.gz": (None, 51_947_742),
}


def _ensure_dataset_acquired(
    dataset_dir: Path,
    base_url: str,
    files: Mapping[str, tuple[str | None, int]],
    analysis_tar_name: str,
) -> None:
    dataset_dir.mkdir(parents=True, exist_ok=True)
    for dest_name, (remote_name, expected_size) in files.items():
        url = f"{base_url}/{remote_name or dest_name}"
        download_file(url, dataset_dir / dest_name, expected_size_bytes=expected_size)

    analysis_dir = dataset_dir / "analysis"
    if not analysis_dir.is_dir():
        with tarfile.open(dataset_dir / analysis_tar_name) as tf:
            tf.extractall(dataset_dir, filter="data")


def ensure_spatial_dataset_acquired(project_root: Path) -> Path:
    dataset_dir = project_root / "data" / "external" / "spatial" / DATASET_DIR_NAME
    _ensure_dataset_acquired(
        dataset_dir,
        _JANESICK_BASE_URL,
        _JANESICK_FILES,
        "Xenium_FFPE_Human_Breast_Cancer_Rep1_analysis.tar",
    )
    return dataset_dir


def ensure_second_spatial_dataset_acquired(project_root: Path) -> Path:
    dataset_dir = project_root / "data" / "external" / "spatial" / SECOND_DATASET_DIR_NAME
    _ensure_dataset_acquired(
        dataset_dir,
        _OLIVEIRA_BASE_URL,
        _OLIVEIRA_FILES,
        "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_analysis.tar.gz",
    )
    return dataset_dir


def build_spatial_dataset_acquisition_summary(project_root: Path) -> dict:
    dataset_dir = ensure_spatial_dataset_acquired(project_root)
    cells_path = dataset_dir / "Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.parquet"
    return _verify_spatial_dataset_acquisition(dataset_dir, cells_path)


def build_second_spatial_dataset_acquisition_summary(project_root: Path) -> dict:
    """Acquisition verification for the second, independent Xenium
    dataset (colorectal cancer, de Oliveira et al. 2025 -- `data/
    external/spatial/Xenium_Oliveira_ColorectalCancer_P1/README.md`),
    identical in role to the breast-cancer verification above."""
    dataset_dir = ensure_second_spatial_dataset_acquired(project_root)
    cells_path = dataset_dir / "Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_cells.parquet"
    return _verify_spatial_dataset_acquisition(dataset_dir, cells_path)


def _verify_spatial_dataset_acquisition(dataset_dir: Path, cells_path: Path) -> dict:
    clusters_path = (
        dataset_dir / "analysis" / "clustering" / "gene_expression_graphclust" / "clusters.csv"
    )

    if not dataset_dir.is_dir():
        raise PipelineError(f"'{dataset_dir}' not found -- dataset has not been acquired.")
    for path in [cells_path, clusters_path]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found -- dataset acquisition is incomplete.")

    checksum_results = verify_checksums(dataset_dir)
    failed = [f for f, ok in checksum_results.items() if not ok]
    if failed:
        raise PipelineError(
            f"Checksum verification failed for {failed} in '{dataset_dir}' -- re-download required, not silently trusted."
        )

    cells = pd.read_parquet(cells_path)
    clusters = pd.read_csv(clusters_path)

    return {
        "dataset_dir": str(dataset_dir),
        "n_files_verified": len(checksum_results),
        "n_cells": len(cells),
        "n_clusters": int(clusters["Cluster"].nunique()),
    }
