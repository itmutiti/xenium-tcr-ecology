"""Detect Xenium output file roles per sample and flag missing mandatory
assets (`02_raw_data_ingestion/04_inventory_xenium_files.py`).

Mandatory/optional split verified against GSE300147's
filelist.txt (not from generic Xenium documentation): every one of
the 18 samples has exactly 6 mandatory files (cell_boundaries.parquet.gz,
cell_feature_matrix.h5, cells.parquet.gz, morphology.ome.tif.gz,
nucleus_boundaries.parquet.gz, transcripts.parquet.gz); 14 of 18 additionally
have a 7th, differently-named standalone raw microscopy .tif.gz that is
optional and inconsistently present.
"""

from __future__ import annotations

from pathlib import Path

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter

MANDATORY_ROLES = {
    "cell_boundaries": "_cell_boundaries.parquet.gz",
    "cell_feature_matrix": "_cell_feature_matrix.h5",
    "cells": "_cells.parquet.gz",
    "morphology": "_morphology.ome.tif.gz",
    "nucleus_boundaries": "_nucleus_boundaries.parquet.gz",
    "transcripts": "_transcripts.parquet.gz",
}

INVENTORY_FIELDS = ["gsm_accession", "role", "filename", "present", "size_bytes"]


def classify_role(filename: str) -> str | None:
    for role, suffix in MANDATORY_ROLES.items():
        if filename.endswith(suffix):
            return role
    if filename.endswith(".tif.gz") and not filename.endswith("_morphology.ome.tif.gz"):
        return "raw_microscopy_optional"
    return None


def build_xenium_file_inventory(staged_root: Path, output_path: Path, project_root: Path) -> dict:
    if not staged_root.is_dir():
        raise PipelineError(
            f"Staged data directory not found: '{staged_root}'. Run `02_raw_data_ingestion/03_extract_archive_safely.py` first."
        )

    sample_dirs = sorted(p for p in staged_root.iterdir() if p.is_dir())
    if not sample_dirs:
        raise PipelineError(f"No sample directories found under '{staged_root}'.")

    if output_path.exists():
        output_path.unlink()
    writer = InventoryWriter(output_path, project_root=project_root, fields=INVENTORY_FIELDS)

    incomplete_samples = []
    for sample_dir in sample_dirs:
        gsm = sample_dir.name
        files_by_role: dict[str, Path] = {}
        for f in sample_dir.iterdir():
            role = classify_role(f.name)
            if role is not None:
                files_by_role.setdefault(role, f)

        missing = [role for role in MANDATORY_ROLES if role not in files_by_role]
        if missing:
            incomplete_samples.append((gsm, missing))

        for role in MANDATORY_ROLES:
            f = files_by_role.get(role)
            writer.write_row(
                gsm_accession=gsm,
                role=role,
                filename=f.name if f else "",
                present=f is not None,
                size_bytes=f.stat().st_size if f else "",
            )
        if "raw_microscopy_optional" in files_by_role:
            f = files_by_role["raw_microscopy_optional"]
            writer.write_row(
                gsm_accession=gsm,
                role="raw_microscopy_optional",
                filename=f.name,
                present=True,
                size_bytes=f.stat().st_size,
            )

    if incomplete_samples:
        raise PipelineError(
            f"{len(incomplete_samples)} sample(s) missing mandatory files: {incomplete_samples}"
        )

    return {
        "samples_scanned": len(sample_dirs),
        "samples_complete": len(sample_dirs) - len(incomplete_samples),
    }
