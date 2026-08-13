"""Standardise staged per-GSM directories into canonical per-section
directories (`02_raw_data_ingestion/05_standardise_sample_directory_layout.py`), symlinked back into data/staged/ rather than copied.

Maps GSM accession -> section_id via metadata/sample_manifest.tsv (Phase
1.01), and renames the canonical role file names (stripping GEO's embedded,
inconsistent timestamp strings) so SpatialData Import onward can read a predictable
layout regardless of original filename.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter
from xenium_tcr_ecology.ingest.xenium_inventory import MANDATORY_ROLES, classify_role

CANONICAL_FILENAMES = {
    "cell_boundaries": "cell_boundaries.parquet.gz",
    "cell_feature_matrix": "cell_feature_matrix.h5",
    "cells": "cells.parquet.gz",
    "morphology": "morphology.ome.tif.gz",
    "nucleus_boundaries": "nucleus_boundaries.parquet.gz",
    "transcripts": "transcripts.parquet.gz",
}

STANDARDIZATION_FIELDS = ["gsm_accession", "section_id", "role", "destination_path", "status"]


def _load_gsm_to_section(sample_manifest_path: Path) -> dict[str, str]:
    if not sample_manifest_path.is_file():
        raise PipelineError(
            f"Missing '{sample_manifest_path}'. Run `01_project_setup_and_governance/01_build_sample_manifest.py` first."
        )
    with sample_manifest_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["gsm_accession"]: row["section_id"] for row in reader}


def standardize_layout(
    staged_root: Path,
    standardised_root: Path,
    sample_manifest_path: Path,
    project_root: Path,
) -> dict:
    if not staged_root.is_dir():
        raise PipelineError(
            f"Staged data directory not found: '{staged_root}'. Run `02_raw_data_ingestion/03_extract_archive_safely.py` first."
        )

    gsm_to_section = _load_gsm_to_section(sample_manifest_path)
    standardised_root.mkdir(parents=True, exist_ok=True)

    inventory_path = (
        project_root
        / "results"
        / "tables"
        / "02_raw_data_ingestion"
        / "standardization_inventory.tsv"
    )
    if inventory_path.exists():
        inventory_path.unlink()
    writer = InventoryWriter(
        inventory_path, project_root=project_root, fields=STANDARDIZATION_FIELDS
    )

    sample_dirs = sorted(p for p in staged_root.iterdir() if p.is_dir())
    unmapped = [d.name for d in sample_dirs if d.name not in gsm_to_section]
    if unmapped:
        raise PipelineError(
            f"{len(unmapped)} staged GSM(s) have no section_id in '{sample_manifest_path}': {unmapped}"
        )

    sections_written = 0
    for sample_dir in sample_dirs:
        gsm = sample_dir.name
        section_id = gsm_to_section[gsm]
        section_dir = standardised_root / section_id
        section_dir.mkdir(parents=True, exist_ok=True)

        files_by_role = {}
        for f in sample_dir.iterdir():
            role = classify_role(f.name)
            if role in MANDATORY_ROLES:
                files_by_role[role] = f

        for role, canonical_name in CANONICAL_FILENAMES.items():
            source = files_by_role.get(role)
            if source is None:
                continue  # already validated as an error by `02_raw_data_ingestion/04_inventory_xenium_files.py`; don't re-fail here
            link_path = section_dir / canonical_name
            if link_path.exists() or link_path.is_symlink():
                link_path.unlink()
            # Relative, not absolute -- an absolute target bakes in this
            # host's repository path and breaks under any bind mount at a
            # different path (Docker's /workspace, Apptainer's Mode 1),
            # even though data/staged/ and data/standardised/ are always
            # at the same fixed relative position to each other.
            relative_target = Path(os.path.relpath(source.resolve(), start=link_path.parent))
            link_path.symlink_to(relative_target)
            writer.write_row(
                gsm_accession=gsm,
                section_id=section_id,
                role=role,
                destination_path=link_path,
                status="linked",
            )
        sections_written += 1

    return {"sections_standardised": sections_written}
