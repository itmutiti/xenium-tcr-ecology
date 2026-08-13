"""Attach clinical/technical metadata from metadata/sample_manifest.tsv onto
each section's AnnData (`03_spatialdata_import/04_attach_clinical_and_technical_metadata.py`), validating the join key is unique on
both sides before merging -- silently allowing a many-to-many join here
would duplicate cells or corrupt every downstream analysis."""

from __future__ import annotations

import csv
from pathlib import Path

import anndata as ad

from xenium_tcr_ecology.infra.exceptions import PipelineError

METADATA_FIELDS = [
    "patient_id",
    "gsm_accession",
    "run_number",
    "included_in_primary_hnscc_cohort",
    "is_technical_replicate",
    "hpv_p16_positive",
    "p16_ihc_status",
    "recurrence_status",
    "smoking_pack_years",
    "tumour_resection_site",
]

# csv.DictReader (stdlib) reads every TSV field as a raw string -- with no
# cast, boolean-like fields end up stored in obs as the literal strings
# "True"/"False" (pandas then broadcasts this to a `category` dtype), not
# a native bool column. This was not caught until `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` needed to use
# `is_technical_replicate` as an actual Python boolean mask -- silently
# present in every downstream artifact since `03_spatialdata_import/04_attach_clinical_and_technical_metadata.py` first ran. The R-side
# equivalent of this same hazard (arrow/read.delim reading these same
# fields back out of parquet/TSV) was already found and fixed in `04_quality_control/08_assess_replicate_concordance.R` --
# this is the Python-side counterpart of the same root cause.
BOOLEAN_METADATA_FIELDS = [
    "included_in_primary_hnscc_cohort",
    "is_technical_replicate",
    "hpv_p16_positive",
]


def _load_manifest(sample_manifest_path: Path) -> dict[str, dict]:
    if not sample_manifest_path.is_file():
        raise PipelineError(
            f"Missing '{sample_manifest_path}'. Run `01_project_setup_and_governance/01_build_sample_manifest.py` first."
        )
    with sample_manifest_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    section_ids = [r["section_id"] for r in rows]
    if len(section_ids) != len(set(section_ids)):
        raise PipelineError(
            f"'{sample_manifest_path}' has duplicate section_id values -- join key must be unique."
        )

    return {r["section_id"]: r for r in rows}


def attach_metadata(anndata_path: Path, sample_manifest_path: Path, output_path: Path) -> dict:
    if not anndata_path.is_file():
        raise PipelineError(
            f"'{anndata_path}' not found. Run `03_spatialdata_import/03_create_anndata_expression_objects.py` first."
        )

    adata = ad.read_h5ad(anndata_path)
    if "section_id" not in adata.obs.columns:
        raise PipelineError(f"'{anndata_path}': obs has no 'section_id' column to join on.")

    section_ids_in_data = adata.obs["section_id"].unique()
    if len(section_ids_in_data) != 1:
        raise PipelineError(
            f"'{anndata_path}' contains more than one section_id: {section_ids_in_data}."
        )
    section_id = section_ids_in_data[0]

    manifest = _load_manifest(sample_manifest_path)
    if section_id not in manifest:
        raise PipelineError(f"section_id '{section_id}' not found in '{sample_manifest_path}'.")

    row = manifest[section_id]
    for field in METADATA_FIELDS:
        value = row[field]
        if field in BOOLEAN_METADATA_FIELDS:
            value = value == "True"
        adata.obs[field] = value

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    return {"section_id": section_id, "n_cells": adata.n_obs, "patient_id": row["patient_id"]}
