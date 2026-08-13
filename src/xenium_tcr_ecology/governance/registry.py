"""Compile config/governance/analysis_registry_input.yaml into
governance/analysis_registry.tsv (`01_project_setup_and_governance/06_create_analysis_registry.py`), enforcing the HPV
single-contrast cap while doing so."""

from __future__ import annotations

from pathlib import Path

import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter
from xenium_tcr_ecology.infra.validation import validate_records

REGISTRY_FIELDS = [
    "analysis_id",
    "phase",
    "script",
    "hypothesis",
    "unit_of_analysis",
    "exclusion_criteria",
    "primary_endpoint",
    "multiplicity_family",
    "taxonomy_version",
    "status",
    "registered_by",
    "registered_date",
]

MAX_HPV_PRIMARY_CONTRASTS = 2


def load_registry_input(path: Path) -> dict:
    if not path.is_file():
        raise PipelineError(f"Required analysis-registry input not found: '{path}'.")
    return yaml.safe_load(path.read_text())


def compile_analysis_registry(
    input_path: Path,
    output_path: Path,
    project_root: Path,
    registered_by: str,
    registered_date: str,
) -> dict:
    data = load_registry_input(input_path)
    analyses = data["analyses"]

    validate_records(
        analyses,
        required_fields=[
            "analysis_id",
            "phase",
            "hypothesis",
            "unit_of_analysis",
            "primary_endpoint",
            "multiplicity_family",
        ],
        source_path=input_path,
        min_records=1,
    )

    hpv_entries = [
        a for a in analyses if a["multiplicity_family"].lower().startswith("primary (hpv)")
    ]
    if len(hpv_entries) > MAX_HPV_PRIMARY_CONTRASTS:
        raise PipelineError(
            f"{len(hpv_entries)} entries registered under multiplicity_family 'primary (HPV)...', "
            f"exceeding the cap of {MAX_HPV_PRIMARY_CONTRASTS} (config/reproducibility_policy.yaml "
            "gates.hpv_single_contrast)."
        )

    if output_path.exists():
        output_path.unlink()
    writer = InventoryWriter(output_path, project_root=project_root, fields=REGISTRY_FIELDS)
    for a in analyses:
        writer.write_row(
            analysis_id=a["analysis_id"],
            phase=a["phase"],
            script=a.get("script", ""),
            hypothesis=a["hypothesis"],
            unit_of_analysis=a["unit_of_analysis"],
            exclusion_criteria=a.get("exclusion_criteria", ""),
            primary_endpoint=a["primary_endpoint"],
            multiplicity_family=a["multiplicity_family"],
            taxonomy_version=a.get("taxonomy_version", "n/a"),
            status=a.get("status", "planned"),
            registered_by=registered_by,
            registered_date=registered_date,
        )

    primary_families = {
        a["multiplicity_family"]
        for a in analyses
        if a["multiplicity_family"].lower().startswith("primary")
    }
    return {
        "analyses_registered": len(analyses),
        "hpv_primary_contrasts_reserved": len(hpv_entries),
        "distinct_primary_families": len(primary_families),
    }
