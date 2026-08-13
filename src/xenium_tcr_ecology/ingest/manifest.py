"""Compile config/geo/sample_manifest_input.yaml into
metadata/sample_manifest.tsv (`01_project_setup_and_governance/01_build_sample_manifest.py`).

Derives patient_id/section_id keys and the is_technical_replicate /
included_in_primary_hnscc_cohort / hpv_positive flags from the raw GEO
fields, rather than requiring the human-curated input to redundantly state
them -- fewer places for the input and the derived manifest to silently
disagree.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter
from xenium_tcr_ecology.infra.validation import validate_records

MANIFEST_FIELDS = [
    "patient_id",
    "section_id",
    "gsm_accession",
    "run_number",
    "specimen_type",
    "included_in_primary_hnscc_cohort",
    "is_technical_replicate",
    "hpv_p16_positive",
    "p16_ihc_status",
    "recurrence_status",
    "smoking_pack_years",
    "tumour_resection_site",
]


def _patient_id(patient_number: int) -> str:
    return f"P{patient_number:02d}"


def _section_id(patient_number: int, run_number: int) -> str:
    return f"{_patient_id(patient_number)}_run{run_number}"


def load_manifest_input(path: Path) -> dict:
    if not path.is_file():
        raise PipelineError(
            f"Required GEO sample-manifest input not found: '{path}'. "
            "This file is manually curated from GSM records, not fabricated -- "
            "see its header comment for provenance and how to (re-)populate it."
        )
    return yaml.safe_load(path.read_text())


def compile_sample_manifest(input_path: Path, output_path: Path, project_root: Path) -> dict:
    """Validate and compile the manifest. Returns a small summary dict
    (counts) so the calling script can sanity-print it and so tests can
    assert on cohort structure without re-parsing the TSV."""
    data = load_manifest_input(input_path)
    samples = data["samples"]

    validate_records(
        samples,
        required_fields=[
            "gsm_accession",
            "patient_number",
            "run_number",
            "specimen_type",
            "p16_ihc_status",
        ],
        source_path=input_path,
        min_records=1,
    )

    # Recomputes uniqueness directly from the raw per-sample records rather
    # than trusting any upstream summary, so a genuine discrepancy would be
    # caught here regardless.
    gsm_ids = [s["gsm_accession"] for s in samples]
    if len(gsm_ids) != len(set(gsm_ids)):
        dupes = [g for g, n in Counter(gsm_ids).items() if n > 1]
        raise PipelineError(f"Duplicate GSM accession(s) in '{input_path}': {dupes}")

    if output_path.exists():
        output_path.unlink()
    writer = InventoryWriter(
        output_path, project_root=project_root, fields=MANIFEST_FIELDS, delimiter="\t"
    )

    patient_run_counts: Counter = Counter()
    for s in samples:
        patient_run_counts[s["patient_number"]] += 1

    hnscc_patients = set()
    replicated_patients = set()
    hnscc_sections = 0
    ameloblastoma_count = 0
    hpv_positive_patients = set()

    for s in samples:
        pnum = s["patient_number"]
        is_hnscc = s["specimen_type"] == "HNSCC tumor"
        is_replicate = patient_run_counts[pnum] > 1
        is_hpv_positive = s["p16_ihc_status"] == "Positive"

        if is_hnscc:
            hnscc_patients.add(pnum)
            hnscc_sections += 1
            if is_replicate:
                replicated_patients.add(pnum)
            if is_hpv_positive:
                hpv_positive_patients.add(pnum)
        else:
            ameloblastoma_count += 1

        writer.write_row(
            patient_id=_patient_id(pnum),
            section_id=_section_id(pnum, s["run_number"]),
            gsm_accession=s["gsm_accession"],
            run_number=s["run_number"],
            specimen_type=s["specimen_type"],
            included_in_primary_hnscc_cohort=is_hnscc,
            is_technical_replicate=is_replicate,
            hpv_p16_positive=is_hpv_positive,
            p16_ihc_status=s["p16_ihc_status"],
            recurrence_status=s["recurrence_status"],
            smoking_pack_years=s["smoking_pack_years"],
            tumour_resection_site=s["tumour_resection_site"],
        )

    return {
        "total_samples": len(samples),
        "hnscc_patients": len(hnscc_patients),
        "hnscc_sections": hnscc_sections,
        "replicated_patients": len(replicated_patients),
        "hpv_positive_patients": len(hpv_positive_patients),
        "ameloblastoma_specimens": ameloblastoma_count,
    }
